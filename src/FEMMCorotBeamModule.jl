module FEMMCorotBeamModule

using LinearAlgebra: norm, Transpose, mul!
using FinEtools
using FinEtools.MatrixUtilityModule: complete_lt!
using FinEtools.IntegDomainModule: IntegDomain
import FinEtoolsDeforLinear.MatDeforElastIsoModule: MatDeforElastIso
using ..FESetL2BeamModule: FESetL2Beam, initial_local_frame!


"""
    FEMMCorotBeam{S<:FESetL2Beam, F<:Function} <: AbstractFEMM

Class for co-rotational beam finite element modeling machine.
"""
mutable struct FEMMCorotBeam{S<:FESetL2Beam, F<:Function} <: AbstractFEMM
    integdomain::IntegDomain{S, F} # integration domain data
    material::MatDeforElastIso # material object
    # The attributes below are buffers used in various operations.
    _ecoords0::FFltMat
    _ecoords1::FFltMat
    _edisp1::FFltMat
    _evel1::FFltMat
    _evel1f::FFltMat
    _dofnums::FIntMat
    _F0::FFltMat
    _Ft::FFltMat
    _FtI::FFltMat
    _FtJ::FFltMat
    _Te::FFltMat
    _tempelmat1::FFltMat
    _tempelmat2::FFltMat
    _tempelmat3::FFltMat
    _elmat::FFltMat
    _elmatTe::FFltMat
    _elmato::FFltMat
    _elvec::FFltVec
    _elvecf::FFltVec
    _aN::FFltMat
    _dN::FFltVec
    _DN::FFltMat
    _PN::FFltVec
    _LF::FFltVec
    _RI::FFltMat
    _RJ::FFltMat
    _OS::FFltMat
end

function FEMMCorotBeam(integdomain::IntegDomain{S, F}, material::MatDeforElastIso) where {S<:FESetL2Beam, F<:Function}
    _ecoords0 = fill(0.0, 2, 3); 
    _ecoords1 = fill(0.0, 2, 3)
    _edisp1 = fill(0.0, 2, 3); 
    _evel1 = fill(0.0, 2, 6); 
    _evel1f = fill(0.0, 2, 6)
    _dofnums = zeros(FInt, 1, 12); 
    _F0 = fill(0.0, 3, 3); 
    _Ft = fill(0.0, 3, 3); 
    _FtI = fill(0.0, 3, 3); 
    _FtJ = fill(0.0, 3, 3)
    _Te = fill(0.0, 12, 12)
    _tempelmat1 = fill(0.0, 12, 12); 
    _tempelmat2 = fill(0.0, 12, 12); 
    _tempelmat3 = fill(0.0, 12, 12)
    _elmat = fill(0.0, 12, 12);    
    _elmatTe = fill(0.0, 12, 12);    
    _elmato = fill(0.0, 12, 12)
    _elvec = fill(0.0, 12);    
    _elvecf = fill(0.0, 12)
    _aN = fill(0.0, 6, 12)
    _dN = fill(0.0, 6)
    _DN = fill(0.0, 6, 6)
    _PN = fill(0.0, 6)
    _LF = fill(0.0, 12)
    _RI = fill(0.0, 3, 3);    
    _RJ = fill(0.0, 3, 3);    
    _OS = fill(0.0, 3, 3)
    return FEMMCorotBeam(integdomain, material,
     _ecoords0, _ecoords1, _edisp1, _evel1, _evel1f, 
     _dofnums, 
     _F0, _Ft, _FtI, _FtJ, _Te,
     _tempelmat1, _tempelmat2, _tempelmat3, _elmat, _elmatTe, _elmato, 
     _elvec, _elvecf, 
     _aN, _dN, _DN, _PN, _LF, 
     _RI, _RJ, _OS)
end

function _transfmat!(Te, Ft)
    Te[1:3, 1:3] = Te[4:6, 4:6] = Te[7:9, 7:9] = Te[10:12, 10:12] = Ft
    return Te
end


"""
    Type of the mass matrix formulation.
"""
const MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA=0;
const MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA=1;
const MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA=2;
const MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA=3;

"""
    current_local_frame!(Ft, xt, FtI, FtJ)

Compute the current local frame from the positions of the nodes and the
orientation matrices (frames) at the nodes.
"""
function current_local_frame!(Ft, xt, FtI, FtJ)
    Ft = fill(0.0, 3, 3);
    # Compute the element frame in configuration t
    Ft[:,1] = (xt[2,:]-xt[1,:]);
    Lt = norm(@view Ft[:,1]);
    Ft[:,1] /= Lt;
    x1x2_vectort = FtI[:,2]+FtJ[:,2];
    #     Ft(:,3)=skewmat(Ft(:,1))*x1x2_vectort; Ft(:,3)=Ft(:,3)/norm(Ft(:,3));
    #     Ft(:,2)=skewmat(Ft(:,3))*Ft(:,1);
    #     Ft(:,3)=skewmat(Ft(:,1))*x1x2_vectort; # In the interest of speed,
    #     replace with below explicit rewrite
    Ft[:,3] .= (-Ft[3,1]*x1x2_vectort[2]+Ft[2,1]*x1x2_vectort[3],
                Ft[3,1]*x1x2_vectort[1]-Ft[1,1]*x1x2_vectort[3],
               -Ft[2,1]*x1x2_vectort[1]+Ft[1,1]*x1x2_vectort[2]);
    Ft[:,3] /= norm(@view Ft[:,3]);
    #     Ft(:,2)=skewmat(Ft(:,3))*Ft(:,1); # In the interest of speed,
    #     replace with below explicit rewrite
    Ft[:,2] .= (-Ft[3,3]*Ft[2,1]+Ft[2,3]*Ft[3,1],
                Ft[3,3]*Ft[1,1]-Ft[1,3]*Ft[3,1],
               -Ft[2,3]*Ft[1,1]+Ft[1,3]*Ft[2,1]);
    return Lt, Ft
end


"""
    local_frame_and_def!(Ft, dN, F0, FtI, FtJ, x0, x1x2_vector, xt, RI, RJ)

Compute the current length of the element, the current element frame,
and the natural deformations

# Arguments
x0= array of node coordinates, one node per row, in initial configuration
x1x2_vector= vector that lies in the x1-x2 local element coordinate
     plane, in initial configuration
xt= array of node coordinates, one node per row, in initial configuration
RI,RJ=nodal rotation (orthogonal) matrix

# Outputs
`Lt`= current length of the element,
`Ft`= current element frame (orthogonal rotation matrix) whose columns are unit
     vectors: they are centered halfway between the current locations of
     the nodes, vector 1 points from node I to node J, vector 3 is
     orthogonal to the sum of the nodal cross-section frame vectors 2
`dN`= vector of natural deformations; dN(1)= total change in length
     between configurations 0 and t; dN(2)= symmetric bending;
     dN(3)= anti-symmetric bending;    dN(4)= symmetric bending
     dN(5)= anti-symmetric bending; dN(6)=total axial torsion angle.

"""
function local_frame_and_def!(Ft, dN, F0, FtI, FtJ, x0, x1x2_vector, xt, RI, RJ)
    # This is the element frame in the configuration t=0
    L0, F0 = initial_local_frame!(F0, x0, x1x2_vector)

    # The nodal cross-section frames are rotated by the nodal rotation matrices
    mul!(FtI, RI, F0);
    mul!(FtJ, RJ, F0);

    # The element frame in the configuration t
    Lt, Ft = current_local_frame!(Ft, xt, FtI, FtJ)

    # Components of FtI,FtJ in the element frame
    LFtI = Ft'*FtI;
    LFtJ = Ft'*FtJ;

    # C  COMPUTE NET DEFORMATIONS- DELTA L,TH1IJ,TH2I,TH2J,TH3I,TH3J
    dN = fill(0.0, 6);
    # Total change in length between configurations 0 and t
    dN[1] = Lt-L0;
    # Total axial TORSION ANGLE
    dN[6] = (LFtJ[3,2]/LFtJ[2,2]-LFtI[3,2]/LFtI[2,2] -LFtJ[2,3]/LFtJ[3,3]+LFtI[2,3]/LFtI[3,3])/2;
    # Total rotation  angles of nodal cross-sections relative to the element frame
    TH2I =-LFtI[3,1]/LFtI[1,1];
    TH2J =-LFtJ[3,1]/LFtJ[1,1];
    TH3I = LFtI[2,1]/LFtI[1,1];
    TH3J = LFtJ[2,1]/LFtJ[1,1];
    dN[2] = TH3I-TH3J; # symmetric bending
    dN[3] = TH3I+TH3J; # anti-symmetric bending
    dN[4] = -TH2I+TH2J; # symmetric bending
    dN[5] = -TH2I-TH2J; # anti-symmetric bending
    return Lt, Ft, dN
end


"""
    local_cartesian_to_natural!(aN, L)

Compute transformation from local Cartesian displacements to natural deformations
of a beam element.

Matrix defined in Eq (4.8) of COMPUTER  METHODS  IN APPLIED  MECHANICS  AND  ENGINEERING   14  (1978)  401-451
ON LARGE  DISPLACEMENT-SMALL   STRAIN  ANALYSIS   OF  STRUCTURES
WITH  ROTATIONAL   DEGREES  OF  FREEDOM.
J.H.  ARGYRIS,   P.C.  DUNNE  and  D.W. SCHARPF

# Arguments
`L`= current length of the element

# Outputs
`aN`= transformation matrix to take Cartesian (local) displacement increments in the
     element frame and to produce increments of natural deformations;
     see `local_frame_and_def!` for the definition of the natural deformations
"""
function local_cartesian_to_natural!(aN, L)
    fill!(aN, 0.0)
    aN[1, 1] = -1; aN[1, 7] = +1
    aN[2, 6] = +1; aN[2, 12] = -1
    aN[3, 2] = 2/L; aN[3, 6] = +1; aN[3, 8] = -2/L; aN[3, 12] = +1
    aN[4, 5] = -1; aN[4, 11] = +1
    aN[5, 3] = 2/L; aN[5, 5] = -1; aN[5, 9] = -2/L; aN[5, 11] = -1
    aN[6, 4] = -1; aN[6, 10] = +1
    # aN=[[ -1,   0,   0,  0,  0, 0, 1,    0,    0, 0,  0,  0]
    #     [  0,   0,   0,  0,  0, 1, 0,    0,    0, 0,  0, -1]
    #     [  0, 2/L,   0,  0,  0, 1, 0, -2/L,    0, 0,  0,  1]
    #     [  0,   0,   0,  0, -1, 0, 0,    0,    0, 0,  1,  0]
    #     [  0,   0, 2/L,  0, -1, 0, 0,    0, -2/L, 0, -1,  0]
    #     [  0,   0,   0, -1,  0, 0, 0,    0,    0, 1,  0,  0]];
        return aN
end


"""
    local_geometric_stiffness!(SM, A, I2, I3, PN, L)

Compute the local geometric stiffness matrix.

# Arguments
`A`= cross-sectional area,
`I2`, `I3`=central moment of inertia of the cross-section about the x2 and x3
coordinate axis,
`PN`= vector of natural forces; see natural_forces() for definitions
`L`= current length of the element,

# Outputs
`SM` = local geometric stiffness matrix, 12 x 12

This geometric stiffness matrix this consistent with relationship between
the natural deformations and the natural forces that assumes there is
only a linear constitutive link: no non-constitutive effects (bowing
etc.) are included. This form of the geometric matrix was derived by
Krenk.
@BOOK{Krenk:2009,
  AUTHOR =       {S. Krenk},
  TITLE =        {Non-linear Modeling and Analysis of Solids and Structures },
  PUBLISHER =    {Cambridge University Press},
  YEAR =         {2009},
  isbn =         {9780521830546}
}
"""
function local_geometric_stiffness!(SM, A, I2, I3, PN, L)
    N = PN[1];
    S_2 = -2*PN[3]/L;
    S_3 = -2*PN[5]/L;
    M_1 = PN[6];
    M_2I = PN[4]+PN[5];
    M_2J = PN[4]-PN[5];
    M_3I = -(PN[2]+PN[3]);
    M_3J = -(PN[2]-PN[3]);
    SM[1:3, 1:3] = [0   -S_2/L   -S_3/L;
                    -S_2/L   N/L   0;
                    -S_3/L   0   N/L];
    SM[1:3, 4:6] = [0   0   0;
                    -M_2I/L   M_1/L   0;
                    -M_3I/L   0   M_1/L];
    SM[1:3, 7:9] = [0   +S_2/L   +S_3/L;
                    +S_2/L   -N/L   0;
                    +S_3/L   0   -N/L];
    SM[1:3, 10:12] = [0   0   0;
                    M_2J/L   -M_1/L   0;
                    M_3J/L   0   -M_1/L];
    SM[4:6, 4:6] = [0   M_3I/2   -M_2I/2;
                    M_3I/2   0   0;
                    -M_2I/2   0   0];
    SM[4:6, 7:9] = [0   M_2I/L   M_3I/L;
                    0   -M_1/L   0;
                    0   0   -M_1/L];
    SM[4:6, 10:12] = [0   0   0;
                    0   0   M_1/2;
                    0   -M_1/2   0];
    SM[7:9, 7:9] = [0   -S_2/L   -S_3/L;
                    -S_2/L   N/L   0;
                    -S_3/L   0   N/L];
    SM[7:9, 10:12] = [0   0   0;
                    -M_2J/L   M_1/L   0;
                    -M_3J/L   0   M_1/L];
    SM[10:12, 10:12] = [0   -M_3J/2   M_2J/2;
                        -M_3J/2   0   0;
                        M_2J/2   0   0];
    complete_lt!(SM)
    return SM
end

function local_mass_CONSISTENT_NO_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    # C
    # C  CONSISTENT MASS MATRIX including ROTATIONAL MASSES
    # C  Formulation of the (3.38), (3.39) equation from Dykstra's thesis
    fill!(MM, 0.0)
    c1 = (rho*A*L)::Float64
    MM[1, 1] = c1 * (1/3)
    MM[1, 7] = c1 * (1/6)
    MM[2, 2] = c1 * (13/35)
    MM[2, 6] = c1 * (11*L/210)
    MM[2, 8] = c1 * (9/70)
    MM[2, 12] = c1 * (-13*L/420)
    MM[3, 3] = c1 * (13/35)
    MM[3, 5] = c1 * (-11*L/210)
    MM[3, 9] = c1 * (9/70)
    MM[3, 11] = c1 * (13*L/420)
    MM[4, 4] = c1 * (I1/3/A)
    MM[4, 10] = c1 * (I1/6/A)
    MM[5, 5] = c1 * (L^2/105)
    MM[5, 9] = c1 * (-13*L/420)
    MM[5, 11] = c1 * (-L^2/140)
    MM[6, 6] = c1 * (L^2/105)
    MM[6, 8] = c1 * (13*L/420)
    MM[6, 12] = c1 * (-L^2/140)
    MM[7, 7] = c1 * (1/3)
    MM[8, 8] = c1 * (13/35)
    MM[8, 12] = c1 * (-11*L/210)
    MM[9, 9] = c1 * (13/35)
    MM[9, 11] = c1 * (11*L/210)
    MM[10, 10] = c1 * (I1/3/A)
    MM[11, 11] = c1 * (L^2/105)
    MM[12, 12] = c1 * (L^2/105)
    # 1       2        3        4       5         6       7        8       9       10      11       12
    complete_lt!(MM)
    return MM
end

function local_mass_CONSISTENT_WITH_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    # C
    # C  CONSISTENT MASS MATRIX including ROTATIONAL MASSES
    # C  Formulation of the (3.38), (3.39) equation from Dykstra's thesis
    fill!(MM, 0.0)
    c1 = (rho*A*L)::Float64
    MM[1, 1] = c1 * (1/3)
    MM[1, 7] = c1 * (1/6)
    MM[2, 2] = c1 * (13/35)
    MM[2, 6] = c1 * (11*L/210)
    MM[2, 8] = c1 * (9/70)
    MM[2, 12] = c1 * (-13*L/420)
    MM[3, 3] = c1 * (13/35)
    MM[3, 5] = c1 * (-11*L/210)
    MM[3, 9] = c1 * (9/70)
    MM[3, 11] = c1 * (13*L/420)
    MM[4, 4] = c1 * (I1/3/A)
    MM[4, 10] = c1 * (I1/6/A)
    MM[5, 5] = c1 * (L^2/105)
    MM[5, 9] = c1 * (-13*L/420)
    MM[5, 11] = c1 * (-L^2/140)
    MM[6, 6] = c1 * (L^2/105)
    MM[6, 8] = c1 * (13*L/420)
    MM[6, 12] = c1 * (-L^2/140)
    MM[7, 7] = c1 * (1/3)
    MM[8, 8] = c1 * (13/35)
    MM[8, 12] = c1 * (-11*L/210)
    MM[9, 9] = c1 * (13/35)
    MM[9, 11] = c1 * (11*L/210)
    MM[10, 10] = c1 * (I1/3/A)
    MM[11, 11] = c1 * (L^2/105)
    MM[12, 12] = c1 * (L^2/105)
    c2 = (rho/L)::Float64
    MM[2, 2] += c2 * (6/5*I2)
    MM[2, 6] += c2 * (L/10*I2)
    MM[2, 8] += c2 * (-6/5*I2)
    MM[2, 12] += c2 * (L/10*I2)
    MM[3, 3] += c2 * (6/5*I3)
    MM[3, 5] += c2 * (-L/10*I3)
    MM[3, 9] += c2 * (-6/5*I3)
    MM[3, 11] += c2 * (-L/10*I3)
    MM[5, 5] += c2 * (2*L^2/15*I3)
    MM[5, 9] += c2 * (L/10*I3)
    MM[5, 11] += c2 * (-L^2/30*I3)
    MM[6, 6] += c2 * (2*L^2/15*I2)
    MM[6, 8] += c2 * (-L/10*I2)
    MM[6, 12] += c2 * (-L^2/30*I2)
    MM[8, 8] += c2 * (6/5*I2)
    MM[8, 12] += c2 * (-L/10*I2)
    MM[9, 9] += c2 * (6/5*I3)
    MM[9, 11] += c2 * (L/10*I3)
    MM[11, 11] += c2 * (2*L^2/15*I3)
    MM[12, 12] += c2 * (2*L^2/15*I2)
    # 1       2        3        4       5         6       7        8       9       10      11       12
    complete_lt!(MM)
    return MM
end

function local_mass_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    # C
    # C  LUMPED DIAGONAL MASS MATRIX WITH ROTATIONAL MASSES
    # C
    HLM  = A*rho*L/2.;
    HLI1 = rho*I1* L/2.;
    HLI2 = rho*I2* L/2.;
    HLI3 = rho*I3* L/2.;
    CA = HLM;
    CB = HLI1;
    CC = HLI2;
    CD = HLI3;
    fill!(MM, 0.0);
    MM[1,1]    = MM[1,1]    + CA;
    MM[2,2]    = MM[2,2]    + CA;
    MM[3,3]    = MM[3,3]    + CA;
    MM[4,4]    = MM[4,4]    + CB;
    MM[5,5]    = MM[5,5]    + CC;
    MM[6,6]    = MM[6,6]    + CD;
    MM[7,7]    = MM[7,7]    + CA;
    MM[8,8]    = MM[8,8]    + CA;
    MM[9,9]    = MM[9,9]    + CA;
    MM[10,10]  = MM[10,10]  + CB;
    MM[11,11]  = MM[11,11]  + CC;
    MM[12,12]  = MM[12,12]  + CD;
    return MM
end

function local_mass_LUMPED_DIAGONAL_NO_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    # C
    # C  LUMPED DIAGONAL ISOTROPIC MASS MATRIX WITHOUT ROTATIONAL MASSES
    # C
    HLM  = A*rho*L/2.;
    CA = HLM;
    CB = 0.0;
    CC = 0.0;
    CD = 0.0;
    fill!(MM, 0.0);
    MM[1,1]    = MM[1,1]    + CA;
    MM[2,2]    = MM[2,2]    + CA;
    MM[3,3]    = MM[3,3]    + CA;
    MM[7,7]    = MM[7,7]    + CA;
    MM[8,8]    = MM[8,8]    + CA;
    MM[9,9]    = MM[9,9]    + CA;
    return MM
end


"""
    local_mass!(MM, A, I1, I2, I3, rho, L, mass_type)

Mass matrix of the beam.

# Arguments
`A`= cross-sectional area,
`I1`=central moment of inertia of the cross-section about the x1 axis,
`I2`, `I3`=central moment of inertia of the cross-section about the x2 and x3
coordinate axis,
`rho`=mass density,
`L`= initial length of the element,

# Outputs
`MM` = local mass matrix, 12 x 12
In the element frame the mass matrix is constant.
"""
function local_mass!(MM, A, I1, I2, I3, rho, L, mass_type)
    if (mass_type == MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA)
        return local_mass_CONSISTENT_WITH_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    elseif (mass_type == MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA)
        return local_mass_CONSISTENT_NO_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    elseif (mass_type == MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA)
        return local_mass_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    elseif (mass_type == MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA)
        return local_mass_LUMPED_DIAGONAL_NO_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    end
end

function local_mass_original!(MM, A, I1, I2, I3, rho, L, mass_type)
    if (mass_type == MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA)
        # C
        # C  CONSISTENT MASS MATRIX including ROTATIONAL MASSES
        # C  Formulation of the (3.38), (3.39) equation from Dykstra's thesis
        MM .= (rho*A*L)*[
            1/3      0        0        0        0        0       1/6       0        0        0        0        0
            0     13/35       0        0        0    11*L/210     0      9/70       0        0        0    -13*L/420
            0        0     13/35       0    -11*L/210    0        0        0       9/70      0    13*L/420     0
            0        0        0     I1/3/A      0        0        0        0        0      I1/6/A     0        0
            0        0   -11*L/210     0     L^2/105     0        0        0    -13*L/420    0    -L^2/140     0
            0    11*L/210     0        0        0     L^2/105     0    13*L/420     0        0        0    -L^2/140
            1/6      0        0        0        0        0       1/3       0        0        0        0        0
            0      9/70       0        0        0     13*L/420    0      13/35      0        0        0   -11*L/210
            0        0       9/70    0     -13*L/420     0        0        0      13/35      0    11*L/210     0
            0        0        0     I1/6/A      0        0        0        0        0      I1/3/A     0        0
            0        0    13*L/420     0    -L^2/140     0        0        0     11*L/210     0     L^2/105    0
            0   -13*L/420     0        0        0    -L^2/140     0   -11*L/210     0         0        0    L^2/105];
        MM .+= (rho/L)*[
            0       0        0        0       0         0       0       0        0       0       0        0
            0    6/5*I2      0        0       0     L/10*I2     0    -6/5*I2     0       0       0    L/10*I2
            0       0     6/5*I3      0   -L/10*I3      0       0       0     -6/5*I3    0   -L/10*I3     0
            0       0        0        0       0         0       0       0        0       0       0        0
            0       0    -L/10*I3     0  2*L^2/15*I3    0       0       0     L/10*I3    0   -L^2/30*I3   0
            0    L/10*I2     0        0       0   2*L^2/15*I2   0    -L/10*I2    0       0       0    -L^2/30*I2
            0       0        0        0       0         0       0       0        0       0       0        0
            0   -6/5*I2      0        0       0     -L/10*I2    0    6/5*I2      0       0       0     -L/10*I2
            0       0     -6/5*I3     0    L/10*I3      0       0       0      6/5*I3    0    L/10*I3     0
            0       0        0        0       0         0       0       0        0       0       0        0
            0       0     -L/10*I3    0  -L^2/30*I3     0       0       0     L/10*I3    0   2*L^2/15*I3  0
            0   L/10*I2      0        0       0    -L^2/30*I2   0   -L/10*I2     0       0       0   2*L^2/15*I2];
    elseif (mass_type == MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA)
        # C
        # C  CONSISTENT MASS MATRIX excluding ROTATIONAL MASSES
        # C  Formulation of the (3.38) equation from Dykstra's thesis, no rotational inertia
        MM .= (rho*A*L)*[
            1/3      0        0        0        0        0       1/6       0        0        0        0        0
            0     13/35       0        0        0    11*L/210     0      9/70       0        0        0    -13*L/420
            0        0     13/35       0    -11*L/210    0        0        0       9/70      0    13*L/420     0
            0        0        0     I1/3/A      0        0        0        0        0      I1/6/A     0        0
            0        0   -11*L/210     0     L^2/105     0        0        0    -13*L/420    0    -L^2/140     0
            0    11*L/210     0        0        0     L^2/105     0    13*L/420     0        0        0    -L^2/140
            1/6      0        0        0        0        0       1/3       0        0        0        0        0
            0      9/70       0        0        0     13*L/420    0      13/35      0        0        0   -11*L/210
            0        0       9/70    0     -13*L/420     0        0        0      13/35      0    11*L/210     0
            0        0        0     I1/6/A      0        0        0        0        0      I1/3/A     0        0
            0        0    13*L/420     0    -L^2/140     0        0        0     11*L/210     0     L^2/105    0
            0   -13*L/420     0        0        0    -L^2/140     0   -11*L/210     0         0        0    L^2/105];
    elseif (mass_type == MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA)
        # C
        # C  LUMPED DIAGONAL MASS MATRIX WITH ROTATIONAL MASSES
        # C
        HLM  = A*rho*L/2.;
        HLI1 = rho*I1* L/2.;
        HLI2 = rho*I2* L/2.;
        HLI3 = rho*I3* L/2.;
        CA = HLM;
        CB = HLI1;
        CC = HLI2;
        CD = HLI3;
        fill!(MM, 0.0);
        MM[1,1]    = MM[1,1]    + CA;
        MM[2,2]    = MM[2,2]    + CA;
        MM[3,3]    = MM[3,3]    + CA;
        MM[4,4]    = MM[4,4]    + CB;
        MM[5,5]    = MM[5,5]    + CC;
        MM[6,6]    = MM[6,6]    + CD;
        MM[7,7]    = MM[7,7]    + CA;
        MM[8,8]    = MM[8,8]    + CA;
        MM[9,9]    = MM[9,9]    + CA;
        MM[10,10]  = MM[10,10]  + CB;
        MM[11,11]  = MM[11,11]  + CC;
        MM[12,12]  = MM[12,12]  + CD;
    elseif (mass_type == MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA)
        # C
        # C  LUMPED DIAGONAL ISOTROPIC MASS MATRIX WITHOUT ROTATIONAL MASSES
        # C
        HLM  = A*rho*L/2.;
        CA = HLM;
        CB = 0.0;
        CC = 0.0;
        CD = 0.0;
        fill!(MM, 0.0);
        MM[1,1]    = MM[1,1]    + CA;
        MM[2,2]    = MM[2,2]    + CA;
        MM[3,3]    = MM[3,3]    + CA;
        MM[7,7]    = MM[7,7]    + CA;
        MM[8,8]    = MM[8,8]    + CA;
        MM[9,9]    = MM[9,9]    + CA;
    end
    return MM
end


"""
    local_stiffness!(SM, E, G, A, I2, I3, J, A2s, A3s, L, aN, DN)

Compute the local elastic stiffness matrix.

# Arguments
`E`, `G`= Young's and shear modulus,
`A`= cross-sectional area,
`I2`, `I3`=central moment of inertia of the cross-section about the x2 and x3
coordinate axis,
`J`=St Venant torsion constant,
`L`= current length of the element,

# Outputs
`SM` = local stiffness matrix, 12 x 12
"""
function local_stiffness!(SM, E, G, A, I2, I3, J, A2s, A3s, L, aN, DN)
    local_cartesian_to_natural!(aN, L);
    natural_stiffness!(DN, E, G, A, I2, I3, J, A2s, A3s, L);
    SM .= aN'*DN*aN;
    return SM
end

"""
    natural_forces!(PN, E, G, A, I2, I3, J, A2s, A3s, L, dN, DN)

Compute the natural forces from the natural deformations.

# Argument
`E`, `G`= Young's and shear modulus,
`A`= cross-sectional area,
`I2`, `I3`=central moment of inertia of the cross-section about the x2 and x3
coordinate axis,
`J`=St Venant torsion constant,
`L`= current length of the element,
`dN`= column vector of natural deformations; see local_frames()

# Outputs
`PN` = column vector of natural forces;
     `PN[1]`= axial force;
     `PN[2]`= symmetric bending moment in the plane x1-x2;
     `PN[3]`= anti-symmetric bending bending moment in the plane x1-x2;
     `PN[4]`= symmetric bending bending moment in the plane x1-x3;
     `PN[5]`= anti-symmetric bending bending moment in the plane x1-x3;
     `PN[6]`= axial torque.
"""
function natural_forces!(PN, E, G, A, I2, I3, J, A2s, A3s, L, dN, DN)
    natural_stiffness!(DN, E, G, A, I2, I3, J, A2s, A3s, L);
    #     Natural forces
    mul!(PN, DN, dN);
    # Note that the non-constitutive stiffness due to pre-existing internal forces is currently omitted
end

"""
    natural_stiffness_Bernoulli!(DN, E, G, A, I2, I3, J, A2s, A3s, L)

Compute the natural stiffness matrix.

function DN = natural_stiffness(self, E, G, A, I2, I3, J, L)

# Arguments
- `E`, `G`= Young's and shear modulus,
- `A`= cross-sectional area,
- `I2`, `I3`= central moment of inertia of the cross-section about the x2 and x3
  coordinate axis,
- `J`= St Venant torsion constant,
- `A2s` = effective area for shear in the direction of x2 (ignored)
- `A3s` = effective area for shear in the direction of x3 (ignored)
- `L`= current length of the element,

# Outputs
- `DN` = 6 x 6 natural stiffness matrix
"""
function natural_stiffness_Bernoulli!(DN, E, G, A, I2, I3, J, A2s, A3s, L)
    fill!(DN, 0.0)
    DN[1, 1] = E*A/L
    DN[2, 2] = E*I3/L
    DN[3, 3] = 3*E*I3/L
    DN[4, 4] = E*I2/L
    DN[5, 5] = 3*E*I2/L
    DN[6, 6] = G*J/L
    return DN
end

"""
    natural_stiffness_Timoshenko!(DN, E, G, A, I2, I3, J, A2s, A3s, L)

Compute the natural stiffness matrix.

function DN = natural_stiffness(self, E, G, A, I2, I3, J, L)

# Arguments
- `E`, `G`= Young's and shear modulus,
- `A`= cross-sectional area,
- `I2`, `I3`= central moment of inertia of the cross-section about the x2 and x3
  coordinate axis,
- `J`= St Venant torsion constant,
- `A2s` = effective area for shear in the direction of x2
- `A3s` = effective area for shear in the direction of x3
- `L`= current length of the element,

# Outputs
- `DN` = 6 x 6 natural stiffness matrix
"""
function natural_stiffness_Timoshenko!(DN, E, G, A, I2, I3, J, A2s, A3s, L)
    fill!(DN, 0.0)
    DN[1, 1] = E*A/L
    DN[2, 2] = E*I3/L
    Phi3 = 12 * E * I3 / (G * A2s * L^2)
    DN[3, 3] = 3*E*I3/L/(1+Phi3)
    DN[4, 4] = E*I2/L
    Phi2 = 12 * E * I2 / (G * A3s * L^2)
    DN[5, 5] = 3*E*I2/L/(1+Phi2)
    DN[6, 6] = G*J/L
    return DN
end

function natural_stiffness!(DN, E, G, A, I2, I3, J, A2s, A3s, L)
    if A2s == Inf || A3s == Inf
        return natural_stiffness_Bernoulli!(DN, E, G, A, I2, I3, J, A2s, A3s, L)
    else
        return natural_stiffness_Timoshenko!(DN, E, G, A, I2, I3, J, A2s, A3s, L)
    end
end

"""
    local_forces!(FL, PN, L, aN)

Compute forces through which the element acts on the nodes in the
local coordinate system.

# Arguments
`PN` = column vector of natural forces;
`L`= current length of the element,
`aN`= transformation matrix to take Cartesian (local) displacement increments in the
     element frame and to produce increments of natural deformations;
     see `local_frame_and_def!` for the definition of the natural deformations

# Outputs
FL = vector of forces acting on the nodes in the local coordinate system
"""
function local_forces!(FL, PN, L, aN)
    local_cartesian_to_natural!(aN, L);
    mul!(FL, Transpose(aN), PN);
    return FL
end

"""
    mass(self::FEMMCorotBeam,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the consistent mass matrix

This is a general routine for the abstract linear-deformation  FEMM.
"""
function mass(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    dN = self._dN
    rho = massdensity(self.material)
    A, I1, I2, I3, x1x2_vector = fes.A, fes.I1, fes.I2, fes.I3, fes.x1x2_vector
    startassembly!(assembler, prod(size(elmat)) * count(fes), nalldofs(dchi), nalldofs(dchi))
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        R1I[:] .= Rfield1.values[fes.conn[i][1], :];
        R1J[:] .= Rfield1.values[fes.conn[i][2], :];
        fill!(elmat,  0.0); # Initialize element matrix
        L1, Ft, dN = local_frame_and_def!(Ft, dN, F0, FtI, FtJ, ecoords0, x1x2_vector[i], ecoords1, R1I, R1J);
        _transfmat!(Te, Ft)
        L0 = norm(ecoords0[2,:]-ecoords0[1,:]); 
        local_mass!(elmat, A[i], I1[i], I2[i], I3[i], rho, L0, mass_type);
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function mass(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return mass(self, assembler, geom0, u1, Rfield1, dchi; mass_type = mass_type);
end

"""
    gyroscopic(self::FEMMCorotBeam,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the quadratic-inertial-term (gyroscopic) mass matrix

This is a general routine for the abstract linear-deformation  FEMM.
"""
function gyroscopic(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, v1::NodalField{T}, dchi::NodalField{TI}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    dN = self._dN
    evel1, evel1f = self._evel1, self._evel1f
    Ge1, Ge2, Ge = self._tempelmat1, self._tempelmat2, self._tempelmat3
    OmegaTilde = self._elmato
    OS = self._OS
    rho = massdensity(self.material)
    A, I1, I2, I3, x1x2_vector = fes.A, fes.I1, fes.I2, fes.I3, fes.x1x2_vector
    startassembly!(assembler, prod(size(elmat)) * count(fes), nalldofs(dchi), nalldofs(dchi))
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        R1I[:] .= Rfield1.values[fes.conn[i][1], :];
        R1J[:] .= Rfield1.values[fes.conn[i][2], :];
        fill!(elmat,  0.0); # Initialize element matrix
        L1, Ft, dN = local_frame_and_def!(Ft, dN, F0, FtI, FtJ, ecoords0, x1x2_vector[i], ecoords1, R1I, R1J);
        _transfmat!(Te, Ft)
        L0 = norm(ecoords0[2,:]-ecoords0[1,:]); 
        local_mass!(elmat, A[i], I1[i], I2[i], I3[i], rho, L0, mass_type);
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gathervalues_asmat!(v1, evel1, fes.conn[i]);
        evel1f[:, 1:3] = evel1[:,1:3]*Ft
        evel1f[:, 4:6] = evel1[:,4:6]*Ft
        Omega = (evel1f[1,4]+evel1f[2,4])/2*Ft[:,1] + (evel1f[1,3]-evel1f[2,3])/L1*Ft[:,2] + (evel1f[2,2]-evel1f[1,2])/L1*Ft[:,3];
        skewmat!(OS, Omega);
        _transfmat!(OmegaTilde, OS)
        mul!(Ge1, OmegaTilde, elmat)
        mul!(Ge2, elmat, OmegaTilde)
        @. Ge = Ge1 - Ge2;
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, Ge, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function gyroscopic(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, v1::NodalField{T}, dchi::NodalField{TI}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return gyroscopic(self, assembler, geom0, u1, Rfield1, v1, dchi; mass_type = mass_type);
end

"""
    stiffness(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the material stiffness matrix.
"""
function stiffness(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    aN, dN, DN = self._aN, self._dN, self._DN
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)::Float64
    A, I2, I3, J, A2s, A3s, x1x2_vector = fes.A, fes.I2, fes.I3, fes.J, fes.A2s, fes.A3s, fes.x1x2_vector
    startassembly!(assembler, prod(size(elmat)) * count(fes), nalldofs(dchi), nalldofs(dchi))
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        R1I[:] .= Rfield1.values[fes.conn[i][1], :];
        R1J[:] .= Rfield1.values[fes.conn[i][2], :];
        fill!(elmat,  0.0); # Initialize element matrix
        L1, Ft, dN = local_frame_and_def!(Ft, dN, F0, FtI, FtJ, ecoords0, x1x2_vector[i], ecoords1, R1I, R1J);
        _transfmat!(Te, Ft)
        local_stiffness!(elmat, E, G, A[i], I2[i], I3[i], J[i], A2s[i], A3s[i], L1, aN, DN);
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function stiffness(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi);
end


"""
    geostiffness(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the geometric stiffness matrix.
"""
function geostiffness(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    aN, dN, DN, PN = self._aN, self._dN, self._DN, self._PN
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)
    A, I2, I3, J, A2s, A3s, x1x2_vector = fes.A, fes.I2, fes.I3, fes.J, fes.A2s, fes.A3s, fes.x1x2_vector
    startassembly!(assembler, prod(size(elmat)) * count(fes), nalldofs(dchi), nalldofs(dchi))
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        R1I[:] .= Rfield1.values[fes.conn[i][1], :];
        R1J[:] .= Rfield1.values[fes.conn[i][2], :];
        fill!(elmat,  0.0); # Initialize element matrix
        L1, Ft, dN = local_frame_and_def!(Ft, dN, F0, FtI, FtJ, ecoords0, x1x2_vector[i], ecoords1, R1I, R1J);
        _transfmat!(Te, Ft)
        natural_forces!(PN, E, G, A[i], I2[i], I3[i], J[i], A2s[i], A3s[i], L1, dN, DN)
        local_geometric_stiffness!(elmat, A[i], I2[i], I3[i], PN, L1)
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function geostiffness(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return geostiffness(self, assembler, geom0, u1, Rfield1, dchi);
end

"""
    restoringforce(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the vector of the restoring elastic forces
"""
function restoringforce(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysvecAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    aN, dN, DN, PN, LF = self._aN, self._dN, self._DN, self._PN, self._LF
    elvec = self._elvec
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)
    A, I2, I3, J, A2s, A3s, x1x2_vector = fes.A, fes.I2, fes.I3, fes.J, fes.A2s, fes.A3s, fes.x1x2_vector
    startassembly!(assembler, nalldofs(dchi));
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        R1I[:] .= Rfield1.values[fes.conn[i][1], :];
        R1J[:] .= Rfield1.values[fes.conn[i][2], :];
        fill!(elmat,  0.0); # Initialize element matrix
        L1, Ft, dN = local_frame_and_def!(Ft, dN, F0, FtI, FtJ, ecoords0, x1x2_vector[i], ecoords1, R1I, R1J);
        _transfmat!(Te, Ft)
        natural_forces!(PN, E, G, A[i], I2[i], I3[i], J[i], A2s[i], A3s[i], L1, dN, DN)
        local_forces!(LF, PN, L1, aN)
        mul!(elvec, Te, -LF)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elvec, dofnums); 
    end # Loop over elements
    return makevector!(assembler);
end

function restoringforce(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {T<:Number, TI<:Number}
    assembler = SysvecAssembler();
    return restoringforce(self, assembler, geom0, u1, Rfield1, dchi);
end

"""
    distribloads_global(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}, fi) where {T<:Number}
    
Compute the load vector due to distributed loads.

Compute the global load vector corresponding to applied distributed
load. Here it means force per unit length of the beam,
in the configuration u1, Rfield1. These are only forces, not moments.

!!! note

    The force intensity must be uniform across the entire element. The
    force intensity is given in the global coordinates.
"""
function distribloads_global(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}, fi) where {ASS<:AbstractSysvecAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    aN, dN, DN, PN = self._aN, self._dN, self._DN, self._PN
    elvec, elvecf = self._elvec, self._elvecf
    Lforce = fill(0.0, 3)
    ignore = fill(0.0 , 0, 0)
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)
    A, I2, I3, J, x1x2_vector = fes.A, fes.I2, fes.I3, fes.J, fes.x1x2_vector
    startassembly!(assembler, nalldofs(dchi))
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        R1I[:] .= Rfield1.values[fes.conn[i][1], :];
        R1J[:] .= Rfield1.values[fes.conn[i][2], :];
        fill!(elmat,  0.0); # Initialize element matrix
        L1, Ft, dN = local_frame_and_def!(Ft, dN, F0, FtI, FtJ, ecoords0, x1x2_vector[i], ecoords1, R1I, R1J);
        _transfmat!(Te, Ft)
        L0 = norm(ecoords0[2,:]-ecoords0[1,:]); 
        force = updateforce!(fi, ignore, ignore, i, 0); # retrieve the applied load
        Lforce = Ft' * force # local distributed load components
        elvecf[1] = Lforce[1]*L0/2;
        elvecf[2] = Lforce[2]*L0/2;
        elvecf[3] = Lforce[3]*L0/2;
        elvecf[4] = 0;
        elvecf[5] = -Lforce[3]*L0^2/12;
        elvecf[6] = +Lforce[2]*L0^2/12;
        elvecf[7] = Lforce[1]*L0/2;
        elvecf[8] = Lforce[2]*L0/2;
        elvecf[9] = Lforce[3]*L0/2;
        elvecf[10] = 0;
        elvecf[11] = +Lforce[3]*L0^2/12;
        elvecf[12] = -Lforce[2]*L0^2/12;
        mul!(elvec, Te, elvecf)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elvec, dofnums); 
    end # Loop over elements
    return makevector!(assembler);
end

function distribloads_global(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}, fi) where {T<:Number, TI<:Number}
    assembler = SysvecAssembler();
    return distribloads_global(self, assembler, geom0, u1, Rfield1, dchi, fi);
end

end # module
