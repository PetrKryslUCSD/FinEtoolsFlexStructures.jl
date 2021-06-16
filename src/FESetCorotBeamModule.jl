module FESetCorotBeamModule

using FinEtools
import FinEtools: cat, subset
using FinEtools.MatrixUtilityModule: complete_lt! 
using ..CrossSectionModule: AbstractCrossSectionType
using LinearAlgebra: norm, Transpose, mul!

mutable struct FESetL2CorotBeam{CT} <: AbstractFESet1Manifold{2}
    conn::Array{NTuple{2, FInt}, 1};
    label::FIntVec; 
    crosssection::CT
    A::FFltVec
    I1::FFltVec
    I2::FFltVec
    I3::FFltVec
    J::FFltVec
    A2s::FFltVec
    A3s::FFltVec
    x1x2_vector::Vector{FFltVec}
    dimensions::Vector{FFltVec}
end

function FESetL2CorotBeam(conn::FIntMat, crosssection::CT) where {CT}
    par = crosssection.parameters(0.0)
    N = size(conn, 1)
    _A = fill(par.A, N)
    _I1 = fill(par.I1, N)
    _I2 = fill(par.I2, N)
    _I3 = fill(par.I3, N)
    _J = fill(par.J, N)
    _A2s = fill(par.A2s, N) 
    _A3s = fill(par.A3s, N)
    _x1x2_vector = [par.x1x2_vector for i in 1:N]
    _dimensions = [par.dimensions for i in 1:N]
    self = FESetL2CorotBeam(NTuple{2, FInt}[], FInt[], crosssection, _A, _I1, _I2, _I3, _J, _A2s, _A3s, _x1x2_vector, _dimensions)
    self = fromarray!(self, conn)
    setlabel!(self, 0)
    return self
end

function FESetL2CorotBeam(conn::FIntMat, crosssection::CT, _A, _I1, _I2, _I3, _J, _A2s, _A3s, _x1x2_vector, _dimensions) where {CT}
    dummy = FESetL2(conn)
    setlabel!(dummy, 0)
    return FESetL2CorotBeam(dummy.conn, dummy.label, crosssection, _A, _I1, _I2, _I3, _J, _A2s, _A3s, _x1x2_vector, _dimensions)
end

"""
    Type of the mass matrix formulation.
"""
const MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA=0;
const MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA=1;
const MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA=2;
const MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA=3;

"""
    cat(self::T,  other::T) where {T<:FESetL2CorotBeam}

Concatenate two sets of beam elements.
"""
function cat(self::T,  other::T) where {T<:FESetL2CorotBeam}
    @assert self.crosssection === other.crosssection "Cannot concatenate sets with distinct cross-sections"
    result = deepcopy(self)
    result.conn = vcat(self.conn, other.conn);
    setlabel!(result, vcat(self.label, other.label))
    result.A = vcat(self.A, other.A);
    result.I1 = vcat(self.I1, other.I1);
    result.I2 = vcat(self.I2, other.I2);
    result.I3 = vcat(self.I3, other.I3);
    result.J = vcat(self.J, other.J);
    result.A2s = vcat(self.A2s, other.A2s);
    result.A3s = vcat(self.A3s, other.A3s);
    result.x1x2_vector = vcat(self.x1x2_vector, other.x1x2_vector);
    result.dimensions = vcat(self.dimensions, other.dimensions);
    return result
end

"""
    subset(self::T, L::FIntVec) where {T<:FESetL2CorotBeam}

Subset of a beam-element set
"""
function subset(self::T, L::FIntVec) where {T<:FESetL2CorotBeam}
    result = deepcopy(self)
    result.conn = deepcopy(self.conn[L])
    result.label = deepcopy(self.label[L])
    result.A = self.A[L] 
    result.I1 = self.I1[L]
    result.I2 = self.I2[L]
    result.I3 = self.I3[L]
    result.J = self.J[L] 
    result.A2s = self.A2s[L] 
    result.A3s = self.A3s[L] 
    result.x1x2_vector = self.x1x2_vector[L]
    result.dimensions = self.dimensions[L]
    return  result
end

function initial_local_frame!(F0, x0, x1x2_vector)
    # This is the element frame in the configuration t=0
    F0[:,1] = (x0[2,:]-x0[1,:]); 
    L0 = norm(@view F0[:,1]);
    F0[:,1] /= L0;
    #     F0(:,3)=skewmat(F0(:,1))*x1x2_vector;
    F0[:,3] .= (-F0[3,1]*x1x2_vector[2]+F0[2,1]*x1x2_vector[3],
                 F0[3,1]*x1x2_vector[1]-F0[1,1]*x1x2_vector[3],
                -F0[2,1]*x1x2_vector[1]+F0[1,1]*x1x2_vector[2]);
    F0[:,3] /= norm(@view F0[:,3]);
    #     F0(:,2)=skewmat(F0(:,3))*F0(:,1);
    F0[:,2] .= (-F0[3,3]*F0[2,1]+F0[2,3]*F0[3,1],
                 F0[3,3]*F0[1,1]-F0[1,3]*F0[3,1],
                -F0[2,3]*F0[1,1]+F0[1,3]*F0[2,1]);
    return L0, F0
end

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
    DN[5, 5] = 3*E*I2/L/(1+Phi3)
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


end # module
