"""
Module for construction of the linear algebra matrix and vector quantities for
the linear beam with eccentric connections to nodes.
"""
module FEMMLinBeamModule

using LinearAlgebra: norm, Transpose, mul!
using LinearAlgebra
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.MatrixUtilityModule: complete_lt!
using FinEtools.IntegDomainModule: IntegDomain
import FinEtoolsDeforLinear.MatDeforElastIsoModule: MatDeforElastIso
using ..FESetL2BeamModule: FESetL2Beam, initial_local_frame!
import FinEtools.FEMMBaseModule: inspectintegpoints

"""
    mutable struct FEMMLinBeam{ID<:IntegDomain{S} where {S<:FESetL2}, M<:MatDeforElastIso} <: AbstractFEMM

Type for linear beam finite element modeling machine.

The beam can be connected to the nodes with given eccentricities (the transverse
eccentricity is uniform along its length). Only linear kinematics is
implemented at the moment. The beam stiffness can be either shear-flexible
(Timoshenko), or shear-rigid (Bernoulli-Euler). The local beam stiffness is
expressed analytically: no numerical integration is involved.
"""
mutable struct FEMMLinBeam{ID<:IntegDomain{S} where {S<:FESetL2}, M<:MatDeforElastIso} <: AbstractFEMM
    integdomain::ID # integration domain data
    material::M # material object
    eccentricities::FFltMat
    # The attributes below are buffers used in various operations.
    _ecoords0::FFltMat
    _dofnums::FIntMat
    _F0::FFltMat
    _Ft::FFltMat
    _FtI::FFltMat
    _FtJ::FFltMat
    _Te::FFltMat
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
end

"""
    FEMMLinBeam(
        integdomain::ID,
        material::M,
        uniform_eccentricity,
    ) where {ID<:IntegDomain{S} where {S<:FESetL2}, M<:MatDeforElastIso}

Constructor.

Supply the integration domain, material, and the eccentricity parameters
(given in the order: first node, f1 direction, second node, f1 direction, f2
direction, f3 direction).
"""
function FEMMLinBeam(
    integdomain::ID,
    material::M,
    uniform_eccentricity,
) where {ID<:IntegDomain{S} where {S<:FESetL2}, M<:MatDeforElastIso}
    typeof(delegateof(integdomain.fes)) <: FESetL2Beam ||
        error("Expected to delegate to FESetL2Beam")
    _eccentricities = fill(0.0, count(integdomain.fes), 4)
    _eccentricities[:, 1] .= uniform_eccentricity[1] # first node, f1 direction
    _eccentricities[:, 2] .= uniform_eccentricity[2] # second node, f1 direction
    _eccentricities[:, 3] .= uniform_eccentricity[3] # f2 direction
    _eccentricities[:, 4] .= uniform_eccentricity[4] # f3 direction
    _ecoords0 = fill(0.0, 2, 3)
    _dofnums = zeros(FInt, 1, 12)
    _F0 = fill(0.0, 3, 3)
    _Ft = fill(0.0, 3, 3)
    _FtI = fill(0.0, 3, 3)
    _FtJ = fill(0.0, 3, 3)
    _Te = fill(0.0, 12, 12)
    _elmat = fill(0.0, 12, 12)
    _elmatTe = fill(0.0, 12, 12)
    _elmato = fill(0.0, 12, 12)
    _elvec = fill(0.0, 12)
    _elvecf = fill(0.0, 12)
    _aN = fill(0.0, 6, 12)
    _dN = fill(0.0, 6)
    _DN = fill(0.0, 6, 6)
    _PN = fill(0.0, 6)
    _LF = fill(0.0, 12)
    return FEMMLinBeam(
        integdomain,
        material,
        _eccentricities,
        _ecoords0,
        _dofnums,
        _F0,
        _Ft,
        _FtI,
        _FtJ,
        _Te,
        _elmat,
        _elmatTe,
        _elmato,
        _elvec,
        _elvecf,
        _aN,
        _dN,
        _DN,
        _PN,
        _LF,
    )
end

"""
    FEMMLinBeam(
        integdomain::ID,
        material::M,
    ) where {ID<:IntegDomain{S} where {S<:FESetL2}, M<:MatDeforElastIso}

Constructor.

Supply the integration domain and the material. The eccentricities are assumed to be zero.
"""
function FEMMLinBeam(
    integdomain::ID,
    material::M,
) where {ID<:IntegDomain{S} where {S<:FESetL2}, M<:MatDeforElastIso}
    return FEMMLinBeam(integdomain, material, [0.0, 0.0, 0.0, 0.0])
end

function properties(fes)
    d = delegateof(fes)
    # A, I1, I2, I3, J, A2s, A3s, x1x2_vector, dimensions = properties(fes)
    d.A, d.I1, d.I2, d.I3, d.J, d.A2s, d.A3s, d.x1x2_vector, d.dimensions
end

function _transfmat!(Te, Ft)
    Te[1:3, 1:3] = Te[4:6, 4:6] = Te[7:9, 7:9] = Te[10:12, 10:12] = Ft
    return Te
end


"""
    Type of the mass matrix formulation.
"""
const MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA = 0
const MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA = 1
const MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA = 2
const MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA = 3

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
    aN[1, 1] = -1
    aN[1, 7] = +1
    aN[2, 6] = +1
    aN[2, 12] = -1
    aN[3, 2] = 2 / L
    aN[3, 6] = +1
    aN[3, 8] = -2 / L
    aN[3, 12] = +1
    aN[4, 5] = -1
    aN[4, 11] = +1
    aN[5, 3] = 2 / L
    aN[5, 5] = -1
    aN[5, 9] = -2 / L
    aN[5, 11] = -1
    aN[6, 4] = -1
    aN[6, 10] = +1
    # aN=[[ -1,   0,   0,  0,  0, 0, 1,    0,    0, 0,  0,  0]
    #     [  0,   0,   0,  0,  0, 1, 0,    0,    0, 0,  0, -1]
    #     [  0, 2/L,   0,  0,  0, 1, 0, -2/L,    0, 0,  0,  1]
    #     [  0,   0,   0,  0, -1, 0, 0,    0,    0, 0,  1,  0]
    #     [  0,   0, 2/L,  0, -1, 0, 0,    0, -2/L, 0, -1,  0]
    #     [  0,   0,   0, -1,  0, 0, 0,    0,    0, 1,  0,  0]];
    return aN
end

function local_mass_CONSISTENT_NO_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    # C
    # C  CONSISTENT MASS MATRIX including ROTATIONAL MASSES
    # C  Formulation of the (3.38), (3.39) equation from Dykstra's thesis
    fill!(MM, 0.0)
    c1 = (rho * A * L)::Float64
    MM[1, 1] = c1 * (1 / 3)
    MM[1, 7] = c1 * (1 / 6)
    MM[2, 2] = c1 * (13 / 35)
    MM[2, 6] = c1 * (11 * L / 210)
    MM[2, 8] = c1 * (9 / 70)
    MM[2, 12] = c1 * (-13 * L / 420)
    MM[3, 3] = c1 * (13 / 35)
    MM[3, 5] = c1 * (-11 * L / 210)
    MM[3, 9] = c1 * (9 / 70)
    MM[3, 11] = c1 * (13 * L / 420)
    MM[4, 4] = c1 * (I1 / 3 / A)
    MM[4, 10] = c1 * (I1 / 6 / A)
    MM[5, 5] = c1 * (L^2 / 105)
    MM[5, 9] = c1 * (-13 * L / 420)
    MM[5, 11] = c1 * (-L^2 / 140)
    MM[6, 6] = c1 * (L^2 / 105)
    MM[6, 8] = c1 * (13 * L / 420)
    MM[6, 12] = c1 * (-L^2 / 140)
    MM[7, 7] = c1 * (1 / 3)
    MM[8, 8] = c1 * (13 / 35)
    MM[8, 12] = c1 * (-11 * L / 210)
    MM[9, 9] = c1 * (13 / 35)
    MM[9, 11] = c1 * (11 * L / 210)
    MM[10, 10] = c1 * (I1 / 3 / A)
    MM[11, 11] = c1 * (L^2 / 105)
    MM[12, 12] = c1 * (L^2 / 105)
    # 1       2        3        4       5         6       7        8       9       10      11       12
    complete_lt!(MM)
    return MM
end

function local_mass_CONSISTENT_WITH_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    # C
    # C  CONSISTENT MASS MATRIX including ROTATIONAL MASSES
    # C  Formulation of the (3.38), (3.39) equation from Dykstra's thesis
    fill!(MM, 0.0)
    c1 = (rho * A * L)::Float64
    MM[1, 1] = c1 * (1 / 3)
    MM[1, 7] = c1 * (1 / 6)
    MM[2, 2] = c1 * (13 / 35)
    MM[2, 6] = c1 * (11 * L / 210)
    MM[2, 8] = c1 * (9 / 70)
    MM[2, 12] = c1 * (-13 * L / 420)
    MM[3, 3] = c1 * (13 / 35)
    MM[3, 5] = c1 * (-11 * L / 210)
    MM[3, 9] = c1 * (9 / 70)
    MM[3, 11] = c1 * (13 * L / 420)
    MM[4, 4] = c1 * (I1 / 3 / A)
    MM[4, 10] = c1 * (I1 / 6 / A)
    MM[5, 5] = c1 * (L^2 / 105)
    MM[5, 9] = c1 * (-13 * L / 420)
    MM[5, 11] = c1 * (-L^2 / 140)
    MM[6, 6] = c1 * (L^2 / 105)
    MM[6, 8] = c1 * (13 * L / 420)
    MM[6, 12] = c1 * (-L^2 / 140)
    MM[7, 7] = c1 * (1 / 3)
    MM[8, 8] = c1 * (13 / 35)
    MM[8, 12] = c1 * (-11 * L / 210)
    MM[9, 9] = c1 * (13 / 35)
    MM[9, 11] = c1 * (11 * L / 210)
    MM[10, 10] = c1 * (I1 / 3 / A)
    MM[11, 11] = c1 * (L^2 / 105)
    MM[12, 12] = c1 * (L^2 / 105)
    c2 = (rho / L)::Float64
    MM[2, 2] += c2 * (6 / 5 * I2)
    MM[2, 6] += c2 * (L / 10 * I2)
    MM[2, 8] += c2 * (-6 / 5 * I2)
    MM[2, 12] += c2 * (L / 10 * I2)
    MM[3, 3] += c2 * (6 / 5 * I3)
    MM[3, 5] += c2 * (-L / 10 * I3)
    MM[3, 9] += c2 * (-6 / 5 * I3)
    MM[3, 11] += c2 * (-L / 10 * I3)
    MM[5, 5] += c2 * (2 * L^2 / 15 * I3)
    MM[5, 9] += c2 * (L / 10 * I3)
    MM[5, 11] += c2 * (-L^2 / 30 * I3)
    MM[6, 6] += c2 * (2 * L^2 / 15 * I2)
    MM[6, 8] += c2 * (-L / 10 * I2)
    MM[6, 12] += c2 * (-L^2 / 30 * I2)
    MM[8, 8] += c2 * (6 / 5 * I2)
    MM[8, 12] += c2 * (-L / 10 * I2)
    MM[9, 9] += c2 * (6 / 5 * I3)
    MM[9, 11] += c2 * (L / 10 * I3)
    MM[11, 11] += c2 * (2 * L^2 / 15 * I3)
    MM[12, 12] += c2 * (2 * L^2 / 15 * I2)
    # 1       2        3        4       5         6       7        8       9       10      11       12
    complete_lt!(MM)
    return MM
end

function local_mass_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    # C
    # C  LUMPED DIAGONAL MASS MATRIX WITH ROTATIONAL MASSES
    # C
    HLM = A * rho * L / 2.0
    HLI1 = rho * I1 * L / 2.0
    HLI2 = rho * I2 * L / 2.0
    HLI3 = rho * I3 * L / 2.0
    CA = HLM
    CB = HLI1
    CC = HLI2
    CD = HLI3
    fill!(MM, 0.0)
    MM[1, 1] = MM[1, 1] + CA
    MM[2, 2] = MM[2, 2] + CA
    MM[3, 3] = MM[3, 3] + CA
    MM[4, 4] = MM[4, 4] + CB
    MM[5, 5] = MM[5, 5] + CC
    MM[6, 6] = MM[6, 6] + CD
    MM[7, 7] = MM[7, 7] + CA
    MM[8, 8] = MM[8, 8] + CA
    MM[9, 9] = MM[9, 9] + CA
    MM[10, 10] = MM[10, 10] + CB
    MM[11, 11] = MM[11, 11] + CC
    MM[12, 12] = MM[12, 12] + CD
    return MM
end

function local_mass_LUMPED_DIAGONAL_NO_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    # C
    # C  LUMPED DIAGONAL ISOTROPIC MASS MATRIX WITHOUT ROTATIONAL MASSES
    # C
    HLM = A * rho * L / 2.0
    CA = HLM
    CB = 0.0
    CC = 0.0
    CD = 0.0
    fill!(MM, 0.0)
    MM[1, 1] = MM[1, 1] + CA
    MM[2, 2] = MM[2, 2] + CA
    MM[3, 3] = MM[3, 3] + CA
    MM[7, 7] = MM[7, 7] + CA
    MM[8, 8] = MM[8, 8] + CA
    MM[9, 9] = MM[9, 9] + CA
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
    local_cartesian_to_natural!(aN, L)
    natural_stiffness!(DN, E, G, A, I2, I3, J, A2s, A3s, L)
    SM .= aN' * DN * aN
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
    natural_stiffness!(DN, E, G, A, I2, I3, J, A2s, A3s, L)
    #     Natural forces
    mul!(PN, DN, dN)
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
    DN[1, 1] = E * A / L
    DN[2, 2] = E * I3 / L
    DN[3, 3] = 3 * E * I3 / L
    DN[4, 4] = E * I2 / L
    DN[5, 5] = 3 * E * I2 / L
    DN[6, 6] = G * J / L
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
    DN[1, 1] = E * A / L
    DN[2, 2] = E * I3 / L
    Phi3 = 12 * E * I3 / (G * A2s * L^2)
    DN[3, 3] = 3 * E * I3 / L / (1 + Phi3)
    DN[4, 4] = E * I2 / L
    Phi2 = 12 * E * I2 / (G * A3s * L^2)
    DN[5, 5] = 3 * E * I2 / L / (1 + Phi2)
    DN[6, 6] = G * J / L
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
    local_cartesian_to_natural!(aN, L)
    mul!(FL, Transpose(aN), PN)
    return FL
end

function _eccentricitytransformation(
    Te,
    eccentricity_f1_1,
    eccentricity_f1_2,
    eccentricity_f2,
    eccentricity_f3,
)
    Te .= zero(eltype(Te))
    Te[1:3, 1:3] = Te[4:6, 4:6] = Te[7:9, 7:9] = Te[10:12, 10:12] = LinearAlgebra.I(3)
    Te[1, 4:6] .= (0.0, eccentricity_f3, -eccentricity_f2)
    Te[2, 4:6] .= (-eccentricity_f3, 0.0, eccentricity_f1_1)
    Te[3, 4:6] .= (eccentricity_f2, -eccentricity_f1_1, 0.0)
    Te[7, 10:12] .= (0.0, eccentricity_f3, -eccentricity_f2)
    Te[8, 10:12] .= (-eccentricity_f3, 0.0, eccentricity_f1_2)
    Te[9, 10:12] .= (eccentricity_f2, -eccentricity_f1_2, 0.0)
    return Te
end

"""
    mass(self::FEMMLinBeam,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the consistent mass matrix

This is a general routine for the abstract linear-deformation  FEMM.
"""
function mass(
    self::FEMMLinBeam,
    assembler::ASS,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI};
    mass_type = MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA,
) where {ASS<:AbstractSysmatAssembler,T<:Number,TI<:Number}
    fes = self.integdomain.fes
    ecoords0, dofnums = self._ecoords0, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    elmat, elmattemp, Secc = self._elmat, self._elmatTe, self._elmato
    eccentricities = self.eccentricities
    dN = self._dN
    rho = massdensity(self.material)
    A, I1, I2, I3, J, A2s, A3s, x1x2_vector, dimensions = properties(fes)
    startassembly!(
        assembler,
        size(elmat)..., count(fes),
        nalldofs(dchi),
        nalldofs(dchi),
    )
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i])
        fill!(elmat, 0.0) # Initialize element matrix
        L0, F0 = initial_local_frame!(F0, ecoords0, x1x2_vector[i])
        _transfmat!(Te, F0)
        L0 = norm(ecoords0[2, :] - ecoords0[1, :])
        _eccentricitytransformation(Secc, eccentricities[i, :]...)
        local_mass!(elmat, A[i], I1[i], I2[i], I3[i], rho, L0, mass_type)
        # First transform from slave to master (U_s = Secc * U_m)
        mul!(elmattemp, elmat, Secc)
        mul!(elmat, Transpose(Secc), elmattemp)
        # Then transform from the master local coordinates into global coordinates
        mul!(elmattemp, elmat, Transpose(Te))
        mul!(elmat, Te, elmattemp)
        gatherdofnums!(dchi, dofnums, fes.conn[i]) # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)
    end # Loop over elements
    return makematrix!(assembler)
end

function mass(
    self::FEMMLinBeam,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI};
    mass_type = MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA,
) where {T<:Number,TI<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return mass(self, assembler, geom0, u1, Rfield1, dchi; mass_type = mass_type)
end

"""
    stiffness(self::FEMMLinBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the material stiffness matrix.
"""
function stiffness(
    self::FEMMLinBeam,
    assembler::ASS,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {ASS<:AbstractSysmatAssembler,T<:Number,TI<:Number}
    fes = self.integdomain.fes
    ecoords0, dofnums = self._ecoords0, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    elmat, elmattemp, Secc = self._elmat, self._elmatTe, self._elmato
    aN, dN, DN = self._aN, self._dN, self._DN
    eccentricities = self.eccentricities
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)::Float64
    A, I1, I2, I3, J, A2s, A3s, x1x2_vector, dimensions = properties(fes)
    startassembly!(
        assembler,
        size(elmat)..., count(fes),
        nalldofs(dchi),
        nalldofs(dchi),
    )
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i])
        fill!(elmat, 0.0) # Initialize element matrix
        L0, F0 = initial_local_frame!(F0, ecoords0, x1x2_vector[i])
        _transfmat!(Te, F0)
        _eccentricitytransformation(Secc, eccentricities[i, :]...)
        local_stiffness!(elmat, E, G, A[i], I2[i], I3[i], J[i], A2s[i], A3s[i], L0, aN, DN)
        # First transform from slave to master (U_s = Secc * U_m)
        mul!(elmattemp, elmat, Secc)
        mul!(elmat, Transpose(Secc), elmattemp)
        # Then transform from the master local coordinates into global coordinates
        mul!(elmattemp, elmat, Transpose(Te))
        mul!(elmat, Te, elmattemp)
        gatherdofnums!(dchi, dofnums, fes.conn[i]) # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)
    end # Loop over elements
    return makematrix!(assembler)
end

function stiffness(
    self::FEMMLinBeam,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {T<:Number,TI<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi)
end


"""
    distribloads_global(self::FEMMLinBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}, fi) where {T<:Number}
    
Compute the load vector due to distributed loads.

Compute the global load vector corresponding to applied distributed
load. Here it means force per unit length of the beam,
in the configuration u1, Rfield1. These are only forces, not moments.

!!! note

    The force intensity must be uniform across the entire element. The
    force intensity is given in the global coordinates.
"""
function distribloads_global(
    self::FEMMLinBeam,
    assembler::ASS,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
    fi,
) where {ASS<:AbstractSysvecAssembler,T<:Number,TI<:Number}
    fes = self.integdomain.fes
    ecoords0, dofnums = self._ecoords0, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    eccentricities = self.eccentricities
    elmat, elmatTe, Secc = self._elmat, self._elmatTe, self._elmato
    disp_g, disp_f_m, forc_f_s = self._elvec, self._elvecf, self._LF
    forc_f_m = deepcopy(forc_f_s)
    forc_f_s = deepcopy(forc_f_m)
    forc_g = deepcopy(forc_f_m)
    Lforce = fill(0.0, 3)
    ignore = fill(0.0, 0, 0)
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)
    A, I1, I2, I3, J, A2s, A3s, x1x2_vector, dimensions = properties(fes)
    startassembly!(assembler, nalldofs(dchi))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i])
        L0, F0 = initial_local_frame!(F0, ecoords0, x1x2_vector[i])
        _transfmat!(Te, F0)
        _eccentricitytransformation(Secc, eccentricities[i, :]...)
        force = updateforce!(fi, ignore, ignore, i, 0) # retrieve the applied load
        Lforce = F0' * force # local distributed load components
        forc_f_s[1] = Lforce[1] * L0 / 2
        forc_f_s[2] = Lforce[2] * L0 / 2
        forc_f_s[3] = Lforce[3] * L0 / 2
        forc_f_s[4] = 0
        forc_f_s[5] = -Lforce[3] * L0^2 / 12
        forc_f_s[6] = +Lforce[2] * L0^2 / 12
        forc_f_s[7] = Lforce[1] * L0 / 2
        forc_f_s[8] = Lforce[2] * L0 / 2
        forc_f_s[9] = Lforce[3] * L0 / 2
        forc_f_s[10] = 0
        forc_f_s[11] = +Lforce[3] * L0^2 / 12
        forc_f_s[12] = -Lforce[2] * L0^2 / 12
        mul!(forc_f_m, Transpose(Secc), forc_f_s)
        mul!(forc_g, Te, forc_f_m)
        gatherdofnums!(dchi, dofnums, fes.conn[i]) # degrees of freedom
        assemble!(assembler, forc_g, dofnums)
    end # Loop over elements
    return makevector!(assembler)
end

function distribloads_global(
    self::FEMMLinBeam,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
    fi,
) where {T<:Number,TI<:Number}
    assembler = SysvecAssembler()
    return distribloads_global(self, assembler, geom0, u1, Rfield1, dchi, fi)
end


"""
    inspectintegpoints(
        self::FEMM,
        geom::NodalField{FT},
        felist::AbstractVector{IT},
        inspector::F,
        idat,
        quantity = :Cauchy;
        context...,
    ) where {FEMM<:AbstractFEMM, FT, IT, F<:Function}

Inspect integration points.
"""
function inspectintegpoints(
    self::FEMM,
    geom0::NodalField{FT},
    dchi::NodalField{FT},
    dT::NodalField{FT},
    felist::AbstractVector{IT},
    inspector::F,
    idat,
    quantity = :localforces;
    context...,
) where {FEMM<:FEMMLinBeam,FT,IT,F<:Function}
    fes = self.integdomain.fes
    ecoords0, dofnums = self._ecoords0, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    elmat, elmatTe, Secc = self._elmat, self._elmatTe, self._elmato
    disp_g, disp_f_m, forces = self._elvec, self._elvecf, self._LF
    mforces = deepcopy(forces)
    disp_f_s = deepcopy(disp_f_m)
    aN, dN, DN = self._aN, self._dN, self._DN
    eccentricities = self.eccentricities
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)::Float64
    A, I1, I2, I3, J, A2s, A3s, x1x2_vector, dimensions = properties(fes)
    # Loop over  all the elements and all the quadrature points within them
    for ilist = eachindex(felist) # Loop over elements
        i = felist[ilist]
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i])
        gathervalues_asvec!(dchi, disp_g, fes.conn[i])
        L0, F0 = initial_local_frame!(F0, ecoords0, x1x2_vector[i])
        _eccentricitytransformation(Secc, eccentricities[i, :]...)
        _transfmat!(Te, F0)
        elmat .= zero(eltype(elmat))
        local_stiffness!(elmat, E, G, A[i], I2[i], I3[i], J[i], A2s[i], A3s[i], L0, aN, DN)
        # First transform from the global coordinates into master local coordinates
        mul!(disp_f_m, Transpose(Te), disp_g) # master displ in f-frame
        # Then transform from master to slave (U_s = Secc * U_m)
        mul!(disp_f_s, Secc, disp_f_m) # slave displ in f-frame
        # Finally compute the forces acting on the slave nodes
        mul!(forces, elmat, disp_f_s)
        idat = inspector(idat, i, fes.conn[i], ecoords0, forces, nothing)
    end # Loop over elements
    return idat # return the updated inspector data
end

end # module
