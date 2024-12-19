"""
Module for working with composite layups.
"""
module CompositeLayupModule

using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.TransformerModule: TransformerQtEQ

"""
    cartesian_csys(axes)

Create a material Cartesian coordinate system.

- `axes` = tuple of signed labels of the axes. For instance, `(1, 2, 3)` creates
  a coordinate system identical to the global cartesian coordinate system. `
  (2, 1, -3)` creates a coordinate system so that the first material basis
  vector is along the second global basis vector, the second material basis
  vector is along the first global basis vector, and the third material basis
  vector is opposite to the third global basis vector.
"""
function cartesian_csys(axes)
    function cartesian!(
        csmatout::FFltMat,
        XYZ::FFltMat,
        tangents::FFltMat,
        feid::FInt,
        qpid::FInt,
    )
        csmatout[:] .= 0.0
        for j = 1:3
            aj = abs(axes[j])
            sj = sign(axes[j])
            csmatout[aj, j] = sj * 1.0
        end
        return csmatout
    end
    return CSys(3, 3, cartesian!)
end

abstract type AbstractPly end

"""
    Ply{M} <: AbstractPly

Type of a ply.
"""
struct Ply{M} <: AbstractPly
    name::String # name of the ply
    material::M # material of the ply
    thickness::FFlt # thickness of the ply
    angle::FFlt # angle in degrees
    _Dps::FFltMat # plane stress stiffness matrix
    _Dts::FFltMat # transverse shear stiffness matrix
end

"""
    Ply{M} <: AbstractPly

Create a ply.

Provide name, material of the ply, thickness of the ply, an angle between the
first bases vector of the layup coordinate system and the first direction of
the ply material coordinate system.
"""
function Ply(name, material::M, thickness, angle) where {M}
    # First we extract the full three dimensional stiffness matrix of the
    # material
    D = fill(0.0, 6, 6)
    # What if the material changes composition from point to point? What do we
    # do now?
    t::FFlt, dt::FFlt, loc::FFltMat, label::FInt = 0.0, 0.0, [0.0 0.0 0.0], 0
    tangentmoduli!(material, D, t, dt, loc, label)
    # Next we reduce the 3d stiffness to plane stress
    # Both stiffness matrices are expressed on the reference coordinate system
    # for the ply: direction 1 = fiber direction, direction 3 = normal to the
    # ply
    Dps = fill(0.0, 3, 3)
    Dps[1:2, 1:2] =
        D[1:2, 1:2] - (reshape(D[1:2, 3], 2, 1) * reshape(D[3, 1:2], 1, 2)) / D[3, 3]
    ix = [1, 2, 4]
    for i = 1:3
        Dps[3, i] = Dps[i, 3] = D[4, ix[i]]
    end
    # And we also extract the transverse shear matrix
    Dts = fill(0.0, 2, 2)
    ix = [5, 6]
    for i = 1:2
        Dts[i, i] = D[ix[i], ix[i]]
    end
    return Ply(name, material, thickness, FFlt(angle), Dps, Dts)
end

"""
    lamina_material(E1, E2, nu12, G12, G13, G23)

Create a transversely isotropic lamina material with default (zero) mass
density.
"""
function lamina_material(E1, E2, nu12, G12, G13, G23)
    rho = 0.0
    return MatDeforElastOrtho(
        DeforModelRed3D,
        rho,
        E1,
        E2,
        E2,
        nu12,
        nu12,
        0.0,
        G12,
        G13,
        G23,
        0.0,
        0.0,
        0.0,
    )
end

"""
    lamina_material(rho, E1, E2, nu12, G12, G13, G23)

Create a transversely isotropic lamina material.
"""
function lamina_material(rho, E1, E2, nu12, G12, G13, G23)
    return MatDeforElastOrtho(
        DeforModelRed3D,
        rho,
        E1,
        E2,
        E2,
        nu12,
        nu12,
        0.0,
        G12,
        G13,
        G23,
        0.0,
        0.0,
        0.0,
    )
end

"""
    lamina_material(E, nu)

Create an isotropic lamina material with default (zero) mass
density.
"""
function lamina_material(E, nu)
    rho = 0.0
    return MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
end

"""
    lamina_material(rho, E, nu)

Create an isotropic lamina material.
"""
function lamina_material(rho, E, nu)
    return MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
end

"""
    CompositeLayup

Type for composite layup.
"""
struct CompositeLayup
    name::String
    offset::FFlt # offset of the reference surface from the mid surface, negative when the offset is against the normal
    plies::Vector{AbstractPly} # vector of plies; the first ply is at the bottom of the shell (SNEG surface), the last ply is at the top (SPOS)
    csys::CSys # updater of the material orientation matrix
end

"""
    CompositeLayup(name, plies, mcsys)

Create a composite layup.

Provide the name, the array of plies, and the coordinate system that defines the
orientation of the composite layup. The first base spector of this coordinate
system is the reference direction for the layup.
"""
function CompositeLayup(name, plies, mcsys)
    offset = 0.0
    return CompositeLayup(name, offset, plies, mcsys)
end

"""
    thickness(cl::CompositeLayup)

Compute the thickness of the layup (some of the thicknesses of the plies).
"""
function thickness(cl::CompositeLayup)
    return sum(p.thickness for p in cl.plies)
end

"""
    laminate_stiffnesses!(cl::CompositeLayup, A, B, D)

Compute the laminate stiffness matrices, membrane, extension-bending
coupling, and bending.

A, B, and D are the familiar 3x3 matrices: Aij coefficients represent in-plane
stiffness of the laminate, the Dij coefficients represent bending stiffness,
the Bij represent bending-extension coupling.
"""
function laminate_stiffnesses!(cl::CompositeLayup, A, B, D)
    A .= zero(eltype(A))
    B .= zero(eltype(A))
    D .= zero(eltype(A))
    Dps = deepcopy(A)
    Tbar = deepcopy(A)
    tf = TransformerQtEQ(Dps)
    # Transform into the composite layup coordinate system.
    layup_thickness = thickness(cl)
    zs = -layup_thickness / 2 - cl.offset
    for p in cl.plies
        ze = zs + p.thickness
        plane_stress_Tbar_matrix!(Tbar, p.angle / 180 * pi)
        # Transform the plane stress matrix into the layup coordinates
        @. Dps = p._Dps
        Dps = tf(Dps, Tbar)
        # Compute the in-plane stiffness
        @. A += (ze - zs) * Dps
        # Compute the extension-bending coupling stiffness
        @. B += (ze^2 - zs^2) / 2 * Dps
        # Compute the bending stiffness
        @. D += (ze^3 - zs^3) / 3 * Dps
        zs += p.thickness
    end
    return A, B, D
end

"""
    laminate_transverse_stiffness!(cl::CompositeLayup, H)

Compute the laminate transverse stiffness.

Hij represent intralaminar shear stiffness.
"""
function laminate_transverse_stiffness!(cl::CompositeLayup, H)
    H .= zero(eltype(H))
    Dts = deepcopy(H)
    T = deepcopy(H)
    tf = TransformerQtEQ(Dts)
    # Transform into the composite layup coordinate system.
    layup_thickness = sum(p.thickness for p in cl.plies)
    zs = -layup_thickness / 2 - cl.offset
    for p in cl.plies
        ze = zs + p.thickness
        transverse_shear_T_matrix!(T, p.angle / 180 * pi)
        # Transform the plane stress matrix into the layup coordinates
        @. Dts = p._Dts
        Dts = tf(Dts, T)
        # Compute the transverse shear stiffness. The correction factor
        # accounting for the parabolic distribution of the shear stress is due to
        # Vinson, Sierakowski, The Behavior of Structures Composed of Composite
        # Materials, 2008
        @. H += 5 / 4 * (ze - zs - 4 / 3 * (ze^3 - zs^3) / layup_thickness^2) * Dts
        zs += p.thickness
    end
    return H
end

"""
    laminate_inertia!(cl::CompositeLayup)

Compute the laminate inertia.
"""
function laminate_inertia!(cl::CompositeLayup)
    layup_thickness = sum(p.thickness for p in cl.plies)
    zs = -layup_thickness / 2 - cl.offset
    mass_density = 0.0
    moment_of_inertia_density = 0.0
    for p in cl.plies
        rho = massdensity(p.material) # mass density
        ze = zs + p.thickness
        mass_density += (ze - zs) * rho
        moment_of_inertia_density += (ze^3 - zs^3) * rho / 3
        zs += p.thickness
    end
    return mass_density, moment_of_inertia_density
end

"""
    plane_stress_Tbar_matrix!(Tm::Array{T, 2}, angle) where {T}
    plane_stress_Tbar_matrix!(Tm::Array{T, 2}, m, n) where {T}

Compute the transformation matrix of engineering strain components FROM the
LAYUP coordinate system TO the PLY coordinate system.

`angle` = angle (in radians) between the first basis vector of the layup
    coordinate system and the first basis vector of the ply coordinate system
`m`, `n` = cosine and sine of the angle

The nomenclature is from Barbero, Finite element analysis of composite materials
using Abaqus (2013).
"""
function plane_stress_Tbar_matrix!(Tbarm::Array{T,2}, angle) where {T}
    m = cos(angle)
    n = sin(angle)
    return plane_stress_Tbar_matrix!(Tbarm, m, n)
end

function plane_stress_Tbar_matrix!(Tbarm::Array{T,2}, m, n) where {T}
    # We are using here the relation between Tbar and T: Tbar = T^-T
    plane_stress_Tinv_matrix!(Tbarm, m, n)
    # Transpose in place
    temp = Tbarm[2, 1]
    Tbarm[2, 1] = Tbarm[1, 2]
    Tbarm[1, 2] = temp
    temp = Tbarm[3, 1]
    Tbarm[3, 1] = Tbarm[1, 3]
    Tbarm[1, 3] = temp
    temp = Tbarm[3, 2]
    Tbarm[3, 2] = Tbarm[2, 3]
    Tbarm[2, 3] = temp
    return Tbarm
end

"""
    plane_stress_Tinv_matrix!(Tinvm::Array{T, 2}, angle) where {T}
    plane_stress_Tinv_matrix!(Tinvm::Array{T, 2}, m, n) where {T}

Compute the transformation matrix of the stress vector components FROM the
LAYUP coordinate system TO the PLY coordinate system.

The components of the stress vector transform as
```math
\\sigma^\\prime = T^{-1} \\sigma
```
where ``\\sigma`` is the stress vector in the layup coordinate system, and
``sigma^\\prime`` is the stress vector in the ply coordinate system. 

`angle` = angle (in radians) between the first basis vector of the layup
    coordinate system and the first basis vector of the ply coordinate system

The nomenclature is from Barbero, Finite element analysis of composite materials
using Abaqus (2013).
"""
function plane_stress_Tinv_matrix!(Tinvm::Array{T,2}, angle) where {T}
    m = cos(angle)
    n = sin(angle)
    return plane_stress_Tinv_matrix!(Tinvm, m, n)
end

function plane_stress_Tinv_matrix!(Tinvm, m, n)
    Tinvm[1, 1] = (m^2)
    Tinvm[1, 2] = (n^2)
    Tinvm[1, 3] = 2 * (m * n)
    Tinvm[2, 1] = (n^2)
    Tinvm[2, 2] = (m^2)
    Tinvm[2, 3] = -2 * (m * n)
    Tinvm[3, 1] = -(m * n)
    Tinvm[3, 2] = (m * n)
    Tinvm[3, 3] = (m^2 - n^2)
    return Tinvm
end

"""
    plane_stress_T_matrix!(Tm::Array{T, 2}, angle) where {T}

Compute the transformation matrix of stress vector components FROM the PLY
coordinate system TO the LAYUP coordinate system.

Their components of the stress vector transform as
```math
\\sigma = T  \\sigma^\\prime
```
where ``\\sigma`` is the stress vector in the layup coordinate system, and
``\\sigma^\\prime`` is the stress vector in the ply coordinate system. 

# Arguments
`angle` = angle (in radians) between the first basis vector of the layup
    coordinate system and the first basis vector of the ply coordinate system

The nomenclature is from Barbero, Finite element analysis of composite materials
using Abaqus (2013).
"""
function plane_stress_T_matrix!(Tm::Array{T,2}, angle) where {T}
    # We are using here the fact that the inverse rotation is given by the
    # negative of the angle
    m = cos(-angle)
    n = sin(-angle)
    return plane_stress_Tinv_matrix!(Tm, m, n)
end

"""
    transverse_shear_T_matrix!(Tm::Array{T, 2}, angle) where {T}
    transverse_shear_T_matrix!(Tm::Array{T, 2}, m, n) where {T}

Compute the transformation matrix for the transverse shear stresses  FROM the
LAYUP coordinate system TO the PLY coordinate system.

`angle` = angle (in radians) between the first basis vector of the layup
    coordinate system and the first basis vector of the ply coordinate system
`m`, `n` = cosine and sine of the angle
"""
function transverse_shear_T_matrix!(Tm::Array{T,2}, angle) where {T}
    m = cos(angle)
    n = sin(angle)
    return transverse_shear_T_matrix!(Tm, m, n)
end

function transverse_shear_T_matrix!(Tm::Array{T,2}, m, n) where {T}
    # Barbero, Introduction, a^T in equation 5.26
    Tm[1, 1] = m
    Tm[1, 2] = -n
    Tm[2, 1] = +n
    Tm[2, 2] = m
    return Tm
end

end # module
