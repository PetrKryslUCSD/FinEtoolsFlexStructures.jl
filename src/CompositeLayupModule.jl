module CompositeLayupModule

using FinEtools
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.TransformerModule: QEQTTransformer

function cartesian_csys(axes)
    function   cartesian!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) 
        csmatout[:] .= 0.0
        for j in 1:3
            aj = abs(axes[j]); sj = sign(axes[j])
            csmatout[aj, j] = sj * 1.0
        end
        return csmatout
    end
    return  CSys(3, 3, cartesian!)
end
     
abstract type AbstractPly end

struct Ply{M} <: AbstractPly
    name::String # name of the ply
    material::M # material of the ply
    thickness::FFlt # thickness of the ply
    angle::FFlt # angle in degrees
    _Dps::FFltMat # plane stress stiffness matrix
    _Dts::FFltMat # transverse shear stiffness matrix
end

function Ply(name, material::M, thickness, angle) where {M}
    # First we extract the full three dimensional stiffness matrix of the
    # material
    D = fill(0.0, 6, 6)
    # What if the material changes composition from point to point? What do we
    # do now?
    t::FFlt, dt::FFlt, loc::FFltMat, label::FInt = 0.0, 0.0, [0.0 0.0 0.0], 0
    tangentmoduli!(material,  D,  t, dt, loc, label)
    # Next we reduce the 3d stiffness to plane stress
    # Both stiffness matrices are expressed on the reference coordinate system
    # for the ply: direction 1 = fiber direction, direction 3 = normal to the
    # ply
    Dps = fill(0.0, 3, 3)
    Dps[1:2, 1:2] = D[1:2, 1:2] -  (reshape(D[1:2,3], 2, 1) * reshape(D[3,1:2], 1, 2))/D[3, 3]
    ix = [1, 2, 4];
    for i = 1:3
        Dps[3,i] = Dps[i,3] = D[4, ix[i]];
    end
    # And we also extract the transverse shear matrix
    Dts = fill(0.0, 2, 2)
    ix = [5, 6];
    for i = 1:2
        Dts[i,i] = D[ix[i], ix[i]];
    end
    return Ply(name, material, thickness, FFlt(angle), Dps, Dts)
end

function lamina_material(E1, E2, nu12, G12, G13, G23)
    rho = 0.0
    return MatDeforElastOrtho(DeforModelRed3D, rho, E1, E2, E2, nu12, 0.0, 0.0, G12, G13, G23, 0.0, 0.0, 0.0)
end

struct CompositeLayup
    name::String
    offset::FFlt # offset of the reference surface from the mid surface, negative when the offset is against the normal
    plies::Vector{AbstractPly} # vector of plies; the first ply is at the bottom of the shell (SNEG surface), the last ply is at the top (SPOS)
    mcsys::CSys # updater of the material orientation matrix
end

function CompositeLayup(name, plies, mcsys)
    offset = 0.0
    return CompositeLayup(name, offset, plies, mcsys)
end

struct CompositeLayups
    layup_list::Vector{CompositeLayup}
end

function laminate_stiffnesses!(cl::CompositeLayup, A, B, C)
    # Aij coefficients represent in-plane stiffness of the laminate, the Cij
    # coefficients represent bending stiffness, the Bij represent
    # bending-extension coupling, and the Hij represent intralaminar shear
    # stiffness.
    A .= zero(eltype(A))
    Dps = deepcopy(A)
    T = deepcopy(A)
    tf = QEQTTransformer(Dps)
    # Transform into the composite layup coordinate system.
    layup_thickness = sum(p.thickness for p in cl.plies)
    zs = -layup_thickness/2 - cl.offset
    for p in cl.plies
        ze = zs + p.thickness
        plane_stress_T_matrix!(T, -p.angle/180*pi)
        # Transform the plane stress matrix into the layup coordinates
        @. Dps = p._Dps
        Dps = tf(Dps, T)
        # Compute the in-plane stiffness
        @. A += (ze - zs) * Dps
        # Compute the extension-bending coupling stiffness
        @. B += (ze^2 - zs^2)/2 * Dps
        # Compute the bending stiffness
        @. C += (ze^3 - zs^3)/3 * Dps
        zs += p.thickness
    end
    return A, B, C
end

function laminate_transverse_stiffness!(cl::CompositeLayup, H)
    # Hij represent intralaminar shear stiffness.
    H .= zero(eltype(H))
    Dts = deepcopy(H)
    T = deepcopy(H)
    tf = QEQTTransformer(Dts)
    # Transform into the composite layup coordinate system.
    layup_thickness = sum(p.thickness for p in cl.plies)
    zs = -layup_thickness/2 - cl.offset
    for p in cl.plies
        ze = zs + p.thickness
        transverse_shear_T_matrix!(T, p.angle/180*pi)
        # Transform the plane stress matrix into the layup coordinates
        @. Dts = p._Dts
        Dts = tf(Dts, T)
        # Compute the transverse shear stiffness The correction factor
        # accounting for the parable distribution of the shear stress is due to
        # Vinson, Sierakowski, The Behavior of Structures Composed of Composite
        # Materials, 2008
        @. H += 5/4*(ze - zs - 4/3*(ze^3-zs^3)/layup_thickness^2) * Dts
        zs += p.thickness
    end
    return H
end

"""
    plane_stress_T_matrix!( Tinvm::Array{T, 2},
     angle) where {T}

Compute the transformation matrix between strain tensor components on the layup
coordinate system into the ply coordinate system.

`angle` = angle between the first basis vector of the layup coordinate system
    and the first basis vector of the ply coordinate system
"""
function  plane_stress_T_matrix!(Tm::Array{T, 2}, angle) where {T}
    m=cos(angle); 
    n=sin(angle); 
    Tm[1, 1] =  (m^2);  Tm[1, 2] = (n^2);   Tm[1, 3] = (2*m*n)
    Tm[2, 1] = (n^2);   Tm[2, 2] = (m^2);   Tm[2, 3] = (-2*m*n)
    Tm[3, 1] = (-m*n);  Tm[3, 2] = (m*n);   Tm[3, 3] = (m*m-n*n)
    return Tm
end

"""
    plane_stress_Tinv_matrix!(Tinvm::Array{T, 2}, angle) where {T}

Compute the transformation matrix between strain tensor components on the ply
coordinate system into the layout coordinate system.

`angle` = angle between the first basis vector of the layup coordinate system
    and the first basis vector of the ply coordinate system
"""
function  plane_stress_Tinv_matrix!(Tinvm::Array{T, 2}, angle) where {T}
    m=cos(angle); 
    n=sin(angle); 
    Tinvm[1, 1] = (m^2); Tinvm[1, 2] = (n^2);  Tinvm[1, 3] = (-2*m*n)
    Tinvm[2, 1] = (n^2); Tinvm[2, 2] = (m^2);  Tinvm[2, 3] = (2*m*n)
    Tinvm[3, 1] = (m*n); Tinvm[3, 2] = (-m*n); Tinvm[3, 3] = (m*m-n*n)
    return Tinvm
end

"""
    plane_stress_T_matrix_eng!( Tinvm::Array{T, 2},
     angle) where {T}

Compute the transformation matrix between strain engineering components on the
layup coordinate system into the ply coordinate system.

`angle` = angle between the first basis vector of the layup coordinate system
    and the first basis vector of the ply coordinate system

The matrix of transformation is `Tme = R * Tm / R`, where `Tm`
is the transformation matrix in tensor components, and `R` is the Reuter matrix,
`R = [1 0 0; 0 1 0; 0 0 2]`.
"""
function  plane_stress_T_matrix_eng!(Tm::Array{T, 2}, angle) where {T}
    m=cos(angle); 
    n=sin(angle); 
    Tm[1, 1] =  (m^2);  Tm[1, 2] = (n^2);   Tm[1, 3] = (m*n)
    Tm[2, 1] = (n^2);   Tm[2, 2] = (m^2);   Tm[2, 3] = (-m*n)
    Tm[3, 1] = (-2*m*n);  Tm[3, 2] = (2*m*n);   Tm[3, 3] = (m*m-n*n)
    return Tm
end

"""
    plane_stress_Tinv_matrix!(Tinvm::Array{T, 2}, angle) where {T}

Compute the transformation matrix between strain engineering components on the
ply coordinate system into the layout coordinate system.

`angle` = angle between the first basis vector of the layup coordinate system
    and the first basis vector of the ply coordinate system

The matrix of transformation is `Tme = R * Tm / R`, where `Tm`
is the transformation matrix in tensor components, and `R` is the Reuter matrix,
`R = [1 0 0; 0 1 0; 0 0 2]`.
"""
function  plane_stress_Tinv_matrix_eng!(Tinvm::Array{T, 2}, angle) where {T}
    m=cos(angle); 
    n=sin(angle); 
    Tinvm[1, 1] = (m^2); Tinvm[1, 2] = (n^2);  Tinvm[1, 3] = (-m*n)
    Tinvm[2, 1] = (n^2); Tinvm[2, 2] = (m^2);  Tinvm[2, 3] = (m*n)
    Tinvm[3, 1] = (2*m*n); Tinvm[3, 2] = (-2*m*n); Tinvm[3, 3] = (m*m-n*n)
    return Tinvm
end

function  transverse_shear_T_matrix!(Tm::Array{T, 2}, angle) where {T}
    m=cos(angle); 
    n=sin(angle); 
    # Barbero, Introduction, a in equation 5.7
    Tm[1, 1] =  m;  Tm[1, 2] = n;   
    Tm[2, 1] = -n;  Tm[2, 2] = m;  
    return Tm
end

end # module