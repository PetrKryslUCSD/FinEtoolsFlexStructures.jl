module CompositeLayupModule

using FinEtools
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed2DStress

abstract type AbstractPly end

struct CompositeLayupPly{M} <: AbstractPly
    name::String
    material::M
    thickness::FFlt
    angle::FFlt
    _Dps::FFltMat # plane stress stiffness matrix
    _Dts::FFltMat # transverse shear stiffness matrix
end

function CompositeLayupPly(name, material::M, thickness, angle) where {M}
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
    
    return CompositeLayupPly(name, material, thickness, angle, Dps, Dts)
end

struct CompositeLayup
    name::String
    plies::Vector{AbstractPly}
end

function plane_stress_stiffness!(cl::CompositeLayup, D)
    # Transform into the composite layup coordinate system.
end

"""
    plane_stress_T_matrix!(::Type{DeforModelRed2DStress},  Tinvm::Array{T, 2},
     angle) where {T}

Compute the transformation matrix between strain tensor components on the layup
coordinate system into the ply coordinate system.

`angle` = angle between the first basis vector of the layup coordinate system
    and the first basis vector of the ply coordinate system
"""
function  plane_stress_T_matrix!(::Type{DeforModelRed2DStress}, Tm::Array{T, 2}, angle) where {T}
    m=cos(angle); 
    n=sin(angle); 
    Tm[1, 1] =  (m^2);  Tm[1, 2] = (n^2);   Tm[1, 3] = (2*m*n)
    Tm[2, 1] = (n^2);   Tm[2, 2] = (m^2);   Tm[2, 3] = (-2*m*n)
    Tm[3, 1] = (-m*n);  Tm[3, 2] = (m*n);   Tm[3, 3] = (m*m-n*n)
    return Tm
end

"""
    plane_stress_Tinv_matrix!(::Type{DeforModelRed2DStress}, Tinvm::Array{T, 2}, angle) where {T}

Compute the transformation matrix between strain tensor components on the ply
coordinate system into the layout coordinate system.

`angle` = angle between the first basis vector of the layup coordinate system
    and the first basis vector of the ply coordinate system
"""
function  plane_stress_Tinv_matrix!(::Type{DeforModelRed2DStress}, Tinvm::Array{T, 2}, angle) where {T}
    m=cos(angle); 
    n=sin(angle); 
    Tinvm[1, 1] = (m^2); Tinvm[1, 2] = (n^2);  Tinvm[1, 3] = (-2*m*n)
    Tinvm[2, 1] = (n^2); Tinvm[2, 2] = (m^2);  Tinvm[2, 3] = (2*m*n)
    Tinvm[3, 1] = (m*n); Tinvm[3, 2] = (-m*n); Tinvm[3, 3] = (m*m-n*n)
    return Tinvm
end

"""
    plane_stress_T_matrix_eng!(::Type{DeforModelRed2DStress},  Tinvm::Array{T, 2},
     angle) where {T}

Compute the transformation matrix between strain engineering components on the
layup coordinate system into the ply coordinate system.

`angle` = angle between the first basis vector of the layup coordinate system
    and the first basis vector of the ply coordinate system

The matrix of transformation is `Tme = R * Tm / R`, where `Tm`
is the transformation matrix in tensor components, and `R` is the Reuter matrix,
`R = [1 0 0; 0 1 0; 0 0 2]`.
"""
function  plane_stress_T_matrix_eng!(::Type{DeforModelRed2DStress}, Tm::Array{T, 2}, angle) where {T}
    m=cos(angle); 
    n=sin(angle); 
    Tm[1, 1] =  (m^2);  Tm[1, 2] = (n^2);   Tm[1, 3] = (m*n)
    Tm[2, 1] = (n^2);   Tm[2, 2] = (m^2);   Tm[2, 3] = (-m*n)
    Tm[3, 1] = (-2*m*n);  Tm[3, 2] = (2*m*n);   Tm[3, 3] = (m*m-n*n)
    return Tm
end

"""
    plane_stress_Tinv_matrix!(::Type{DeforModelRed2DStress}, Tinvm::Array{T, 2}, angle) where {T}

Compute the transformation matrix between strain engineering components on the
ply coordinate system into the layout coordinate system.

`angle` = angle between the first basis vector of the layup coordinate system
    and the first basis vector of the ply coordinate system

The matrix of transformation is `Tme = R * Tm / R`, where `Tm`
is the transformation matrix in tensor components, and `R` is the Reuter matrix,
`R = [1 0 0; 0 1 0; 0 0 2]`.
"""
function  plane_stress_Tinv_matrix_eng!(::Type{DeforModelRed2DStress}, Tinvm::Array{T, 2}, angle) where {T}
    m=cos(angle); 
    n=sin(angle); 
    Tinvm[1, 1] = (m^2); Tinvm[1, 2] = (n^2);  Tinvm[1, 3] = (-m*n)
    Tinvm[2, 1] = (n^2); Tinvm[2, 2] = (m^2);  Tinvm[2, 3] = (m*n)
    Tinvm[3, 1] = (2*m*n); Tinvm[3, 2] = (-2*m*n); Tinvm[3, 3] = (m*m-n*n)
    return Tinvm
end

end # module
