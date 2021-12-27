module CompositeLayupModule

using FinEtools
using LinearAlgebra: norm, Transpose, mul!
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!

struct CompositeLayupPly{M}
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
    # Both stiffness matrices are expressed on the reference coordinate system
    # for the ply: direction 1 = fiber direction, direction 3 = normal to the
    # ply

    return CompositeLayupPly(name, material, thickness, angle, Dps, Dts)
end

    
function _shell_material_stiffness(material)
    # Compute the two material stiffness matrices: the plane stress matrix, and
    # the transverse shear matrix
    D = fill(0.0, 6, 6)
    t::FFlt, dt::FFlt, loc::FFltMat, label::FInt = 0.0, 0.0, [0.0 0.0 0.0], 0
    tangentmoduli!(material,  D,  t, dt, loc, label)
    Dps = fill(0.0, 3, 3)
    Dps[1:2, 1:2] = D[1:2, 1:2] -  (reshape(D[1:2,3], 2, 1) * reshape(D[3,1:2], 1, 2))/D[3, 3]
    ix=[1, 2, 4];
    for i = 1:3
        Dps[3,i] = Dps[i,3] = D[4, ix[i]];
    end
    Dt = fill(0.0, 2, 2)
    ix=[5, 6];
    for i = 1:2
        Dt[i,i] = D[ix[i], ix[i]];
    end
    return Dps, Dt
end

end # module
