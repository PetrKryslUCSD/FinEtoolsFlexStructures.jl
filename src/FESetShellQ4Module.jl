module FESetShellQ4Module

using FinEtools
using FinEtools.MeshQuadrilateralModule
import FinEtools.FESetModule: cat, subset, nodesperelem
using FinEtools.MatrixUtilityModule: complete_lt! 
using LinearAlgebra: norm, Transpose, mul!

struct FESetShellQ4 <: AbstractFESet2Manifold{4}
end

function local_frame!(fes::FESetShellQ4, F0, J0)
    # This is the tangent to the coordinate curve 1
    a = @view J0[:, 1]
    L0 = norm(a);
    F0[:,1] = a/L0;
    # This is the tangent to the coordinate curve 2
    b = @view J0[:, 2]
    #     F0(:,3)=skewmat(F0(:,1))*b;
    # Now compute the normal
    F0[:, 3] .= (-F0[3,1]*b[2]+F0[2,1]*b[3],
                  F0[3,1]*b[1]-F0[1,1]*b[3],
                 -F0[2,1]*b[1]+F0[1,1]*b[2]);
    F0[:, 3] /= norm(@view F0[:, 3]);
    #     F0(:,2)=skewmat(F0(:,3))*F0(:,1);
    F0[:, 2] .= (-F0[3,3]*F0[2,1]+F0[2,3]*F0[3,1],
                  F0[3,3]*F0[1,1]-F0[1,3]*F0[3,1],
                 -F0[2,3]*F0[1,1]+F0[1,3]*F0[2,1]);
    return F0
end

end # module
