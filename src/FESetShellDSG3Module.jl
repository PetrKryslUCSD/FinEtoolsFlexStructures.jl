module FESetShellDSG3Module

using FinEtools
using FinEtools.MeshQuadrilateralModule
import FinEtools.FESetModule: cat, subset, nodesperelem
using FinEtools.MatrixUtilityModule: complete_lt! 
using LinearAlgebra: norm, Transpose, mul!

mutable struct FESetShellDSG3 <: AbstractFESet2Manifold{3}
    conn::Array{NTuple{3, FInt}, 1};
    label::FIntVec; 
    thickness::FFlt
end

function FESetShellDSG3(conn::FIntMat; thickness = zero(FFlt)) 
    @assert size(conn, 2) == 3
    self = FESetShellDSG3(NTuple{3, FInt}[], FInt[], thickness)
    self = fromarray!(self, conn)
    setlabel!(self, 0)
    return self
end

nodesperelem(fes::FESetShellDSG3) = 3

function local_frame!(F0, J0)
    # This is the tangent to the coordinate curve 1
    a = @view J0[:, 1]
    L0 = norm(a);
    F0[:,1] = a/L0;
    # This is the tangent to the coordinate curve 2
    b = J0[:, 2]
    #     F0(:,3)=skewmat(F0(:,1))*b;
    b = @view J0[:, 2]
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
