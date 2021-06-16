module FESetShellQ4SRIModule

using FinEtools
import FinEtools: cat, subset
using FinEtools.MatrixUtilityModule: complete_lt! 
using LinearAlgebra: norm, Transpose, mul!

mutable struct FESetShellQ4SRI <: AbstractFESet2Manifold{4}
    conn::Array{NTuple{4, FInt}, 1};
    label::FIntVec; 
    dimensions::Vector{FFltVec}
end

function FESetShellQ4SRI(conn::FIntMat) 
    N = size(conn, 1)
    self = FESetL2CorotBeam(NTuple{4, FInt}[], FInt[])
    self = fromarray!(self, conn)
    setlabel!(self, 0)
    return self
end

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
