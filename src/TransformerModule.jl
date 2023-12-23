module TransformerModule

using LinearAlgebra: norm, Transpose, mul!
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict

"""
    TransformerQtEQ

TransformerQtEQ of element matrices.

A callable object: computes `Q^T * E * Q`, where `E` the element stiffness
matrix, and `Q` is the transformation matrix. Both are assumed to be square.

Buffers the intermediate result. Hence no allocation is incurred.
"""
struct TransformerQtEQ{T}
    _buff::Matrix{T}
    function TransformerQtEQ(Q::Matrix{T}) where {T}
        _buff = fill(zero(eltype(Q)), size(Q)...)
        return new{T}(_buff)
    end
end

"""
    (o::TransformerQtEQ)(E, Q)

Perform the transformation on the matrix `E` with the transformation matrix
`Q`: `Ebar = Q^T * E * Q`.
"""
(o::TransformerQtEQ)(E, Q) =
    let
        @assert size(o._buff) == size(Q)
        @assert size(E) == size(Q)
        mul!(o._buff, E, Q)
        mul!(E, Transpose(Q), o._buff)
        return E
    end


"""
    Layup2ElementAngle{T}

A callable object: computes the cosine and sine of the angle between the element frame and the layup reference direction.
"""
struct Layup2ElementAngle{T}
    M::Matrix{T}
    M2D::Matrix{T}
end

function Layup2ElementAngle(T = FFlt)
    Layup2ElementAngle(fill(zero(T), 3, 3), fill(zero(T), 2, 2))
end

function (o::Layup2ElementAngle)(E_G, lcsmat)
    mul!(o.M, E_G', lcsmat)
    o.M2D[1, 1] = o.M[1, 1]
    o.M2D[1, 2] = o.M[1, 2]
    o.M2D[2, 1] = o.M[2, 1]
    o.M2D[2, 2] = o.M[2, 2]
    o.M2D[:, 1] ./= norm(@view o.M2D[:, 1])
    o.M2D[:, 2] ./= norm(@view o.M2D[:, 2])
    m = (o.M2D[1, 1] + o.M2D[2, 2]) / 2
    n = (o.M2D[1, 2] - o.M2D[2, 1]) / 2
    sn = (n >= 0.0 ? +1.0 : -1.0)
    n = sn * sqrt(1 - m^2)
    return m, n
end

end # module

