module TransformerModule

using LinearAlgebra: norm, Transpose, mul!
using FinEtools

"""
    Transformer

Transformer of element matrices.

A callable object: computes Q^T * E * T.
Buffers the intermediate result. Hence no allocation is incurred.
"""
struct Transformer
    _buff::FFltMat
    function Transformer(T)
        _buff = fill(0.0, size(T)...); 
        return new(_buff)
    end
end

"""
    (o::Transformer)(elmat, T)

Perform the transformation on the matrix `elmat` with the transformation matrix
`T`.
"""
(o::Transformer)(elmat, T) = begin
    @assert size(o._buff) == size(T)
    @assert size(elmat) == size(T)
    mul!(o._buff, elmat, T)
    mul!(elmat, Transpose(T), o._buff)
    return elmat
end

end # module
