module TransformerModule

using LinearAlgebra: norm, Transpose, mul!
using FinEtools

"""
    QTEQTransformer

QTEQTransformer of element matrices.

A callable object: computes `Q^T * E * Q`, where `E` the element stiffness
matrix, and `Q` is the transformation matrix.

Buffers the intermediate result. Hence no allocation is incurred.
"""
struct QTEQTransformer
    _buff::FFltMat
    function QTEQTransformer(Q)
        _buff = fill(0.0, size(Q)...); 
        return new(_buff)
    end
end

"""
    (o::QTEQTransformer)(E, Q)

Perform the transformation on the matrix `E` with the transformation matrix
`Q`: `Ebar = Q^T * E * Q`.
"""
(o::QTEQTransformer)(E, Q) = begin
    @assert size(o._buff) == size(Q)
    @assert size(E) == size(Q)
    mul!(o._buff, E, Q)
    mul!(E, Transpose(Q), o._buff)
    return E
end



"""
    QEQTTransformer

QEQTTransformer of element matrices.

A callable object: computes `Q * E * Q^T`, where `E` the element stiffness
matrix, and `Q` is the transformation matrix.

Buffers the intermediate result. Hence no allocation is incurred.
"""
struct QEQTTransformer
    _buff::FFltMat
    function QEQTTransformer(Q)
        _buff = fill(0.0, size(Q)...); 
        return new(_buff)
    end
end

"""
    (o::QEQTTransformer)(E, Q)

Perform the transformation on the matrix `E` with the transformation matrix
`Q`: `Ebar = Q^T * E * Q`.
"""
(o::QEQTTransformer)(E, Q) = begin
    @assert size(o._buff) == size(Q)
    @assert size(E) == size(Q)
    mul!(o._buff, E, Transpose(Q))
    mul!(E, Q, o._buff)
    return E
end

end # module
