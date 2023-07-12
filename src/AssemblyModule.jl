"""
    AssemblyModule

Module for assemblers  of system matrices and vectors.
"""
module AssemblyModule

__precompile__(true)

using FinEtools
import FinEtools.AssemblyModule: startassembly!, assemble!, makematrix!
using SparseMatricesCSR
using  LinearAlgebra: diag

"""
    SysmatAssemblerSparseCSRSymm{T<:Number} <: AbstractSysmatAssembler

Assembler for a **symmetric square** matrix  assembled from symmetric square
matrices. The matrix is stored in the Compressed Sparse Row (CSR) format.
"""
mutable struct SysmatAssemblerSparseCSRSymm{IT, T, MBT<:AbstractVector{T}, IJT<:AbstractVector{IT}} <: AbstractSysmatAssembler
    buffer_length::IT
    matbuffer::MBT
    rowbuffer::IJT
    colbuffer::IJT
    buffer_pointer::IT
    row_nalldofs::IT
    col_nalldofs::IT
    nomatrixresult::Bool
    force_init::Bool
end

"""
    SysmatAssemblerSparseCSRSymm(zero::T=0.0) where {T<:Number}

Construct blank system matrix assembler for symmetric CSR matrices. The matrix
entries are of type `T`.

# Example

This is how a symmetric sparse matrix is assembled from two square dense matrices.
```
	a = SysmatAssemblerSparseCSRSymm(0.0)                                                        
	startassembly!(a, 5, 5, 3, 7, 7)    
	m = [0.24406   0.599773    0.833404  0.0420141                                             
		0.786024  0.00206713  0.995379  0.780298                                              
		0.845816  0.198459    0.355149  0.224996]                                     
	assemble!(a, m'*m, [5 2 1 4], [5 2 1 4])        
	m = [0.146618  0.53471   0.614342    0.737833                                              
		 0.479719  0.41354   0.00760941  0.836455                                              
		 0.254868  0.476189  0.460794    0.00919633                                            
		 0.159064  0.261821  0.317078    0.77646                                               
		 0.643538  0.429817  0.59788     0.958909]                                   
	assemble!(a, m'*m, [2 3 1 5], [2 3 1 5])                                        
	A = makematrix!(a) 
```
# See also
[]
"""
function SysmatAssemblerSparseCSRSymm(z::T, nomatrixresult = false) where {T}
    return SysmatAssemblerSparseCSRSymm(0, T[z], Int[0], Int[0], 0, 0, 0, nomatrixresult, false)
end

function SysmatAssemblerSparseCSRSymm()
    return SysmatAssemblerSparseCSRSymm(zero(Float64))
end

"""
    startassembly!(self::SysmatAssemblerSparseCSRSymm{T},
      elem_mat_dim::FInt, ignore1::FInt, elem_mat_nmatrices::FInt,
      ndofs::FInt, ignore2::FInt) where {T<:Number}

Start the assembly of a symmetric square CSR global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.
- `elem_mat_nrows`= number of rows in typical element matrix,
- `elem_mat_ncols`= number of columns in a typical element matrix,
- `elem_mat_nmatrices`= number of element matrices,
- `ndofs_row`= Total number of equations in the row direction,
- `ndofs_col`= Total number of equations in the column direction.

The values stored in the buffers are initially undefined!
"""
function startassembly!(self::SysmatAssemblerSparseCSRSymm{T},
    expected_ntriples::IT,
    row_nalldofs::IT,
    col_nalldofs::IT;
    force_init = false
    ) where {T, IT<:Integer}
    # Only resize the buffers if the pointer is less than 1; otherwise the
    # buffers are already initialized and in use.
    if self.buffer_pointer < 1
        self.buffer_length = expected_ntriples
        resize!(self.rowbuffer, self.buffer_length)
        resize!(self.colbuffer, self.buffer_length)
        resize!(self.matbuffer, self.buffer_length)
        self.buffer_pointer = 1
        self.row_nalldofs = row_nalldofs
        self.col_nalldofs = col_nalldofs
    end
    # Leave the buffers uninitialized, unless the user requests otherwise
    self.force_init = force_init
    if self.force_init
        self.rowbuffer .= 1
        self.colbuffer .= 1
        self.matbuffer .= zero(T)
    end
    return self
end

"""
    assemble!(self::SysmatAssemblerSparseCSRSymm{T}, mat::FMat{T},  dofnums::FIntVec, ignore::FIntVec) where {T<:Number}

Assemble a square symmetric matrix.

`dofnums` are the row degree of freedom numbers, the column degree of freedom
number input is ignored (the row and column numbers are assumed to be the same).
"""
function assemble!(
    self::SysmatAssemblerSparseCSRSymm,
    mat::MBT,
    dofnums_row::CIT,
    dofnums_col::CIT,
) where {MBT, CIT}
    # Assembly of a rectangular matrix.
    # The method assembles a rectangular matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    p = self.buffer_pointer
    if p + ncolumns * nrows >= self.buffer_length
        self = _resize(self, ncolumns * nrows * 1000)
    end
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j in 1:ncolumns
        dj = dofnums_col[j]
        dj < 1 && error("Column degree of freedom < 1")
        dj > self.col_nalldofs && error("Column degree of freedom > size")
        for i in 1:nrows
            di = dofnums_row[i]
            di < 1 && error("Row degree of freedom < 1")
            di > self.row_nalldofs && error("Row degree of freedom > size")
            self.matbuffer[p] = mat[i, j] # serialized matrix
            self.rowbuffer[p] = dofnums_row[i]
            self.colbuffer[p] = dofnums_col[j]
            p = p + 1
        end
    end
    self.buffer_pointer = p
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparseCSRSymm)

Make a sparse symmetric square matrix.
"""
function makematrix!(self::SysmatAssemblerSparseCSRSymm)
    @assert length(self.rowbuffer) >= self.buffer_pointer - 1
    @assert length(self.colbuffer) >= self.buffer_pointer - 1
    # We have the option of retaining the assembled results, but not
    # constructing the sparse matrix.
    if self.nomatrixresult
        # Dummy (zero) sparse matrix is returned. The entire result of the
        # assembly is preserved in the assembler buffers. The ends of the
        # buffers are filled with legal (ignorable) values.
        self.rowbuffer[self.buffer_pointer:end] .= 1
        self.colbuffer[self.buffer_pointer:end] .= 1
        self.matbuffer[self.buffer_pointer:end] .= 0.0
        return spzeros(self.row_nalldofs, self.col_nalldofs)
    end
    # The sparse matrix is constructed and returned. The  buffers used for
    # the assembly are cleared.
    S = SparseMatricesCSR.sparsecsr(
        view(self.rowbuffer, 1:self.buffer_pointer-1),
        view(self.colbuffer, 1:self.buffer_pointer-1),
        view(self.matbuffer, 1:self.buffer_pointer-1),
        self.row_nalldofs,
        self.col_nalldofs,
        )
    # Get ready for more assembling
    self.buffer_pointer = 1
    # Construct the blocks of the matrix
    return S
end


# """
#     SysmatAssemblerVecOfMatSymm{T<:Number} <: AbstractSysmatAssembler

# Assembler for a vector of matrices.

# The global matrix is represented by a collection of constituent matrices, which
# in this case are not assembled, but rather kept in a vector. The global matrix
# is assumed to be symmetric.
# """
# mutable struct SysmatAssemblerVecOfMatSymm{T<:Number} <: AbstractSysmatAssembler
#     # Type for assembling of a sparse global matrix from elementwise matrices.
#     elem_mat_dim::FInt
#     buffer_length::FInt;
#     buffer::Vector{Tuple{Vector{FInt}, Matrix{T}}};
#     buffer_pointer::FInt;
#     row_nalldofs::IT
#     col_nalldofs::IT
#     nomatrixresult::Bool
#     force_init::Bool
# end

# """
#     SysmatAssemblerVecOfMatSymm(zero::T=0.0) where {T<:Number}

# Construct blank system matrix assembler for vector of matrices. The matrix
# entries have elements of type `T`.

# # Example

# This is how a symmetric sparse matrix is assembled from two square dense matrices.
# ```
#     a = SysmatAssemblerVecOfMatSymm(0.0)
#     startassembly!(a, 5, 5, 3, 7, 7)
#     m = [0.24406   0.599773    0.833404  0.0420141
#         0.786024  0.00206713  0.995379  0.780298
#         0.845816  0.198459    0.355149  0.224996]
#     assemble!(a, m'*m, [5 2 1 4], [5 2 1 4])
#     m = [0.146618  0.53471   0.614342    0.737833
#          0.479719  0.41354   0.00760941  0.836455
#          0.254868  0.476189  0.460794    0.00919633
#          0.159064  0.261821  0.317078    0.77646
#          0.643538  0.429817  0.59788     0.958909]
#     assemble!(a, m'*m, [2 3 1 5], [2 3 1 5])
#     A = makematrix!(a)
# ```
# # See also
# []
# """
# function SysmatAssemblerVecOfMatSymm(zero::T=0.0) where {T<:Number}
#     return SysmatAssemblerVecOfMatSymm{T}(0,0,[],0,0)
# end

# """
#     startassembly!(self::SysmatAssemblerVecOfMatSymm{T}, ignore0::FInt, ignore1::FInt, elem_mat_nmatrices::FInt, ndofs::FInt, ignore2::FInt) where {T<:Number}

# Start the assembly of a symmetric square CSR global matrix.

# The method makes buffers for matrix assembly. It must be called before
# the first call to the method `assemble!`.
# - `elem_mat_dim`= dimension of the element matrix,
# - `ignore1`= ignored,
# - `elem_mat_nmatrices`= number of element matrices,
# - `ndofs`= Total number of equations in the row and column direction.

# The values stored in the buffers are initially undefined!
# """
# function startassembly!(self::SysmatAssemblerVecOfMatSymm{T}, elem_mat_dim::FInt, ignore1::FInt, elem_mat_nmatrices::FInt, ndofs::FInt, ignore2::FInt) where {T<:Number}
#     self.elem_mat_dim = elem_mat_dim
#     self.buffer_length = elem_mat_nmatrices;
#     self.buffer = Array{Tuple{Vector{FInt}, Matrix{T}}, 1}(undef, self.buffer_length);
#     self.buffer_pointer = 1;
#     self.ndofs = ndofs;
#     return self
# end

# """
#     assemble!(self::SysmatAssemblerVecOfMatSymm{T}, mat::FMat{T}, dofnums::FIntVec, ignore::FIntVec) where {T<:Number}

# Assemble a square symmetric matrix.

# `dofnums` are the row degree of freedom numbers, the column degree of freedom
# number input is ignored (the row and column numbers are assumed to be the same).
# """
# function assemble!(self::SysmatAssemblerVecOfMatSymm{T}, mat::FMat{T},  dofnums::Union{FIntVec, FIntMat}, ignore::Union{FIntVec, FIntMat}) where {T<:Number}
#     p = self.buffer_pointer
#     @assert p <= self.buffer_length
#     @assert self.elem_mat_dim == size(mat, 1)
#     @assert self.elem_mat_dim == size(mat, 2)
#     @assert self.elem_mat_dim == length(dofnums)
#     d, m = deepcopy(vec(dofnums)), deepcopy(mat)
#     fnz = 1
#     @inbounds for j in 1:length(d)
#         if d[j] == 0
#             d[j] = fnz
#             @inbounds for k in 1:length(d)
#                 m[j, k] = 0.0
#             end
#             @inbounds for k in 1:length(d)
#                 m[k, j] = 0.0
#             end
#         else
#             fnz = d[j]
#         end
#     end
#     self.buffer[p] = (d, m)
#     self.buffer_pointer=p+1;
#     return self
# end

# """
#     makematrix!(self::SysmatAssemblerVecOfMatSymm)

# This function is dummy: the contents is preserved in the assembler itself.
# """
# function makematrix!(self::SysmatAssemblerVecOfMatSymm)
#     self.buffer_length = self.buffer_pointer-1
#     return nothing
# end

end
