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
mutable struct SysmatAssemblerSparseCSRSymm{T<:Number} <: AbstractSysmatAssembler
    # Type for assembling of a sparse global matrix from elementwise matrices.
    buffer_length:: FInt;
    matbuffer::Vector{T};
    rowbuffer::Vector{FInt};
    colbuffer::Vector{FInt};
    buffer_pointer:: FInt;
    ndofs:: FInt;
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
function SysmatAssemblerSparseCSRSymm(zero::T=0.0) where {T<:Number}
    return SysmatAssemblerSparseCSRSymm{T}(0,[zero],[0],[0],0,0)
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
function startassembly!(self::SysmatAssemblerSparseCSRSymm{T}, elem_mat_dim::FInt, ignore1::FInt, elem_mat_nmatrices::FInt, ndofs::FInt, ignore2::FInt) where {T<:Number}
    self.buffer_length = elem_mat_nmatrices*elem_mat_dim^2;
    self.rowbuffer = Array{FInt, 1}(undef, self.buffer_length);
    self.colbuffer = Array{FInt, 1}(undef, self.buffer_length);
    self.matbuffer = Array{T, 1}(undef, self.buffer_length);
    self.buffer_pointer = 1;
    self.ndofs = ndofs;
    return self
end

"""
    assemble!(self::SysmatAssemblerSparseCSRSymm{T}, mat::FMat{T},  dofnums::FIntVec, ignore::FIntVec) where {T<:Number}

Assemble a square symmetric matrix.

`dofnums` are the row degree of freedom numbers, the column degree of freedom
number input is ignored (the row and column numbers are assumed to be the same).
"""
function assemble!(self::SysmatAssemblerSparseCSRSymm{T}, mat::FMat{T},  dofnums::Union{FIntVec, FIntMat}, ignore::Union{FIntVec, FIntMat}) where {T<:Number}
    # Assembly of a square symmetric matrix.
    # The method assembles the square symmetric matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows=length(dofnums); ncolumns=nrows;
    p = self.buffer_pointer
    @assert p+ncolumns*nrows <= self.buffer_length+1
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j in 1:ncolumns
        @inbounds for i in 1:nrows
            self.matbuffer[p] = mat[i,j] # serialized matrix
            self.rowbuffer[p] = dofnums[i];
            self.colbuffer[p] = dofnums[j];
            p=p+1
        end
    end
    self.buffer_pointer=p;
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparseCSRSymm)

Make a sparse symmetric square matrix.
"""
function makematrix!(self::SysmatAssemblerSparseCSRSymm)
    # Make a sparse matrix.
    # The method makes a sparse matrix from the assembly buffers.
    @assert length(self.rowbuffer) >= self.buffer_pointer-1
    @assert length(self.colbuffer) >= self.buffer_pointer-1
    # Here we will go through the rows and columns, and whenever the row or
    # the column refer to indexes outside of the limits of the matrix, the
    # corresponding value will be set to 0 and assembled to row and column 1.
    @inbounds for j=1:self.buffer_pointer-1
        r = self.rowbuffer[j]
        if (r > self.ndofs) || (r <= 0)
            self.rowbuffer[j] = 1;
            self.matbuffer[j] = 0.0
        end
        c = self.colbuffer[j]
        if (c > self.ndofs) || (c <= 0)
            self.colbuffer[j] = 1;
            self.matbuffer[j] = 0.0
        end
    end
    # Use the CSR package.
    S = SparseMatricesCSR.sparsecsr(self.rowbuffer[1:self.buffer_pointer-1],
               self.colbuffer[1:self.buffer_pointer-1],
               self.matbuffer[1:self.buffer_pointer-1],
               self.ndofs, self.ndofs); 
    self = SysmatAssemblerSparseCSRSymm(zero(eltype(self.matbuffer))) # get rid of the buffers
    return S
end


"""
    SysmatAssemblerVecOfMatSymm{T<:Number} <: AbstractSysmatAssembler

Assembler for a vector of matrices. 

The global matrix is represented by a collection of constituent matrices, which
in this case are not assembled, but rather kept in a vector. The global matrix
is assumed to be symmetric.
"""
mutable struct SysmatAssemblerVecOfMatSymm{T<:Number} <: AbstractSysmatAssembler
    # Type for assembling of a sparse global matrix from elementwise matrices.
    elem_mat_dim::FInt
    buffer_length::FInt;
    buffer::Vector{Tuple{Vector{FInt}, Matrix{T}}};
    buffer_pointer::FInt;
    ndofs::FInt;
end

"""
    SysmatAssemblerVecOfMatSymm(zero::T=0.0) where {T<:Number}

Construct blank system matrix assembler for vector of matrices. The matrix
entries have elements of type `T`.

# Example

This is how a symmetric sparse matrix is assembled from two square dense matrices.
```
    a = SysmatAssemblerVecOfMatSymm(0.0)                                                        
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
function SysmatAssemblerVecOfMatSymm(zero::T=0.0) where {T<:Number}
    return SysmatAssemblerVecOfMatSymm{T}(0,0,[],0,0)
end

"""
    startassembly!(self::SysmatAssemblerVecOfMatSymm{T}, ignore0::FInt, ignore1::FInt, elem_mat_nmatrices::FInt, ndofs::FInt, ignore2::FInt) where {T<:Number}

Start the assembly of a symmetric square CSR global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.
- `elem_mat_dim`= dimension of the element matrix,
- `ignore1`= ignored,
- `elem_mat_nmatrices`= number of element matrices,
- `ndofs`= Total number of equations in the row and column direction.

The values stored in the buffers are initially undefined!
"""
function startassembly!(self::SysmatAssemblerVecOfMatSymm{T}, elem_mat_dim::FInt, ignore1::FInt, elem_mat_nmatrices::FInt, ndofs::FInt, ignore2::FInt) where {T<:Number}
    self.elem_mat_dim = elem_mat_dim
    self.buffer_length = elem_mat_nmatrices;
    self.buffer = Array{Tuple{Vector{FInt}, Matrix{T}}, 1}(undef, self.buffer_length);
    self.buffer_pointer = 1;
    self.ndofs = ndofs;
    return self
end

"""
    assemble!(self::SysmatAssemblerVecOfMatSymm{T}, mat::FMat{T}, dofnums::FIntVec, ignore::FIntVec) where {T<:Number}

Assemble a square symmetric matrix.

`dofnums` are the row degree of freedom numbers, the column degree of freedom
number input is ignored (the row and column numbers are assumed to be the same).
"""
function assemble!(self::SysmatAssemblerVecOfMatSymm{T}, mat::FMat{T},  dofnums::Union{FIntVec, FIntMat}, ignore::Union{FIntVec, FIntMat}) where {T<:Number}
    p = self.buffer_pointer
    @assert p <= self.buffer_length
    @assert self.elem_mat_dim == size(mat, 1)
    @assert self.elem_mat_dim == size(mat, 2)
    @assert self.elem_mat_dim == length(dofnums)
    d, m = deepcopy(vec(dofnums)), deepcopy(mat)
    fnz = 1
    @inbounds for j in 1:length(d)
        if d[j] == 0
            d[j] = fnz
            @inbounds for k in 1:length(d)
                m[j, k] = 0.0
            end
            @inbounds for k in 1:length(d)
                m[k, j] = 0.0
            end
        else
            fnz = d[j]
        end
    end
    self.buffer[p] = (d, m)
    self.buffer_pointer=p+1;
    return self
end

"""
    makematrix!(self::SysmatAssemblerVecOfMatSymm)

This function is dummy: the contents is preserved in the assembler itself.
"""
function makematrix!(self::SysmatAssemblerVecOfMatSymm)
    self.buffer_length = self.buffer_pointer-1
    return nothing
end

end
