"""
    AssemblyModule
Module for assemblers  of system matrices and vectors.
"""
module AssemblyModule

__precompile__(true)

using FinEtools
import FinEtools.AssemblyModule: startassembly!, assemble!, makematrix!
using SparseArrays
using SparseMatricesCSR
using  LinearAlgebra: diag

"""
    SysmatAssemblerSparseCSRSymm{A<:AbstractSysmatAssembler} <: AbstractSysmatAssembler

Assembler for a **symmetric square** matrix  assembled from symmetric square
matrices. The matrix is stored in the Compressed Sparse Row (CSR) format.
"""
mutable struct SysmatAssemblerSparseCSRSymm{A<:AbstractSysmatAssembler} <: AbstractSysmatAssembler
  a::A
end

function SysmatAssemblerSparseCSRSymm(zero::T=0.0) where {T<:Number}
    return SysmatAssemblerSparseCSRSymm(SysmatAssemblerSparse(zero))
end

function startassembly!(self::SysmatAssemblerSparseCSRSymm,
    expected_ntriples::IT,
    row_nalldofs::IT,
    col_nalldofs::IT;
    force_init = false) where {IT <: Integer}
    self.a = startassembly!(self.a, expected_ntriples, row_nalldofs, col_nalldofs; force_init)
    return self
end

function assemble!(self::SysmatAssemblerSparseCSRSymm,
    mat::MBT,
    dofnums_row::CIT,
    dofnums_col::CIT) where {MBT, CIT}
    assemble!(self.a, mat, dofnums_row, dofnums_col)
    return self
end

function makematrix!(self::SysmatAssemblerSparseCSRSymm)
    S = makematrix!(self.a)
    I, J, V = findnz(S)
    # Use the CSR package.
    S = SparseMatricesCSR.sparsecsr(I, J, V, size(S)...)
    return S
end


end

