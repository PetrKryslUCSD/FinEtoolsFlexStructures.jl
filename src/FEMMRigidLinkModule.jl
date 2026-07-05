"""
Module for manipulation of the rigid link finite element modeling machine (FEMM). 
"""
module FEMMRigidLinkModule

using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.MatrixUtilityModule: add_n1n2t!
using FinEtools.IntegDomainModule: IntegDomain

"""
    mutable struct FEMMRigidLink{ID<:IntegDomain{S} where {S<:FESetP1}} <: AbstractFEMM

Type for linear rigid-link finite element modeling machine.

    
The stiffness matrix is computed as

```math
    K =\begin{bmatrix} 
        C^T\Gamma C & -C^T\Gamma  \\
        -\Gamma C & \Gamma \\
    \end{bmatrix}.
```

Here ``C`` is a matrix computed from the vector  ``r = h e_x``, 
which is the difference between the location of 
the subordinate and the location of the master.

In three dimensions

```math
    C =\begin{bmatrix} 
        1 & \widetilde{r} \\
        0 & 1 \\
    \end{bmatrix}.
```

Here ``\widetilde{r}`` is a skew matrix corresponding to the 
vector ``r``, and ``0`` and ``1`` stand for ``3\times3`` zero 
and identity matrices respectively.

Further, ``\Gamma`` is a diagonal matrix, such that 
the diagonal entries provide penalty on the difference 
between the individual degrees of freedom.

Reference: APPLICATION OF RIGID LINKS  IN STRUCTURAL DESIGN MODELS,
Sergey Yu. Fialko, International Journal for Computational Civil 
and Structural Engineering, 13(3) 119-137 (2017). 

"""
mutable struct FEMMRigidLink{ID<:IntegDomain{S} where {S<:FESetP1}} <: AbstractFEMM
    integdomain::ID # integration domain data
    master::Int
    penalty::Float64
end

"""
    stiffness(self::FEMMRigidLink, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute and assemble the rigid link stiffness matrix. For each pair master-subordinate,
the matrix is computed and assembled.
"""
function stiffness(
    self::FEMMRigidLink,
    assembler::ASS,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {ASS<:AbstractSysmatAssembler,T<:Number,TI<:Number}
    fes = self.integdomain.fes
    sdim = ndofs(geom0)
    xM = reshape(copy(geom0.values[self.master, :]), 1, sdim)
    xS = copy(xM)
    Id = I(sdim)
    if sdim == 2
        C = zeros(T, 3, 3)
        C .= Id
    else
        C = zeros(T, 6, 6)
        C[1:3, 1:3] .= Id
        C[4:6, 4:6] .= Id
    end
    Gamma = self.penalty * I(ndofs(dchi))
    elmat = zeros(T, 2 * ndofs(dchi), 2 * ndofs(dchi))
    Sdofnums = zeros(FInt, 1, ndofs(dchi))
    Mdofnums = zeros(FInt, 1, ndofs(dchi))
    gatherdofnums!(dchi, Mdofnums, (self.master,)) # degrees of freedom
    startassembly!(assembler, size(elmat)..., count(fes), nalldofs(dchi), nalldofs(dchi))
    for i in eachindex(fes) # Loop over "elements" (points representing the subordinates)
        if fes.conn[i][1] != self.master
            gathervalues_asmat!(geom0, xS, fes.conn[i])
            r = xS - xM
            fill!(elmat, 0.0) # Initialize element matrix
            if sdim == 2
                rx, ry = r[1], r[2]
                C[1, 3] = -ry
                C[2, 3] = +rx
            else
                rx, ry, rz = r[1], r[2], r[3]
                C[1, 5] = +rz
                C[1, 6] = -ry
                C[2, 4] = -rz
                C[2, 6] = +rx
                C[3, 4] = +ry
                C[3, 5] = -rx
            end
            elmat .= vcat(
                hcat(transpose(C) * (Gamma * C), -transpose(C) * Gamma),
                hcat(-Gamma * C, Gamma)
            )
            gatherdofnums!(dchi, Sdofnums, fes.conn[i]) # degrees of freedom
            dofnums = hcat(Mdofnums, Sdofnums)
            assemble!(assembler, elmat, dofnums, dofnums)
        end
    end # Loop over elements
    return makematrix!(assembler)
end

function stiffness(
    self::FEMMRigidLink,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {T<:Number,TI<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi)
end

end # module
