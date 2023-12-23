module FEMMPointGroundedSpringModule

using LinearAlgebra: norm, Transpose, mul!
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.IntegDomainModule: IntegDomain
import FinEtoolsDeforLinear.MatDeforElastIsoModule: MatDeforElastIso

"""
    FEMMPointGroundedSpring{S<:AbstractFESet}

Type for finite element modeling machine for grounded springs.
"""
mutable struct FEMMPointGroundedSpring{S<:AbstractFESet} <: AbstractFEMM
    integdomain::IntegDomain{S} # integration domain data
    _stiffmatrix::FFltMat
end

"""
    mass(self::FEMMCorotBeam,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the stiffness matrix
"""
function stiffness(
    self::FEMMPointGroundedSpring,
    assembler::ASS,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {ASS<:AbstractSysmatAssembler,T<:Number,TI<:Number}
    fes = self.integdomain.fes
    dofnums = zeros(FInt, 1, 6)
    startassembly!(assembler, 6 * 6 * count(fes), nalldofs(dchi), nalldofs(dchi))
    for i = 1:count(fes) # Loop over elements
        gatherdofnums!(dchi, dofnums, fes.conn[i]) # degrees of freedom
        assemble!(assembler, self._stiffmatrix, dofnums, dofnums)
    end # Loop over elements
    return makematrix!(assembler)
end

function stiffness(
    self::FEMMPointGroundedSpring,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {T<:Number,TI<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi)
end

end # module
