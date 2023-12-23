module FEMMPointMassModule

using LinearAlgebra: norm, Transpose, mul!
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.IntegDomainModule: IntegDomain
import FinEtoolsDeforLinear.MatDeforElastIsoModule: MatDeforElastIso

"""
    FEMMPointMass{S<:AbstractFESet}

Type for linear deformation finite element modeling machine.
"""
mutable struct FEMMPointMass{S<:AbstractFESet} <: AbstractFEMM
    integdomain::IntegDomain{S} # integration domain data
    _massmatrix::FFltMat
end

"""
    mass(self::FEMMCorotBeam,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the consistent mass matrix

This is a general routine for the abstract linear-deformation  FEMM.
"""
function mass(
    self::FEMMPointMass,
    assembler::ASS,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {ASS<:AbstractSysmatAssembler,T<:Number,TI<:Number}
    fes = self.integdomain.fes
    dofnums = zeros(FInt, 1, 6)
    startassembly!(assembler, 3 * 3 * count(fes), nalldofs(dchi), nalldofs(dchi))
    for i = 1:count(fes) # Loop over elements
        gatherdofnums!(dchi, dofnums, fes.conn[i]) # degrees of freedom
        assemble!(assembler, self._massmatrix, dofnums[1:3], dofnums[1:3])
    end # Loop over elements
    return makematrix!(assembler)
end

function mass(
    self::FEMMPointMass,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {T<:Number,TI<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return mass(self, assembler, geom0, u1, Rfield1, dchi)
end

end # module
