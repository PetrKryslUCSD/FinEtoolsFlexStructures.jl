module FEMMPointMassModule

using LinearAlgebra: norm, Transpose, mul!
using FinEtools
using FinEtools.IntegDomainModule: IntegDomain
import FinEtoolsDeforLinear.MatDeforElastIsoModule: MatDeforElastIso
using ..FESetCorotBeamModule: FESetL2CorotBeam, local_frame_and_def!, local_mass!, local_stiffness!, natural_forces!, local_geometric_stiffness!, local_forces!, MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA


"""
    FEMMPointMass{S<:AbstractFESet}

Class for linear deformation finite element modeling machine.
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
function mass(self::FEMMPointMass, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    dofnums = zeros(FInt, 1, 6); 
    startassembly!(assembler, 3, 3, count(fes), dchi.nfreedofs, dchi.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, self._massmatrix, dofnums[1:3], dofnums[1:3]); 
    end # Loop over elements
    return makematrix!(assembler);
end

function mass(self::FEMMPointMass, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return mass(self, assembler, geom0, u1, Rfield1, dchi);
end

end # module
