module FESetShellT3Module

using FinEtools
using FinEtools.MeshQuadrilateralModule
import FinEtools.FESetModule: cat, subset, nodesperelem
using FinEtools.MatrixUtilityModule: complete_lt! 
using LinearAlgebra: norm, Transpose, mul!

"""
    struct FESetShellT3 <: AbstractFESet2Manifold{3}

Type of a finite element set for the three-node shell.
"""
struct FESetShellT3 <: AbstractFESet2Manifold{3} 
end

end # module
