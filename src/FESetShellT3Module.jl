"""
Module for a 3-node shell element set.
"""
module FESetShellT3Module

using FinEtools

"""
    struct FESetShellT3 <: AbstractFESet2Manifold{3}

Type of a finite element set for the three-node shell.
"""
struct FESetShellT3 <: AbstractFESet2Manifold{3} end

end # module
