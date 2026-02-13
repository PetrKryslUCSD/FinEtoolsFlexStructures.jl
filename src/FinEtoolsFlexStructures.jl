"""
FinEtoolsFlexStructures (C) 2020-2024, Petr Krysl

FinEtools used for

- Simulations of large-displacement response of three-dimensional flexible-beam
structures. Linear static analysis, modal analysis, linear buckling analysis.
Nonlinear statics and dynamics; 
- Simulations of shell structures. Linear static analysis, modal analysis,
explicit dynamic analysis. Shells can be homogeneous or layered (laminated,
composite).
"""
module FinEtoolsFlexStructures

# __precompile__(true)

using FinEtools

include("RotUtilModule.jl")
include("TransformerModule.jl")
include("AssemblyModule.jl")

include("CrossSectionModule.jl")
include("CompositeLayupModule.jl")

include("FESetL2BeamModule.jl")
include("FESetShellQ4Module.jl")
include("FESetShellT3Module.jl")

include("MeshFrameMemberModule.jl")

include("FEMMCorotBeamModule.jl")
include("FEMMLinBeamModule.jl")
include("FEMMRITBeamModule.jl")

include("FEMMPointMassModule.jl")

include("FEMMPointGroundedSpringModule.jl")

include("FEMMShellT3FFModule.jl")
include("FEMMShellT3FFCompModule.jl")

include("FEMMCorotTrussModule.jl")

include("FEMMShellQ4RNTModule.jl")

# Enable LSP look-up in test modules.
if false include("../test/runtests.jl") end

end # module
