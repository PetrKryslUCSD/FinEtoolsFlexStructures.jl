module FinEtoolsFlexStructures

__precompile__(true)

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
# include("FEMMRITBeamModule.jl")

include("FEMMPointMassModule.jl")

include("FEMMPointGroundedSpringModule.jl")

include("FEMMShellT3FFModule.jl")
include("FEMMShellT3FFCompModule.jl")

end # module
