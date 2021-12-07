module FinEtoolsFlexStructures

__precompile__(true)

using FinEtools

include("RotUtilModule.jl")
include("VisUtilModule.jl")
include("TransformerModule.jl")

include("CrossSectionModule.jl")

include("FESetCorotBeamModule.jl")
include("FESetShellQ4Module.jl")
include("FESetShellT3Module.jl")

include("MeshFrameMemberModule.jl")

include("FEMMCorotBeamModule.jl")

include("FEMMPointMassModule.jl")

include("FEMMPointGroundedSpringModule.jl")

include("FEMMShellT3FFModule.jl")
include("FEMMShellT3FFAModule.jl")
include("FEMMShellT3FFCModule.jl")

end # module
