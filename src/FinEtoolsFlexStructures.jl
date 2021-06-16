module FinEtoolsFlexStructures

__precompile__(true)

using FinEtools

include("RotUtilModule.jl")
include("VisUtilModule.jl")

include("CrossSectionModule.jl")
include("FESetCorotBeamModule.jl")
include("MeshFrameMemberModule.jl")
include("FEMMCorotBeamModule.jl")

include("FEMMPointMassModule.jl")

include("FEMMPointGroundedSpringModule.jl")

include("FESetShellQ4SRIModule.jl")

end # module
