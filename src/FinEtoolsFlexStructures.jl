module FinEtoolsFlexStructures

__precompile__(true)

using FinEtools

include("RotUtilModule.jl")
include("VisUtilModule.jl")

include("CrossSectionModule.jl")

include("FESetCorotBeamModule.jl")
include("FESetShellQ4Module.jl")
include("FESetShellT3Module.jl")

include("MeshFrameMemberModule.jl")

include("FEMMCorotBeamModule.jl")
include("FEMMPointMassModule.jl")
include("FEMMPointGroundedSpringModule.jl")
include("FEMMShellQ4SRIModule.jl")
include("FEMMShellDSG3Module.jl")
include("FEMMShellDSG3IModule.jl")
include("FEMMShellCSDSG3Module.jl")
include("FEMMShellT3Module.jl")

end # module
