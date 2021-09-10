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

include("FEMMShellQ4SRIModule.jl")
include("FEMMShellT3DSGOModule.jl")
include("FEMMShellT3DSGICModule.jl")
include("FEMMShellT3DSGModule.jl")
include("FEMMShellT3DSGMTModule.jl")
include("FEMMShellT3DSGAModule.jl")
include("FEMMShellCSDSG3Module.jl")
include("FEMMShellIsoPModule.jl")

end # module
