module mcompshell0
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using Test

function test()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    # From Barbero's Finite Element Analysis using Abaqus ... book Example 3.1
    ax = ay = 2000*phun("mm")
    nx = ny = 3
    # ASFD/9310
    E1 = 133860*phun("MPa")
    E2 = 7706*phun("MPa")
    G12 = 4306*phun("MPa")
    G13 = G12
    nu12 = 0.301;
    nu23 = 0.396
    G23 = 2760*phun("MPa")
    npairs = 2
    thickness = 10*phun("mm");
    mcsys = CM.cartesian_csys((1, 2, 3))
    CM = CompositeLayupModule
    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[]
    for p in npairs
        push!(plies, CM.Ply("ply_0_$p", mater, thickness, 0))
        push!(plies, CM.Ply("ply_90_$p", mater, thickness, 90))
    end
    cl = CM.CompositeLayup("example_3.1", plies, mcsys)
    fens, fes = T3block(ax,ay,nx,ny);
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), [cl])
    true
end
end
using .mcompshell0
mcompshell0.test()
