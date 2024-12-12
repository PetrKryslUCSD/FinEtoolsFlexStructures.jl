"""
From Barbero's Finite Element Analysis using Abaqus ... book Example 3.1

Consider a simply supported square plate ax = ay = 2000 mm, laminated
with AS4D/9310 (Table 3.1) in a [0/90]n configuration. The total laminate thickness is
t = 10 mm for all values of n.
The plate is loaded in compression with an edge load
Nx = 1 N/mm and (Ny = Nxy = Mx = My = Mxy = 0). Tabulate the center deflection
perpendicular to the plate surface when the number of layers is n = 1, 5, 10, 15, and 20.
Calculate the A, B, D, H matrices and enter them into Abaqus. Use symmetry to model
1/4 of the plate.
"""
module barbero_3_1_examples

using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule

    ax = ay = 2000 * phun("mm")
    nx = ny = 8
    # ASFD/9310
    E1 = 133860 * phun("MPa")
    E2 = 7706 * phun("MPa")
    G12 = 4306 * phun("MPa")
    G13 = G12
    nu12 = 0.301
    nu23 = 0.396
    G23 = 2760 * phun("MPa")
    npairs = 1
    thickness = 10 * phun("mm")
    tolerance = ax / nx / 100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[]
    for p = 1:npairs
        push!(plies, CM.Ply("ply_0_$p", mater, thickness / npairs / 2, 0))
        push!(plies, CM.Ply("ply_90_$p", mater, thickness / npairs / 2, 90))
    end
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("example_3.1", plies, mcsys)

    fens, fes = T3block(ax, ay, nx, ny)
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)

    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), layup)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # Pin one of the corners
    l1 = selectnode(fens; box = Float64[ax ax 0 0 -Inf Inf], inflate = tolerance)
    for i in [1, 2, 3]
        setebc!(dchi, l1, true, i)
    end
    # Roller at the other
    l1 = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], inflate = tolerance)
    for i in [2]
        setebc!(dchi, l1, true, i)
    end
    # Simple support
    l1 = connectednodes(meshboundary(fes))
    for i in [3]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, geom0, u0, Rfield0, dchi)

    # Edge load
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, box = Float64[-Inf Inf 0 0 0 0], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(Float64[0, 0.1*phun("MPa")*thickness, 0, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 1)
    l1 = selectelem(fens, bfes, box = Float64[-Inf Inf ay ay 0 0], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(Float64[0, -0.1*phun("MPa")*thickness, 0, 0, 0, 0])
    F += distribloads(lfemm, geom0, dchi, fi, 1)

    # Solve
    solve_blocked!(dchi, K, F)
    for i = 1:3
        @show  maximum(dchi.values[:, i]) ./phun("mm")
        @show minimum(dchi.values[:, i]) ./phun("mm")
    end
    vtkwrite("plate-$npairs-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])
    @test maximum(dchi.values[:, 3]) ./ phun("mm") ≈ 0.22004349767718365
    true
end

function allrun()
    println("#####################################################")
    println("# test ")
    test()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
