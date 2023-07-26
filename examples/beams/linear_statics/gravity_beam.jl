

"""
gravity loading on a straight beam
"""
module mgravitybeam2

using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMLinBeamModule
using FinEtoolsFlexStructures.FEMMLinBeamModule: FEMMLinBeam
stiffness = FEMMLinBeamModule.stiffness
distribloads_global = FEMMLinBeamModule.distribloads_global
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function test(nel = 2)
    E = 12.0*phun("GPa")
    nu = 0.3
    # Section Properties
    b = Thickness = 0.2*phun("m")
    h = Depth = 0.4*phun("m")
    L = 10.0*phun("m")
    k_s = 5/6 # shear correction factor
    w = 2000.0*phun("kg/m^3") * 10.0*phun("m/sec^2")

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s) # Timoshenko
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0]) # Bernoulli

    xyz = [[-L/2 0 0]; [L/2 0 0]]
    fens, fes = frame_member(xyz, nel, cs)

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    leftl = selectnode(fens; box = initbox!([], xyz[1, :]), inflate = L/1000)
    rightl = selectnode(fens; box = initbox!([], xyz[2, :]), inflate = L/1000)
    for i in [1,2,3,4,5] # pin
        setebc!(dchi, leftl, true, i)
        setebc!(dchi, rightl, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMLinBeam(IntegDomain(fes, GaussRule(1, 1), b * h), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    q = fill(0.0, 6); q[2] = w
    fi = ForceIntensity(q)
    @show F = distribloads(femm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)
    @show extrema(dchi.values[:, 2])

    scaling = 1e1
    dchi.values .*= scaling
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[-L -L -L]; [L L L]]),
        plot_nodes(fens),
        plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    pl = render(plots)

    true
end
test(4)
end # module

nothing

