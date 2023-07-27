
nothing



"""
gravity loading on a straight beam
Consistent loading with distribloads_global
"""
module mgravitybeam4

using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMLinBeamModule
using FinEtoolsFlexStructures.FEMMLinBeamModule: FEMMLinBeam
stiffness = FEMMLinBeamModule.stiffness
distribloads_global = FEMMLinBeamModule.distribloads_global
inspectintegpoints = FEMMLinBeamModule.inspectintegpoints
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
    I = b * h^3 / 12
    ecc = 0.10*phun("m")
    L = 10.0*phun("m")
    k_s = 5/6 # shear correction factor
    w = 2000.0*phun("kg/m^3") * 10.0*phun("m/sec^2")
    uniform_eccentricity = [0.0, 0.0, 0.0, ecc]
    # uniform_eccentricity = [0.0, 0.0, 0.0*phun("m"), 0.0]

    # Cross-sectional properties
    # cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s) # Timoshenko
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0]) # Bernoulli

    xyz = [[L/2 0 0]; [-L/2 0 0]]
    fens, fes = frame_member(xyz, nel, cs)
    @show count(fes)

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
    femm = FEMMLinBeam(IntegDomain(fes, GaussRule(1, 1), b * h), material, uniform_eccentricity)

    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    q = fill(0.0, 3); q[2] = w * b * h
    fi = ForceIntensity(q)
    F = distribloads_global(femm, geom0, u0, Rfield0, dchi, fi);

    # Solve the static problem
    solve!(dchi, K, F)
# @show dchi.values

    K_df = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:df]
    F_d =  K_df * gathersysvec(dchi, :f)
    @show H =  F_d[dchi.dofnums[leftl[1], 1]-nfreedofs(dchi)]
    M0 = H * ecc
    @show deflw = (5 * w * b * h * L^4) / (384 * E * I)
    @show deflH = (M0 * L^2) / (8 * E * I)
    @show deflw - deflH
    @show extrema(dchi.values[:, 2])

    function inspector(idat, i, conn, ecoords, elvecf, loc)
        if conn[1] == leftl[1] || conn[2] == leftl[1]
            @info "element at support $(xyz[1, :])"
            @show elvecf
        end
    end

    inspectintegpoints(femm, geom0, dchi, NodalField([1.0]), 1:count(fes), inspector, nothing)

    scaling = 2e2
    dchi.values .*= scaling
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[-L -L -L]; [L L L]]),
        plot_nodes(fens),
        plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    pl = render(plots)

    true
end
test(8)
end # module

nothing



"""
gravity loading on a straight beam
Alternative definition of the cross section orientation
"""
module mgravitybeam5

using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMLinBeamModule
using FinEtoolsFlexStructures.FEMMLinBeamModule: FEMMLinBeam
stiffness = FEMMLinBeamModule.stiffness
distribloads_global = FEMMLinBeamModule.distribloads_global
inspectintegpoints = FEMMLinBeamModule.inspectintegpoints
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
    # Section Properties: The beam is being bent in the f1-f3 plane
    b = 0.4*phun("m")
    h = 0.2*phun("m")
    I = b^3 * h / 12
    ecc = 0.10*phun("m")
    L = 10.0*phun("m")
    k_s = 5/6 # shear correction factor
    w = 2000.0*phun("kg/m^3") * 10.0*phun("m/sec^2")
    uniform_eccentricity = [0.0, 0.0, ecc, 0.0]
    # uniform_eccentricity = [0.0, 0.0, 0.0*phun("m"), 0.0]

    # Cross-sectional properties
    # cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s) # Timoshenko
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 1.0, 0.0]) # Bernoulli

    xyz = [[L/2 0 0]; [-L/2 0 0]]
    fens, fes = frame_member(xyz, nel, cs)
    @show count(fes)

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
    femm = FEMMLinBeam(IntegDomain(fes, GaussRule(1, 1), b * h), material, uniform_eccentricity)

    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    q = fill(0.0, 3); q[2] = w * b * h
        fi = ForceIntensity(q)
        F = distribloads_global(femm, geom0, u0, Rfield0, dchi, fi);

    # Solve the static problem
    solve!(dchi, K, F)
# @show dchi.values

    K_df = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:df]
    F_d =  K_df * gathersysvec(dchi, :f)
    @show H =  F_d[dchi.dofnums[leftl[1], 1]-nfreedofs(dchi)]
    M0 = H * ecc
    @show deflw = (5 * w * b * h * L^4) / (384 * E * I)
    @show deflH = (M0 * L^2) / (8 * E * I)
    @show deflw - deflH
    @show extrema(dchi.values[:, 2])

    function inspector(idat, i, conn, ecoords, elvecf, loc)
        if conn[1] == leftl[1] || conn[2] == leftl[1]
            @info "element at support $(xyz[1, :])"
            @show elvecf
        end
    end

    inspectintegpoints(femm, geom0, dchi, NodalField([1.0]), 1:count(fes), inspector, nothing)

    scaling = 2e2
    dchi.values .*= scaling
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[-L -L -L]; [L L L]]),
        plot_nodes(fens),
        plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    pl = render(plots)

    true
end
test(8)
end # module

nothing

