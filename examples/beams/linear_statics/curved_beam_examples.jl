module curved_beam_examples



using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMLinBeamModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

"""
CURVED BEAM WITH STATIC LOADS
PROBLEM DESCRIPTION
In this example a curved cantilever beam, modeled with shell elements, is
subjected to unit forces at the tip in the in-plane and out-of-plane directions, that
is, the Y and Z directions, respectively. The in-plane and out-of-plane loads are
applied in different load cases. The tip displacements in the direction of the load
are compared with independent hand calculated results.
The geometry, properties and loading are as suggested in MacNeal and Harder
1985. The cantilever beam is bent into a 90° arc. It has a 4.12 inch inner radius
and a 4.32 inch outer radius. Thus it is 0.2 inch wide and approximately 6.63 inch
long at its centerline. The beam is 0.1 inch thick in the Y direction. For modeling
in SAP2000, the curved beam is meshed into six area objects, each subtending a
15° arc.

Material Properties
E = 10,000,000 lb/in2
ν = 0.25
G = 4,000,000 lb/in2
Section Properties
Thickness = 0. 1 in
"""
function test2(direction = 2)
    E = 10000000.0 #  lb/in2
    nu = 0.25
    # Section Properties
    b = Thickness = 0.1 # in
    h = Depth = 0.2 # in
    radius = (4.32 + 4.12) / 2 # in
    k_s = 5/6 # shear correction factor
    force = 1.0 # lb
    # Select the number of elements per leg.
    nel = 4
    if direction == 2
        deflex = 0.0886
    else
    # direction = 3
        deflex = 0.5
    end

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s)

    ang=90
    members = []
    xyz=fill(0.0, 2, 3)
    for i in 1:nel
        xyz[1, :] .= (radius*cos((i-1)/nel*ang/360*2*pi), radius*sin((i-1)/nel*ang/360*2*pi), 0)
        xyz[2, :] .= (radius*cos((i)/nel*ang/360*2*pi), radius*sin((i)/nel*ang/360*2*pi), 0)
        push!(members, frame_member(xyz, 1, cs))
    end
    fens, fes = merge_members(members; tolerance = radius / 1000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    tipl = selectnode(fens; box = Float64[0 0 radius radius 0 0], inflate = radius/1000)
    basel = selectnode(fens; box = Float64[radius radius 0 0 0 0], inflate = radius/1000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, basel, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMLinBeamModule.FEMMLinBeam(IntegDomain(fes, GaussRule(1, 1)), material)
    K = FEMMLinBeamModule.stiffness(femm, geom0, u0, Rfield0, dchi);

    loadbdry = FESetP1(reshape(tipl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = fill(0.0, 6); q[direction] = force
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)
    @show dchi.values[tipl, direction], deflex
    @test norm(dchi.values[tipl, direction] .- deflex) / deflex < 0.05

    scaling = 1e0
    dchi.values .*= scaling
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[-radius -radius -radius]; [radius radius radius]]),
        plot_nodes(fens),
        plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    pl = render(plots)

    true
end


using FinEtoolsFlexStructures.FEMMCorotBeamModule

"""
CURVED BEAM WITH STATIC LOADS
PROBLEM DESCRIPTION
In this example a curved cantilever beam, modeled with shell elements, is
subjected to unit forces at the tip in the in-plane and out-of-plane directions, that
is, the Y and Z directions, respectively. The in-plane and out-of-plane loads are
applied in different load cases. The tip displacements in the direction of the load
are compared with independent hand calculated results.
The geometry, properties and loading are as suggested in MacNeal and Harder
1985. The cantilever beam is bent into a 90° arc. It has a 4.12 inch inner radius
and a 4.32 inch outer radius. Thus it is 0.2 inch wide and approximately 6.63 inch
long at its centerline. The beam is 0.1 inch thick in the Y direction. For modeling
in SAP2000, the curved beam is meshed into six area objects, each subtending a
15° arc.

Material Properties
E = 10,000,000 lb/in2
ν = 0.25
G = 4,000,000 lb/in2
Section Properties
Thickness = 0. 1 in
"""

function test3(direction = 2)
    E = 10000000.0 #  lb/in2
    nu = 0.25
    # Section Properties
    b = Thickness = 0.1 # in
    h = Depth = 0.2 # in
    radius = (4.32 + 4.12) / 2 # in
    k_s = 5/6 # shear correction factor
    force = 1.0 # lb
    # Select the number of elements per leg.
    nel = 4
        if direction == 2
        deflex = 0.0886
    else
        # direction = 3
        deflex = 0.5
    end

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s)

    ang=90
    members = []
    xyz=fill(0.0, 2, 3)
    for i in 1:nel
        xyz[1, :] .= (radius*cos((i-1)/nel*ang/360*2*pi), radius*sin((i-1)/nel*ang/360*2*pi), 0)
        xyz[2, :] .= (radius*cos((i)/nel*ang/360*2*pi), radius*sin((i)/nel*ang/360*2*pi), 0)
        push!(members, frame_member(xyz, 1, cs))
    end
    fens, fes = merge_members(members; tolerance = radius / 1000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    tipl = selectnode(fens; box = Float64[0 0 radius radius 0 0], inflate = radius/1000)
    basel = selectnode(fens; box = Float64[radius radius 0 0 0 0], inflate = radius/1000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, basel, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMCorotBeamModule.FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 1)), material)
    K = FEMMCorotBeamModule.stiffness(femm, geom0, u0, Rfield0, dchi);

    loadbdry = FESetP1(reshape(tipl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = fill(0.0, 6); q[direction] = force
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)
    @show dchi.values[tipl, direction], deflex
    @test norm(dchi.values[tipl, direction] .- deflex) / deflex < 0.05

    scaling = 1e0
    dchi.values .*= scaling
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[-radius -radius -radius]; [radius radius radius]]),
        plot_nodes(fens),
        plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    pl = render(plots)

    true
end

function allrun()
    println("#####################################################")
    println("# test2 ")
    test2(2)
    test2(3)
    println("#####################################################")
    println("# test3 ")
    test3(2)
    test3(3)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
