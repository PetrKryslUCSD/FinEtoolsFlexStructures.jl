##  Modal analysis of the L-frame of Argyris.
# Eigenvector 1 frequency 11.2732  [Hz]
#   Eigenvector 2 frequency 30.5269  [Hz]
module argyr_l_frame_modal_examples

using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMCorotBeamModule
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
mass = FEMMCorotBeamModule.mass
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using VisualStructures: plot_points, plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, default_layout_3d
using PlotlyJS
using JSON

function argyr_l_frame_modal()
    # Parameters:
    E=71240.0;#MPa
    nu=0.31;# Poisson ratio
    rho=5e-9;
    b=3.0; h=30.0; L=240.0; # cross-sectional dimensions and length of each leg in millimeters
    # Choose the mass formulation:
    mass_type=2;

    # Reference frequencies
    reffs = [11.2732, 30.5269]
    neigvs = 2;

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 1.0, 0.0])

    # Select the number of elements per leg.
    n=8;
    members = []
    push!(members, frame_member([0 0 L; L 0 L], n, cs))
    push!(members, frame_member([L 0 L; L 0 0], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    l1 = selectnode(fens; box = [0 0 0 0 L L], tolerance = L / 10000)
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    M = mass(femm, geom0, u0, Rfield0, dchi);

    K_ff = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:ff]
    M_ff = matrix_blocked(M, nfreedofs(dchi), nfreedofs(dchi))[:ff]

    # Solve the eigenvalue problem
    d,v,nconv = eigs(K_ff, M_ff; nev=2*neigvs, which=:SM, explicittransform=:none)
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Natural frequencies: $fs [Hz]")
    println("Reference: $reffs [Hz]")

    # Visualize vibration modes
    scattersysvec!(dchi, v[:, 1])
    update_rotation_field!(Rfield0, dchi)
    tbox = plot_space_box([[-L/2 -L/2 0]; [L/2 L/2 L]])
    plots = cat(tbox, plot_nodes(fens),
        plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
        dims = 1)
    layout = default_layout_3d()
    layout[:scene][:aspectratio] = space_aspectratio(fens.xyz)
    layout[:scene][:aspectmode] = "manual"
    render(plots; layout = layout)

    return true
end # argyr_l_frame_modal

function argyr_l_frame_modal_anim()
    # Parameters:
    E=71240.0;#MPa
    nu=0.31;# Poisson ratio
    rho=5e-9;
    b=3.0; h=30.0; L=240.0; # cross-sectional dimensions and length of each leg in millimeters
    # Choose the mass formulation:
    mass_type=2;
    
    # Reference frequencies
    reffs = [11.2732, 30.5269]
    neigvs = 2;

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 1.0, 0.0])

    # Select the number of elements per leg.
    n = 4;
    members = []
    push!(members, frame_member([0 0 L; L 0 L], n, cs))
    push!(members, frame_member([L 0 L; L 0 0], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    l1 = selectnode(fens; box = [0 0 0 0 L L], tolerance = L/10000)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    M = mass(femm, geom0, u0, Rfield0, dchi);

    K_ff = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:ff]
    M_ff = matrix_blocked(M, nfreedofs(dchi), nfreedofs(dchi))[:ff]

    # Solve the eigenvalue problem
    d,v,nconv = eigs(K_ff, M_ff; nev=2*neigvs, which=:SM, explicittransform=:none)
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Natural frequencies: $fs [Hz]")
    println("Reference: $reffs [Hz]")

    # Visualize vibration modes
    mode = 2
    tbox = plot_space_box([[-L/2 -L/2 0]; [L/2 L/2 1.1*L]])
    tenv0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    plots = cat(tbox, tenv0; dims = 1)
    layout = default_layout_3d()
    layout[:scene][:camera][:eye] = Dict(:x=>1.02, :y=> 1.02, :z=> 0.836)
    layout[:scene][:camera][:center] = Dict(:x=>0.058, :y=>0.065, :z=>-0.122)
    pl = render(plots; layout = layout, title = "Mode $(mode)")
    Rfield1 = deepcopy(Rfield0)
    scale = L/3  /  maximum(v[:, mode])
    for xscale in scale.*sin.(collect(0:1:72).*(2*pi/21))
        scattersysvec!(dchi, xscale.*v[:, mode])
        Rfield1 = deepcopy(Rfield0)
        update_rotation_field!(Rfield1, dchi)
        tenv1 = plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield1.values, facecolor = "rgb(25, 255, 25)");
        plots = cat(tbox, tenv0, tenv1; dims = 1)
        react!(pl, plots, pl.plot.layout)
        sleep(0.115)
    end
    savejson(pl, "plots.json")

    return true
end # argyr_l_frame_modal_anim

function allrun()
    println("#####################################################")
    println("# argyr_l_frame_modal ")
    argyr_l_frame_modal()
    println("#####################################################")
    println("# argyr_l_frame_modal_anim ")
    argyr_l_frame_modal_anim()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
