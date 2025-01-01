"""
Analysis of Geometrically 
Nonlinear Structures 
Second Edition 

by 
Robert Levy 
Technion-Israel Institute of Technology, 
Haifa, Israel 

and 
William R. Spillers 
New Jersey Institute of Technology, 

Two bar buckling example on page 12
"""
module two_bars_examples

using FinEtools
using FinEtools.AssemblyModule: SysmatAssemblerFFBlock, SysvecAssemblerFBlock
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMCorotTrussModule
using FinEtoolsFlexStructures.FEMMCorotTrussModule: FEMMCorotTruss
stiffness = FEMMCorotTrussModule.stiffness
mass = FEMMCorotTrussModule.mass
geostiffness = FEMMCorotTrussModule.geostiffness
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function force_1(direction = 1, visualize = false)
    # Parameters:
    E = 1000000.0;
    nu = 0.3;
    L =   30.0; # Length of the bar
    b =  0.5; # width
    h =  4.0; # height
    magn = -0.2056167583560*1e4/L^2;
    neigvs = 2;

    reffs = [48.5475, 124.1839]

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, -1.0, 0.0])

    # Select the number of elements per leg.
    n = 1;
    members = []
    push!(members, frame_member([0 0 0; 0 L 0; ], n, cs))
    push!(members, frame_member([0 L 0; L L 0; ], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000);
    @show count(fens), count(fes)
    
    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    @show l1 = selectnode(fens; box = Float64[0 0 0 0 0 0], tolerance = L/10000)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[L L L L 0 0], tolerance = L/10000)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    for i in eachindex(fens) # no out of plane displacements, no rotations
        setebc!(dchi, [i], true, 3)
        setebc!(dchi, [i], true, 4)
        setebc!(dchi, [i], true, 5)
        setebc!(dchi, [i], true, 6)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    @show dchi

    # Assemble the global discrete system
    # massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    # vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # @show fens, fes
    femm = FEMMCorotTruss(IntegDomain(fes, GaussRule(1, 2)), material)
    # @show femm4#
    @show K = stiffness(femm, geom0, u0, Rfield0, dchi);
    elbown = selectnode(fens; box = Float64[L L L L 0 0], tolerance = L/10000)
    loadbdry = FESetP1(reshape(elbown, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = zeros(6, 1)
    q[direction] = -magn*b*h/magn_scale
    fi = ForceIntensity(Float64[0, -magn*b*h/magn_scale, 0, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    @show fr = freedofs(dchi)
    @show K[fr, fr]
    U_f = K[fr, fr] \ F[fr]
    scattersysvec!(dchi, U_f)
    @show  dchi.values[elbown, :]

    # Compute the geometric stiffness
    u0.values = dchi.values[:,1:3]
    update_rotation_field!(Rfield0, dchi)
    Kg = geostiffness(femm, massem, geom0, u0, Rfield0, dchi);

    # Solve the eigenvalue problem
    d,v,nconv = eigs(-Kg[fr, fr], K[fr, fr]; nev=neigvs, which=:LM, explicittransform=:none)
    fs = 1.0 ./ (d) ./ magn_scale
    println("Buckling factors: $fs [ND]")
    println("Reference: $reffs [ND]")

    # Visualize vibration modes
    if visualize
        scattersysvec!(dchi, 10*v[:, 1])
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[-L/2 -L/2 0]; [L/2 L/2 1.1*L]]),
                    plot_nodes(fens),
                    plot_solid(fens, fes;
                               x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
                    dims = 1)
        pl = render(plots)
        save_to_json(pl, "plot.json")
        # @show p.plot.data
    end
    
    true
end # force_1

function force_1_trimmed(direction = 2, visualize = false)
    # Parameters:
    E = 1000000.0;
    nu = 0.3;
    L =   30.0; # Length of the bar
    b =  0.5; # width
    h =  0.5; # height
    P = 10
    neigvs = 2;

    reffs = [Inf, E*b*h]

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, -1.0, 0.0])

    # Select the number of elements per leg.
    n = 1;
    members = []
    push!(members, frame_member([-L 0 0; -L L 0; ], n, cs))
    push!(members, frame_member([-L L 0; 0 L 0; ], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 3))

    # Apply EBC's
    @show l1 = selectnode(fens; box = Float64[-L -L 0 0 0 0], tolerance = L/10000)
    for i in [1,2,3,]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 0 L L 0 0], tolerance = L/10000)
    for i in [1,2,3,]
        setebc!(dchi, l1, true, i)
    end
    for i in eachindex(fens) # no out of plane displacements, no rotations
        setebc!(dchi, [i], true, 3)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    @show dchi

    # Assemble the global discrete system
    # massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    # vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # @show fens, fes
    femm = FEMMCorotTruss(IntegDomain(fes, GaussRule(1, 2)), material)
    # @show femm4#
    @show K = stiffness(femm, geom0, u0, Rfield0, dchi);
    @show elbown = selectnode(fens; box = Float64[-L -L L L 0 0], tolerance = L/10000)
    loadbdry = FESetP1(reshape(elbown, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = zeros(3)
    q[direction] = P
    fi = ForceIntensity(q);
    @show F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    @show fr = freedofs(dchi)
    @show K[fr, fr]
    @show U_f = K[fr, fr] \ F[fr]
    scattersysvec!(dchi, U_f)
    @show  dchi.values[elbown, :]

    # Compute the geometric stiffness
    u0.values = dchi.values[:,1:3]
    # update_rotation_field!(Rfield0, dchi)
    Kg = geostiffness(femm, geom0, u0, Rfield0, dchi);

    # Solve the eigenvalue problem
    d,v = eigen(Matrix(Kg[fr, fr]), Matrix(K[fr, fr]))
    @show Matrix(K[fr, fr])
    @show Matrix(Kg[fr, fr])
    fs = P ./ (d) 
    println("Buckling factors: $fs [ND]")
    println("Reference: $(reffs) [ND]")

    # Visualize vibration modes
    if visualize
        scattersysvec!(dchi, 10*v[:, 1])
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[-L/2 -L/2 0]; [L/2 L/2 1.1*L]]),
                    plot_nodes(fens),
                    plot_solid(fens, fes;
                               x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
                    dims = 1)
        pl = render(plots)
        save_to_json(pl, "plot.json")
        # @show p.plot.data
    end
    
    true
end # force_1

function allrun()
    println("#####################################################")
    println("# force_1_trimmed ")
    force_1_trimmed()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
