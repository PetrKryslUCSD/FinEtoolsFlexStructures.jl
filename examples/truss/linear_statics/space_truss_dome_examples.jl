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

Space truss dome in section 2.4.2
"""
module space_truss_dome_examples

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

function solve(visualize=false)
    # Parameters:
    E = 30_000_000.0 * phun("psi")
    nu = 0.3
    x = Float64[
        0.0 0.0 .32346000e1 # 1
        .49212500e1 .85239000e1 .24472000e1 # 2
        -.49212500e1 .85239000e1 .24472000e1 # 3
        -.98425000e1 0.0 .24472000e1 # 4
        -.49212500e1 -.85239000e1 .24472000e1 # 5
        .49212500e1 -.85239000e1 .24472000e1 # 6
        .98425000e1 0.0 .24472000e1 # 7
        0.0 .19685000e02 0.0 # 8
        -.17047200e02 .98425000e1 0.0 # 9
        -.17047200e02 -.98425000e1 0.0 # 10
        0.0 -.19685000e02 0.0 # 11
        .17047200e02 -.98425000e1 0.0 # 12
        .17047200e02 .98425000e1 0.0 # 13
    ] * phun("in")
    area = 0.0155 * phun("in^2")
    P = 220.46 * phun("lbf")
    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> sqrt(area), s -> sqrt(area), s -> [0.0, 0.0, -1.0])

    # Select the number of elements per leg.
    n = 1
    members = []
    for j in (
        [1, 2],
        [1, 3],
        [1, 4],
        [1, 5],
        [1, 6],
        [1, 7],
        [6, 12],
        [7, 12],
        [7, 13],
        [2, 13],
        [2, 8],
        [3, 8],
        [3, 9],
        [4, 9],
        [4, 10],
        [5, 10],
        [5, 11],
        [6, 11],
        [2, 3],
        [3, 4],
        [4, 5],
        [5, 6],
        [6, 7],
        [7, 2]
    )
        push!(members, frame_member(x[j, :], n, cs))
    end
    fens, fes = merge_members(members; tolerance=1 / 1000)

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 3))

    # Apply EBC's
    l1 = [8, 9, 10, 11, 12, 13]
    for i in [1, 2, 3,]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)
    
    # Assemble the global discrete system
    # massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    # vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    femm = FEMMCorotTruss(IntegDomain(fes, GaussRule(1, 2)), material)
    
    K = stiffness(femm, geom0, u0, Rfield0, dchi)
    loadbdry = FESetP1(reshape([1], 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = zeros(3)
    q[3] = -P
    fi = ForceIntensity(q)
    @show F = distribloads(lfemm, geom0, dchi, fi, 3)

    # Solve the static problem
    fr = freedofs(dchi)
    U_f = K[fr, fr] \ F[fr]
    scattersysvec!(dchi, U_f)
    @show dchi.values[1, :] ./ phun("in")
    

    if visualize
        scaling = 1e1
        dchi.values .*= scaling
        radius = 20 * phun("in")
        plots = cat(plot_space_box([[-radius -radius -radius]; [radius radius radius]]),
                    plot_nodes(fens),
                    plot_solid(fens, fes;
                               x = geom0.values, u = dchi.values[:, 1:3], );
                    dims = 1)
        pl = render(plots)
    end
    
    true
end # force_1

function allrun()
    println("#####################################################")
    println("# solve ")
    solve(true)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
