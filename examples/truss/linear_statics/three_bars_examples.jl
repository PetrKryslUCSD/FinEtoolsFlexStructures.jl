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

Three bar example on page 32
"""
module three_bars_examples

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

function solve(direction=2, visualize=false)
    # Parameters:
    E = 30000000.0
    nu = 0.3
    x = Float64[10 20 0
        0 20 0
        0 10 0
        10 0 0]
        
    area = .65700000E+00
    ux = -0.333259380e-02
    uy = -0.159162100e-02
    FORCE = [
        -0.656854250+04
        -0.485281370+04
        -0.156854250+04
    ]





    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> area, s -> 1, s -> [0.0, 0.0, -1.0])

    # Select the number of elements per leg.
    n = 1
    members = []
    push!(members, frame_member(x[[1, 2], :], n, cs))
    push!(members, frame_member(x[[1, 3], :], n, cs))
    push!(members, frame_member(x[[1, 4], :], n, cs))
    fens, fes = merge_members(members; tolerance=1 / 10000)

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 3))

    # Apply EBC's
    l1 = [2, 3, 4]
    for i in [1, 2, 3,]
        setebc!(dchi, l1, true, i)
    end
    l1 = [1]
    for i in [ 3,]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)
    @show dchi

    # Assemble the global discrete system
    # massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    # vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # @show fens, fes
    femm = FEMMCorotTruss(IntegDomain(fes, GaussRule(1, 2)), material)
    # @show femm4#
    @show K = stiffness(femm, geom0, u0, Rfield0, dchi)
    loadbdry = FESetP1(reshape([1], 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = zeros(3)
    q[1] = -1.0e4
    q[2] = -1.0e4/2
    fi = ForceIntensity(q)
    @show F = distribloads(lfemm, geom0, dchi, fi, 3)

    # Solve the static problem
    fr = freedofs(dchi)
    K[fr, fr]
    @show U_f = K[fr, fr] \ F[fr]
    scattersysvec!(dchi, U_f)
    @show dchi.values[1, :]
    @show ux, uy

    # Compute the geometric stiffness
    u0.values = dchi.values[:, 1:3]
    

    true
end # force_1

function allrun()
    println("#####################################################")
    println("# solve ")
    solve()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
