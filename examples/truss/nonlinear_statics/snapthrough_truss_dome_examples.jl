"""
Analysis of Geometrically Nonlinear Structures Second Edition 

by Robert Levy Technion-Israel Institute of Technology, Haifa, Israel 

and William R. Spillers New Jersey Institute of Technology, 

Space truss dome in section 2.4.2

Vertical deflection at the crown: -.20641184e+00 in (linear analysis)

In nonlinear analysis we consider loading which is just one tenth of the force
considered for the linear analysis. We reach  displacement obtained in linear
analysis for just about sixty percent of the loading (loading parameter ~6). The
limit point is considerably below 1.0. Therefore the structure would not be able
withstand buckling.
"""
module snapthrough_truss_dome_examples

using FinEtools
using FinEtools.AssemblyModule: SysmatAssemblerFFBlock, SysvecAssemblerFBlock
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtoolsFlexStructures.FEMMCorotTrussModule
using FinEtoolsFlexStructures.FEMMCorotTrussModule: FEMMCorotTruss
stiffness = FEMMCorotTrussModule.stiffness
mass = FEMMCorotTrussModule.mass
geostiffness = FEMMCorotTrussModule.geostiffness
restoringforce = FEMMCorotTrussModule.restoringforce
import LinearAlgebra: dot
import LinearAlgebra: norm
using Arpack
using LinearAlgebra
using SparseArrays
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json
using PlotlyJS

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
    P = 220.46 / 10 * phun("lbf")
    maxstep = 15
    fixed_step = true
    # maxstep = 30
    # fixed_step = true
    optit = 5
    maxit = 3
    rtol = 1e-5

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
    Rfield1 = deepcopy(Rfield0) # not used by the truss mechanics, but it needs to be present
    dchi = NodalField(zeros(size(fens.xyz, 1), 3))

    # Apply EBC's
    l1 = [8, 9, 10, 11, 12, 13]
    for i in [1, 2, 3,]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)
    fr = freedofs(dchi)

    # Assemble the global discrete system
    # massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    # vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    femm = FEMMCorotTruss(IntegDomain(fes, GaussRule(1, 2)), material)

    loadbdry = FESetP1(reshape([1], 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = zeros(3)
    q[3] = -P 
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3)
    F = F[fr]

    # Auxiliary Variables 
    u1 = deepcopy(u0)
    step_length_scaling = 1.0
    dchiv = fill(0.0, length(fr))
    dchiv1 = fill(0.0, length(fr))
    dchivprev = fill(0.0, length(fr))
    dellamprev = 0.0
    lam = 0.1
    dellam1 = 0.0

    tip_displacement = Float64[0.0]
    actual_lams = Float64[0.0]

    step = 1
    while step < maxstep
        u0.values[:] .= u1.values[:]
        applyebc!(dchi) # Apply boundary conditions

        println("Step: $step")
        println("===========================================================")
        if step == 1
            K = stiffness(femm, geom0, u1, Rfield1, dchi) + geostiffness(femm, geom0, u1, Rfield1, dchi)
            dchiv1 .= K[fr, fr] \ F
            @show dellam1 = lam
        else
            dchiv1 .= step_length_scaling * dchivprev
            dellam1 = step_length_scaling * dellamprev
            lam += dellam1
        end
        dchi = scattersysvec!(dchi, dchiv1)
        u1.values[:] .= u0.values[:] .+ dchi.values[:]   # increment displacement

        dchiv .= dchiv1
        dellam = dellam1

        error = 1.0
        iter = 0
        while error > rtol
            iter += 1
            Fr = restoringforce(femm, geom0, u1, Rfield1, dchi)[fr]       # Internal forces
            K = stiffness(femm, geom0, u1, Rfield1, dchi) + geostiffness(femm, geom0, u1, Rfield1, dchi)
            d1 = K[fr, fr] \ F
            d2 = K[fr, fr] \ (lam .* F + Fr)

            @show dellamcorr = -dot(dchiv1, d2) / dot(dchiv1, d1)
            dchivcorr = dellamcorr * d1 + d2

            dellam += dellamcorr
            lam += dellamcorr

            @show lam

            dchiv += dchivcorr

            error = norm(lam * F + Fr) / norm(lam * F)
            @info "Iter $iter, balance error $error"

            if !fixed_step
                step_length_scaling = (1 / 2)^((iter - optit) / 4)
            end

            dchivprev = dchiv
            dellamprev = dellam

            dchi = scattersysvec!(dchi, dchivcorr)
            u1.values[:] .+= dchi.values[:]   # increment displacement

            if iter > maxit
                break
            end
        end

        # Store output from step
        push!(actual_lams, lam)
        push!(tip_displacement, u1.values[1, 3])

        step += 1
    end

    tip_displacement = vec(tip_displacement)
    plots = cat(
        # scatter(;x=vec(tip_displacement) ./ phun("mm"), y=analytical.(tip_displacement), mode="lines", line_color = "rgb(0, 0, 0)", name="Analytical"),
        scatter(; x=tip_displacement ./ phun("in"), y=actual_lams, mode="markers", line_color="rgb(255, 0, 0)", name="FEA"),
        dims=1)
    layout = Layout(width=700, height=500,
        xaxis=attr(title="Tip displacement [in]"),
        yaxis=attr(title="Load parameter [ND]"))
    pl = plot(plots, layout)
    display(pl)

    if visualize
        scaling = 1e1
        dchi.values .*= scaling
        radius = 20 * phun("in")
        plots = cat(plot_space_box([[-radius -radius -radius]; [radius radius radius]]),
            plot_nodes(fens),
            plot_solid(fens, fes;
                x=geom0.values, u=u1.values[:, 1:3],);
            dims=1)
        pl = render(plots)
    end

    true
end # force_1

function allrun()
    println("#####################################################")
    println("# solve ")
    solve(false)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
