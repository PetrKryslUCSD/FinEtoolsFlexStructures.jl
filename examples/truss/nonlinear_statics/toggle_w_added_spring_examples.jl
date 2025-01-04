"""
Two-bar toggle.
Path-tracing simulation.
"""
module toggle_w_added_spring_examples

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
    E = 5e6
    nu = 0.3
    W = 10.0 * phun("m")
    H = 0.5 * phun("m")
    x = Float64[
        0 H 0
        -W 0 0
        W 0 0
    ] 
    area = 1.0
    P = 100
    
    maxstep = 13
    fixed_step = false
    # maxstep = 30
    # fixed_step = true
    optit = 2
    maxit = 3
    rtol = 1e-5
    analytical(w) = E * area * (H + w) * (2 * (H + w) * w - w^2) / (2 * W^3)

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> sqrt(area), s -> sqrt(area), s -> [0.0, 0.0, -1.0])

    # Select the number of elements per leg.
    n = 1
    members = []
    for j in (
        [1, 2],
        [1, 3],
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
    pinned = vcat(selectnode(fens; box = boundingbox(x[2, :]), tolerance = H/100),
                  selectnode(fens; box = boundingbox(x[3, :]), tolerance = H/100))
    for i in [1, 2, 3,]
        setebc!(dchi, pinned, true, i)
    end
    loaded = selectnode(fens; box = boundingbox(x[1, :]), tolerance = H/100)
    for i in [1, 3,]
        setebc!(dchi, loaded, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)
    fr = freedofs(dchi)
    
    femm = FEMMCorotTruss(IntegDomain(fes, GaussRule(1, 2)), material)
    
    K = stiffness(femm, geom0, u0, Rfield0, dchi)
    loadbdry = FESetP1(reshape(loaded, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = zeros(3)
    q[2] = -P
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
    lam = 1.0
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
            dellam1 = lam
        else
            dchiv1 .= step_length_scaling * dchivprev
            dellam1 = step_length_scaling * dellamprev
            lam += dellam1
        end
        dchi = scattersysvec!(dchi, dchiv1); 
        u1.values[:]  .= u0.values[:] .+ dchi.values[:];   # increment displacement
        
        dchiv .=  dchiv1
        dellam = dellam1
        
        error = 1.0
        iter = 0
        while error > rtol
            iter += 1
            Fr = restoringforce(femm, geom0, u1, Rfield1, dchi)[fr]       # Internal forces
            K = stiffness(femm, geom0, u1, Rfield1, dchi) + geostiffness(femm, geom0, u1, Rfield1, dchi)
            d1 = K[fr, fr] \ F
            d2 = K[fr, fr] \ (lam .* F + Fr)
                        
            dellamcorr = -dot(dchiv1, d2) / dot(dchiv1, d1)
            dchivcorr = dellamcorr*d1 + d2
            
            dellam += dellamcorr
            lam += dellamcorr
      
            dchiv += dchivcorr

            error  = norm(lam*F + Fr) / norm(lam*F)
            @info "Iter $iter, balance error $error"

            if !fixed_step
                step_length_scaling = (1 / 2)^((iter - optit) / 4)
            end
            
            dchivprev = dchiv
            dellamprev  = dellam

            dchi = scattersysvec!(dchi, dchivcorr); 
            u1.values[:]  .+= dchi.values[:];   # increment displacement
            
            if iter > maxit
                break
            end
        end

        # Store output from step
        push!(actual_lams, lam)
        push!(tip_displacement, u1.values[loaded[1], 2])

        step += 1
    end
    
    tip_displacement = vec(tip_displacement)
    plots = cat(
    # scatter(;x=vec(tip_displacement) ./ phun("mm"), y=analytical.(tip_displacement), mode="lines", line_color = "rgb(0, 0, 0)", name="Analytical"),
    scatter(;x=tip_displacement ./ phun("mm"), y=P .* actual_lams, mode="markers", line_color = "rgb(255, 0, 0)", name="FEA"),
    dims = 1)
    layout = Layout(width=700, height=500, 
        xaxis=attr(title="Tip displacement [mm]"),
        yaxis=attr(title="Load [N]"))
    pl = plot(plots, layout)
    display(pl)

    if visualize
        scaling = 1e1
        dchi.values .*= scaling
        radius = 20 * phun("in")
        plots = cat(plot_space_box([[-radius -radius -radius]; [radius radius radius]]),
                    plot_nodes(fens),
                    plot_solid(fens, fes;
                               x = geom0.values, u = u1.values[:, 1:3], );
                    dims = 1)
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
