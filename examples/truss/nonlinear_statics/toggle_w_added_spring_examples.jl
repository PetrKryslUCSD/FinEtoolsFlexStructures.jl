"""

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

Base.:/(s::Tuple{Vector{Float64},Float64}, m::Float64) = (s[1] / m, s[2] / m)
Base.:*(s::Tuple{Vector{Float64},Float64}, m::Float64) = (s[1] * m, s[2] * m)
Base.:-(sa::Tuple{Vector{Float64},Float64}, sb::Tuple{Vector{Float64},Float64}) = (sa[1] - sb[1], sa[2] - sb[2])
dot(sola::Tuple{Vector{Float64}}, solb::Tuple{Vector{Float64}}) = (dot(sola[1], solb[1]) + sola[2] * solb[2])
norm(sola::Tuple{Vector{Float64}}, solb::Tuple{Vector{Float64}}) = sqrt(dot(sola, solb))

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
    
    maxstep = 5
    optiter = 5
    maxit = 3
    fixed_step = true
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
    R = fill(0.0, length(fr))
    rtol = 1e-5

    
    maxFactor = 1.0e20
    totalFactor = 1.0
    factor = 1.0
    Da = fill(0.0, length(fr))
    Da1 = fill(0.0, length(fr))
    Daprev = fill(0.0, length(fr))
    Dlamprev = 0.0
    lam = 1.0
    Dlam1 = 0.0

    tip_displacement = Float64[0.0]
    actual_lams = Float64[0.0]
    
    step = 1
    while true
        u0.values[:] .= u1.values[:]
        applyebc!(dchi) # Apply boundary conditions
        
        # Predictor
        println("Step: $step")
        println("========================================")
        if step == 1    
            K = stiffness(femm, geom0, u1, Rfield0, dchi) + geostiffness(femm, geom0, u1, Rfield0, dchi)
            Da1 .= K[fr, fr] \ F
            @show Dlam1  = lam
        else
            Da1   .= factor * Daprev
           @show  Dlam1  = factor * Dlamprev
           @show  lam += Dlam1
        end
        dchi = scattersysvec!(dchi, Da1); 
        u1.values[:]  .= u0.values[:] .+ dchi.values[:];   # increment displacement
        
        Da .=  Da1
        Dlam = Dlam1
        
        error = 1.0
        iter = 0
        while error > rtol
            iter += 1
            Fr = restoringforce(femm, geom0, u1, Rfield0, dchi)[fr]       # Internal forces
            K = stiffness(femm, geom0, u1, Rfield0, dchi) + geostiffness(femm, geom0, u1, Rfield0, dchi)
            @show d1 = K[fr, fr] \ F
            @show lam
            @show d2 = K[fr, fr] \ (lam .* F + Fr)
                        
            @show ddlam = -dot(Da1, d2) / dot(Da1, d1)
            dda = ddlam*d1 + d2
            
            Dlam += ddlam
            @show lam += ddlam
      
            Da += dda
            # a  += dda

            R = lam*F + Fr 

            error  = norm(R) / norm(lam*F)
            println("    Iteration $iter ........... : $error")

            if !fixed_step
                factor = (1 / 2)^(0.25 * (iter - optiter))
                totalFactor *= factor
            end
            factor = if totalFactor > maxFactor
                1.0
            else
                factor
            end

            Daprev = Da
            Dlamprev  = Dlam

            dchi = scattersysvec!(dchi, dda); 
            u1.values[:]  .+= dchi.values[:];   # increment displacement
            
            if iter > maxit
                break
            end
        end
        push!(actual_lams, lam)
        push!(tip_displacement, u1.values[loaded[1], 2])
        step += 1
        if step > maxstep
            break
        end
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
