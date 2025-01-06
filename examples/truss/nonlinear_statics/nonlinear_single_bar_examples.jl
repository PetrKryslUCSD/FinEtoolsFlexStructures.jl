"""

"""
module nonlinear_single_bar_examples

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
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json
using PlotlyJS

function solve(visualize=false)
    # Parameters:
    E = 30_000.0 * phun("MPa")
    nu = 0.3
    W = 2500.0 * phun("mm")
    H = 25.0 * phun("mm")
    x = Float64[
        0.0 0.0 0.0 # 1
        W H 0.0  # 2
    ] 
    area = 5e7 / E
    P = 0.01 * phun("kilo*N")
    maxit = 10
    maxsteps = 10
    analytical(w) = E * area * (H + w) * (2 * (H + w) * w - w^2) / (2 * W^3)

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> sqrt(area), s -> sqrt(area), s -> [0.0, 0.0, -1.0])

    # Select the number of elements per leg.
    n = 1
    members = []
    for j in (
        [1, 2],
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
    pinned = selectnode(fens; box = boundingbox(x[1, :]), tolerance = H/100)
    for i in [1, 2, 3,]
        setebc!(dchi, pinned, true, i)
    end
    loaded = selectnode(fens; box = boundingbox(x[2, :]), tolerance = H/100)
    for i in [1, 3,]
        setebc!(dchi, loaded, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)
    
    # Assemble the global discrete system
    # massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    # vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    femm = FEMMCorotTruss(IntegDomain(fes, GaussRule(1, 2)), material)
    
    K = stiffness(femm, geom0, u0, Rfield0, dchi)
    loadbdry = FESetP1(reshape(loaded, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = zeros(3)
    q[2] = P
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3)
    
    # Solve the static problem
    fr = freedofs(dchi)
    
    # Auxiliary Variables 
    u1 = deepcopy(u0)
    rhs = fill(0.0, length(fr))
    initial_dchi = fill(0.0, length(fr))
    utol = 1e-9*nfreedofs(dchi);

    F = F[fr]

    
    tip_displacement = Float64[]
    actual_load_parameters = Float64[]
    load_parameters = (0.0, -1.3);
    initial_load_parameter_incr = (load_parameters[2] - load_parameters[1]) / 5
    load_parameter = load_parameters[1]
    previous_load_parameter = load_parameter
    step = 0
    while true
        applyebc!(dchi) # Apply boundary conditions
        u1.values[:] = u0.values[:]; # guess
        
        println("Load: $load_parameter")

        converged = false;
        iter = 1;
        while true
            Fr = restoringforce(femm, geom0, u1, Rfield0, dchi);       # Internal forces
            rhs .= load_parameter .* F + Fr[fr];
            K = stiffness(femm, geom0, u1, Rfield0, dchi) + geostiffness(femm, geom0, u1, Rfield0, dchi);
            dchi = scattersysvec!(dchi, (K[fr, fr])\rhs); # Disp. incr
            u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
            if iter  ==  1 
                initial_dchi .= gathersysvec(dchi, DOF_KIND_FREE);
            else
                current_dchi = gathersysvec(dchi, DOF_KIND_FREE);
                @show load_parameter, previous_load_parameter;
                delta_lambda_initial = load_parameter - previous_load_parameter;
                delta_lambda = - dot(initial_dchi, current_dchi) / delta_lambda_initial
                load_parameter += delta_lambda
            end
            print("$iter: ||du||=$(maximum(abs.(dchi.values[:])))\n")
            if maximum(abs.(dchi.values[:])) < utol # convergence check
                converged = true;
                break;
            end
            if (iter > maxit)# bailout for failed convergence
                converged = false;
                break
            end
            iter += 1;
        end
        @show converged

        if !converged
            load_parameter = previous_load_parameter + 0.5 * (load_parameter - previous_load_parameter)
        else
            @show load_parameter, previous_load_parameter
            del_load_parameter = load_parameter - previous_load_parameter
            previous_load_parameter = load_parameter
            if step < 1
                load_parameter = previous_load_parameter + initial_load_parameter_incr
            else
                load_parameter = previous_load_parameter + del_load_parameter
            end

            step += 1
            u0.values[:] = u1.values[:]       # update the displacement
            push!(actual_load_parameters, previous_load_parameter)
            push!(tip_displacement, u1.values[loaded[1], 2])


            if step >= maxsteps
                break
            end
        end

        
    end
    
    tip_displacement = vec(tip_displacement)
    plots = cat(
    scatter(;x=vec(tip_displacement) ./ phun("mm"), y=analytical.(tip_displacement), mode="lines", line_color = "rgb(0, 0, 0)", name="Analytical"),
    scatter(;x=tip_displacement ./ phun("mm"), y=P .* actual_load_parameters, mode="markers", line_color = "rgb(255, 0, 0)", name="FEA"),
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

function solve_riks(visualize=false)
    # Parameters:
    E = 30_000.0 * phun("MPa")
    nu = 0.3
    W = 2500.0 * phun("mm")
    H = 25.0 * phun("mm")
    x = Float64[
        0.0 0.0 0.0 # 1
        W H 0.0  # 2
    ] 
    area = 5e7 / E
    P = 0.01 * phun("kilo*N")
    maxstep = 10
    maxit = 30
    rtol = 1e-9
    disp_mag = H / 5 # vertical coordinate of the crown
    step_length = 3.3
    analytical(w) = -E * area * (H + w) * (2 * (H + w) * w - w^2) / (2 * W^3)

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> sqrt(area), s -> sqrt(area), s -> [0.0, 0.0, -1.0])

    # Select the number of elements per leg.
    n = 1
    members = []
    for j in (
        [1, 2],
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
    pinned = selectnode(fens; box = boundingbox(x[1, :]), tolerance = H/100)
    for i in [1, 2, 3,]
        setebc!(dchi, pinned, true, i)
    end
    loaded = selectnode(fens; box = boundingbox(x[2, :]), tolerance = H/100)
    for i in [1, 3,]
        setebc!(dchi, loaded, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)
    fr = freedofs(dchi)
    
    # Assemble the global discrete system
    # massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    # vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    femm = FEMMCorotTruss(IntegDomain(fes, GaussRule(1, 2)), material)
    
    K = stiffness(femm, geom0, u0, Rfield0, dchi)
    loadbdry = FESetP1(reshape(loaded, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = zeros(3)
    q[2] = -P
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3)
    
    # Solve the static problem
    
    F = F[fr]

    
    # Auxiliary Variables 
    u1 = deepcopy(u0)
    dchiv = fill(0.0, length(fr))
    chivc = fill(0.0, length(fr))
    chivc1 = fill(0.0, length(fr))
    dchivprev = fill(0.0, length(fr))

    lam = 0.0
    lamprev = 0.0

    tip_displacement = Float64[0.0]
    actual_lams = Float64[0.0]

    step = 1
    while step <= maxstep
        u0.values[:] .= u1.values[:]
        applyebc!(dchi) # Apply boundary conditions

        println("Step: $step, starting from load parameter $lam")
        println("===========================================================")

        K = stiffness(femm, geom0, u1, Rfield1, dchi) + geostiffness(femm, geom0, u1, Rfield1, dchi)
        chivc1 .= K[fr, fr] \ F
        @show lam, lamprev
        @show dellam1 = step_length / sqrt(dot(chivc1, chivc1) / disp_mag^2 + 1)
        @show dircoeff = dellam1 * (dot(chivc1, dchivprev) / disp_mag^2 + (lam - lamprev))
        if dircoeff < 0.0
            dellam1 = -dellam1
        end
        lamprev = lam
        lam += dellam1

        dchiv .= dchivprev 
        u1.values[:] .= u0.values[:] .+ dchi.values[:]   # increment displacement7
        
        converged = true
        error = 1.0
        iter = 0
        while true
            iter += 1
            Fr = restoringforce(femm, geom0, u1, Rfield1, dchi)[fr]       # Internal forces
            res = lam * F + Fr
        
            error = norm(res) / norm(F)
            @info "Iter $iter, lam = $lam, bal. error $error"

            if error < rtol
                break
            end

            K = stiffness(femm, geom0, u1, Rfield1, dchi) + geostiffness(femm, geom0, u1, Rfield1, dchi)
            chivc .= K[fr, fr] \ F
            ddchiv = K[fr, fr] \ res

            mu = - (dot(ddchiv, chivc1) / disp_mag^2) / (dot(chivc, chivc1) / disp_mag^2 + 1)
            
            lam += mu

            dchiv .+= ddchiv + mu .* chivc

            dchi = scattersysvec!(dchi, dchiv)
            u1.values[:] .= u0.values[:] .+ dchi.values[:]   # increment displacement

            if iter > maxit
                converged = false
                break
            end
        end
        
        if !converged
            @warn "Nonlinear iteration did not converge"
            break
        end

        dchivprev .= dchiv
        
        # Store output from step
        push!(actual_lams, lam)
        push!(tip_displacement, u1.values[loaded[1], 2])

        step += 1
    end
    
    tip_displacement = vec(tip_displacement)
    plots = cat(
    scatter(;x=vec(tip_displacement) ./ phun("mm"), y=analytical.(tip_displacement), mode="lines", line_color = "rgb(0, 0, 0)", name="Analytical"),
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
    solve_riks(false)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
