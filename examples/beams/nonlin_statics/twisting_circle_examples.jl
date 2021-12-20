module twisting_circle_examples
# Twisting of a circle into shape with two turns
# G. Rebel; Finite Rotation Shell Theory including Drill Rotations and its 
# Finite Element Implementation; PhD Thesis, Delft University Press (1998).

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMCorotBeamModule
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
geostiffness = FEMMCorotBeamModule.geostiffness
mass = FEMMCorotBeamModule.mass
distribloads_global = FEMMCorotBeamModule.distribloads_global
restoringforce = FEMMCorotBeamModule.restoringforce
gyroscopic = FEMMCorotBeamModule.gyroscopic
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Examples.VisUtilModule: plot_points, plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio
using PlotlyJS
using JSON


num3dig(a) = round(a * 1000) / 1000

function twisting_circle()
    # Parameters:
    E = 200000.0;
    nu = 0.3;# Poisson ratio
    rho = 7.8e-9;
    b = 0.6; h = 6.0;  # cross-sectional dimensions 
    radius = 120 # radius of the circle
    utol = 1e-3;
    maxit = 30
    rotmag = 2*pi
    # Select the number of elements around the circumference.
    nel = 64;
    
    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> h, s -> b, s -> [0.0, 0.0, 1.0])

    tolerance = radius/nel/1000;
    fens, fes = frame_member([0 0 0; 2*pi 0 0], nel, cs)
    for i in 1:count(fens)
        a = fens.xyz[i, 1]
        fens.xyz[i, :] .= (radius+radius*cos(a), radius*sin(a), 0)
    end
    fens, fes = mergenodes(fens, fes, tolerance, [1, nel+1])
    
    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, 0.0, E, nu, 0.0)
    
    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))
    
    # Locate some nodes
    box = fill(0.0, 6)
    # The pinned end is
    initbox!(box, [0.0,0.0,0.0])
    clampedn = selectnode(fens; box = box, inflate = tolerance)
    # The force is applied at the elbow
    initbox!(box, [2*radius,0,0])
    torquen = selectnode(fens; box = box, inflate = tolerance)

    # Apply EBC's
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, clampedn, true, i)
    end
    for i in [2, 3, 4, 5, 6]
        setebc!(dchi, torquen, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    
    # Additional fields
    u1 = deepcopy(u0)
    Rfield1 = deepcopy(Rfield0)
    rhs = gathersysvec(dchi);
    TMPv = deepcopy(rhs)
    utol = 1e-13*dchi.nfreedofs;
    
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    
    tbox = plot_space_box([[0 -radius -radius]; [2.5*radius radius radius]])
    tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    plots = cat(tbox,  tshape0; dims = 1)
    pl = render(plots; width = 800, height = 800)
    sleep(0.5)
    # savefig(pl, "twisting_circle-$(0).png", width = 300, height = 300)

    load_parameter_delta = 0.002
    load_parameters = 0.0:load_parameter_delta:1;
    step = 0
    for load_parameter in load_parameters
        dchi.values[:] .= 0.0
        dchi.fixed_values[torquen[1], 4] = rotmag * load_parameter_delta
        applyebc!(dchi) # Apply boundary conditions
        u1.values[:] = u0.values[:]; # guess
        Rfield1.values[:] = Rfield0.values[:]; # guess
        update_rotation_field!(Rfield1, dchi)

        println("Load: $load_parameter")
        iter = 1;
        while true
            Fr = restoringforce(femm, geom0, u1, Rfield1, dchi);       # Internal forces
            rhs = Fr;
            K = stiffness(femm, geom0, u1, Rfield1, dchi) + geostiffness(femm, geom0, u1, Rfield1, dchi);
            dchi = scattersysvec!(dchi, (K)\rhs); # Disp. incr
            dchi.values[torquen[1], 4] = 0.0
            u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
            update_rotation_field!(Rfield1, dchi)
            print("$iter: ||du||=$(maximum(abs.(dchi.values[:])))\n")
            if maximum(abs.(dchi.values[:])) < utol # convergence check
                break;
            end
            if (iter > maxit)# bailout for failed convergence
                error("Possible failed convergence");
            end
            iter += 1;
        end
        u0.values[:] = u1.values[:];       # update the displacement
        Rfield0.values[:] = Rfield1.values[:]; # update the rotations
        
        
        step = step + 1
        if mod(step, 5) == 1
            @show step
            tshape1 = plot_solid(fens, fes; x = geom0.values, u = u1.values, R = Rfield1.values, facecolor = "rgb(125, 15, 15)");
            plots = cat(tbox,  tshape0,  tshape1; dims = 1)
            pl.plot.layout[:title] = "Angle=$(num3dig(rotmag * load_parameter / pi))*pi"
            react!(pl, plots, pl.plot.layout)
            sleep(0.25)
            # savefig(pl, "twisting_circle-$(step).png", width = 300, height = 300)
        end

    end
    
    return true
end # twisting_circle

function twisting_circle_frames()
    # Parameters:
    E = 200000.0;
    nu = 0.3;# Poisson ratio
    rho = 7.8e-9;
    b = 0.6; h = 6.0;  # cross-sectional dimensions 
    radius = 120 # radius of the circle
    utol = 1e-3;
    maxit = 30
    rotmag = 2*pi
    # Select the number of elements around the circumference.
    nel = 64;
    
    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> h, s -> b, s -> [0.0, 0.0, 1.0])

    tolerance = radius/nel/1000;
    fens, fes = frame_member([0 0 0; 2*pi 0 0], nel, cs)
    for i in 1:count(fens)
        a = fens.xyz[i, 1]
        fens.xyz[i, :] .= (radius+radius*cos(a), radius*sin(a), 0)
    end
    fens, fes = mergenodes(fens, fes, tolerance, [1, nel+1])
    
    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, 0.0, E, nu, 0.0)
    
    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))
    
    # Locate some nodes
    box = fill(0.0, 6)
    # The pinned end is
    initbox!(box, [0.0,0.0,0.0])
    clampedn = selectnode(fens; box = box, inflate = tolerance)
    # The force is applied at the elbow
    initbox!(box, [2*radius,0,0])
    torquen = selectnode(fens; box = box, inflate = tolerance)

    # Apply EBC's
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, clampedn, true, i)
    end
    for i in [2, 3, 4, 5, 6]
        setebc!(dchi, torquen, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    
    # Additional fields
    u1 = deepcopy(u0)
    Rfield1 = deepcopy(Rfield0)
    rhs = gathersysvec(dchi);
    TMPv = deepcopy(rhs)
    utol = 1e-13*dchi.nfreedofs;
    
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    
    tbox = plot_space_box([[0 -radius -radius]; [2.5*radius radius radius]])
    tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    plots = cat(tbox,  tshape0; dims = 1)
    pl = render(plots; width = 800, height = 800)
    sleep(0.5)
    savefig(pl, "twisting_circle-$(0).png", width = 300, height = 300)

    load_parameter_delta = 0.002
    load_parameters = 0.0:load_parameter_delta:1;
    step = 0
    for load_parameter in load_parameters
        dchi.values[:] .= 0.0
        dchi.fixed_values[torquen[1], 4] = rotmag * load_parameter_delta
        applyebc!(dchi) # Apply boundary conditions
        u1.values[:] = u0.values[:]; # guess
        Rfield1.values[:] = Rfield0.values[:]; # guess
        update_rotation_field!(Rfield1, dchi)

        println("Load: $load_parameter")
        iter = 1;
        while true
            Fr = restoringforce(femm, geom0, u1, Rfield1, dchi);       # Internal forces
            rhs = Fr;
            K = stiffness(femm, geom0, u1, Rfield1, dchi) + geostiffness(femm, geom0, u1, Rfield1, dchi);
            dchi = scattersysvec!(dchi, (K)\rhs); # Disp. incr
            dchi.values[torquen[1], 4] = 0.0
            u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
            update_rotation_field!(Rfield1, dchi)
            print("$iter: ||du||=$(maximum(abs.(dchi.values[:])))\n")
            if maximum(abs.(dchi.values[:])) < utol # convergence check
                break;
            end
            if (iter > maxit)# bailout for failed convergence
                error("Possible failed convergence");
            end
            iter += 1;
        end
        u0.values[:] = u1.values[:];       # update the displacement
        Rfield0.values[:] = Rfield1.values[:]; # update the rotations
        
        
        step = step + 1
        if mod(step, 5) == 1
            @show step
            tshape1 = plot_solid(fens, fes; x = geom0.values, u = u1.values, R = Rfield1.values, facecolor = "rgb(125, 15, 15)");
            plots = cat(tbox,  tshape0,  tshape1; dims = 1)
            pl.plot.layout[:title] = "Angle=$(num3dig(rotmag * load_parameter / pi))*pi"
            react!(pl, plots, pl.plot.layout)
            sleep(0.5)
            savefig(pl, "twisting_circle-$(step).png", width = 300, height = 300)
        end

    end
   
    return true
end # twisting_circle_frames

function twisting_circle_frames_alt()
    # Parameters:
    E = 200000.0;
    nu = 0.3;# Poisson ratio
    rho = 7.8e-9;
    b = 0.6; h = 6.0;  # cross-sectional dimensions 
    radius = 120 # radius of the circle
    utol = 1e-3;
    maxit = 30
    rotmag = 2*pi
    # Select the number of elements around the circumference.
    nel = 64;
    
    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> h, s -> b, s -> [0.0, 0.0, 1.0])

    tolerance = radius/nel/1000;
    fens, fes = frame_member([0 0 0; 2*pi 0 0], nel, cs)
    for i in 1:count(fens)
        a = fens.xyz[i, 1]
        fens.xyz[i, :] .= (radius+radius*cos(a), radius*sin(a), 0)
    end
    fens, fes = mergenodes(fens, fes, tolerance, [1, nel+1])
    
    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, 0.0, E, nu, 0.0)
    
    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))
    
    # Locate some nodes
    box = fill(0.0, 6)
    # The pinned end is
    initbox!(box, [0.0,0.0,0.0])
    clampedn = selectnode(fens; box = box, inflate = tolerance)
    # The force is applied at the elbow
    initbox!(box, [2*radius,0,0])
    torquen = selectnode(fens; box = box, inflate = tolerance)

    # Apply EBC's
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, clampedn, true, i)
    end
    for i in [2, 3, 4, 5, 6]
        setebc!(dchi, torquen, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    
    # Additional fields
    u1 = deepcopy(u0)
    Rfield1 = deepcopy(Rfield0)
    rhs = gathersysvec(dchi);
    TMPv = deepcopy(rhs)
    utol = 1e-13*dchi.nfreedofs;
    
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    
    tbox = plot_space_box([[0 -radius -radius]; [2.05*radius radius radius]])
    tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    plots = cat(tbox,  tshape0; dims = 1)

    layout = Layout(width="500", height="500", autosize=true, title="title",
        scene=attr(
            xaxis=attr(title="X", gridcolor="rgb(255, 255, 255)",
              zerolinecolor="rgb(255, 255, 255)",
              showbackground=true,
              backgroundcolor="rgb(255, 255, 255)"),
            yaxis=attr(title="Y", gridcolor="rgb(255, 255, 255)",
             zerolinecolor="rgb(255, 255, 255)",
             showbackground=true,
             backgroundcolor="rgb(255, 255, 255)"),
            zaxis=attr(title="Z", gridcolor="rgb(255, 255, 255)",
             zerolinecolor="rgb(255, 255, 255)",
             showbackground=true,
             backgroundcolor="rgb(255, 255, 255)"),
            camera = attr(eye = attr(x = 0.98, y = 0.9, z = 0.1), 
                center = attr(x = 0, y = 0, z = 0),
                up = attr(x = 0, y = 0, z = 1)),
            aspectmode = "manual"),
        showlegend=false,
        )
    pl = plot(plots, layout)
    display(pl)
    sleep(0.5)
    # savejson(pl, "plot.json")
    # return pl

    load_parameter_delta = 0.004
    load_parameters = 0.0:load_parameter_delta:1;
    step = 0
    for load_parameter in load_parameters
        dchi.values[:] .= 0.0
        dchi.fixed_values[torquen[1], 4] = rotmag * load_parameter_delta
        applyebc!(dchi) # Apply boundary conditions
        u1.values[:] = u0.values[:]; # guess
        Rfield1.values[:] = Rfield0.values[:]; # guess
        update_rotation_field!(Rfield1, dchi)

        println("Load: $load_parameter")
        iter = 1;
        while true
            Fr = restoringforce(femm, geom0, u1, Rfield1, dchi);       # Internal forces
            rhs = Fr;
            K = stiffness(femm, geom0, u1, Rfield1, dchi) + geostiffness(femm, geom0, u1, Rfield1, dchi);
            dchi = scattersysvec!(dchi, (K)\rhs); # Disp. incr
            dchi.values[torquen[1], 4] = 0.0
            u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
            update_rotation_field!(Rfield1, dchi)
            print("$iter: ||du||=$(maximum(abs.(dchi.values[:])))\n")
            if maximum(abs.(dchi.values[:])) < utol # convergence check
                break;
            end
            if (iter > maxit)# bailout for failed convergence
                error("Possible failed convergence");
            end
            iter += 1;
        end
        u0.values[:] = u1.values[:];       # update the displacement
        Rfield0.values[:] = Rfield1.values[:]; # update the rotations
        
        step = step + 1
        if mod(step, 5) == 1
            @show step
            tshape1 = plot_solid(fens, fes; x = geom0.values, u = u1.values, R = Rfield1.values, facecolor = "rgb(125, 15, 15)");
            plots = cat(tbox,  tshape0,  tshape1; dims = 1)
            pl.plot.layout[:title] = "Angle=$(num3dig(rotmag * load_parameter / pi))*pi"
            react!(pl, plots, pl.plot.layout)
            sleep(0.5)
            savefig(pl, "twisting_circle-$(step).png", width = 300, height = 300)
        end

    end
    
    return true
end # twisting_circle_frames_alt

function allrun()
    println("#####################################################")
    println("# twisting_circle ")
    twisting_circle()
    println("#####################################################")
    println("# twisting_circle_frames ")
    twisting_circle_frames()
    println("#####################################################")
    println("# twisting_circle_frames_alt ")
    twisting_circle_frames_alt()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
