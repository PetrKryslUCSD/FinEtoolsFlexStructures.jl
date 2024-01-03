module twisting_circle_examples
# Twisting of a circle into shape with two turns

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMCorotBeamModule
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
mass = FEMMCorotBeamModule.mass
distribloads_global = FEMMCorotBeamModule.distribloads_global
restoringforce = FEMMCorotBeamModule.restoringforce
gyroscopic = FEMMCorotBeamModule.gyroscopic
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using VisualStructures: plot_space_box, plot_midline, plot_solid, render, react!, default_layout_3d, save_to_json
using PlotlyJS
using JSON


function twisting_circle()
    # Parameters:
    # The  parameters given in the paper  are in terms of stiffness and mass coefficients. 
    E = 200000.0;
    nu = 0.3;# Poisson ratio
    rho = 7.8e-9;
    b = 0.6; h = 6.0;  # cross-sectional dimensions 
    radius = 120 # radius of the circle
    utol = 1e-3;
    maxit = 30
    Tmag = 700;
    dt = 0.0002;
    tend = 30.0;
    ng = 1/2; nb = 1/4*(1/2+ng)^2;
    # Select the number of elements around the circumference.
    n=32;
    # Choose the mass formulation:
    mass_type=1;

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> h, s -> b, s -> [0.0, 0.0, 1.0])

    @show tolerance = radius/n/1000;
    fens, fes = frame_member([0 0 0; 2*pi 0 0], n, cs)
    for i in 1:count(fens)
        a = fens.xyz[i, 1]
        fens.xyz[i, :] .= (radius+radius*cos(a), radius*sin(a), 0)
    end
    fens, fes = mergenodes(fens, fes, tolerance, [1, n+1])
    
    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

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
    

    # # Apply EBC's
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, clampedn, true, i)
    end
    for i in [2, 3, 5, 6]
        setebc!(dchi, torquen, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Additional fields
    stepdchi = deepcopy(dchi)
    u1 = deepcopy(u0)
    v0 = deepcopy(dchi) # need the numbering of the degrees of freedom
    v1 = deepcopy(dchi)
    a1 = deepcopy(dchi) # zero out the Acceleration
    a0 = deepcopy(dchi) # zero out the Acceleration
    Rfield1 = deepcopy(Rfield0)
    v0v = gathersysvec(v0)
    a0v = gathersysvec(a0)
    vpv = gathersysvec(v0);
    dchipv = gathersysvec(dchi);
    stepdchiv = gathersysvec(dchi);
    rhs = gathersysvec(dchi);
    TMPv = deepcopy(rhs)
    utol = 1e-13*nfreedofs(dchi);

    tbox = plot_space_box([[0 -radius -radius]; [2*radius radius radius]])
    tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    plots = cat(tbox,  tshape0; dims = 1)
    pl = render(plots)
    sleep(0.5)

    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    loadbdry = FESetP1(reshape(torquen, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    t = 0.0; #
    step = 0;
    while (t <= tend)
        t = t + dt;
        # println("Time $(t)"); # pause
        # Initialization
        applyebc!(dchi) # Apply boundary conditions
        u1.values[:] = u0.values[:]; # guess
        Rfield1.values[:] = Rfield0.values[:]; # guess
        stepdchi.values[:] .= 0.0;# Total increment in current step
        a1.values[:] = -(1/nb/dt)*v0.values[:] -(1/2-nb)/nb*a0.values[:];
        v1.values[:] = v0.values[:] + dt*((1-ng)*a0.values[:] + ng*a1.values[:]);
        gathersysvec!(v0, v0v)
        gathersysvec!(a0, a0v)
        dchipv = dt*v0v + (dt^2/2*(1-2*nb))*a0v
        vpv = v0v +(dt*(1-ng))*a0v;
        fi = ForceIntensity(Float64[0, 0, 0, (t < 0.08)*Tmag, 0, 0]);

        iter = 1;
        while true
            F = distribloads(lfemm, vassem, geom0, dchi, fi, 3);
            Fr = restoringforce(femm, vassem, geom0, u1, Rfield1, dchi);       # Internal forces
            @. rhs = F + Fr;
            K = stiffness(femm, massem, geom0, u1, Rfield1, dchi);
            M = mass(femm, massem, geom0, u1, Rfield1, dchi);
            G = gyroscopic(femm, massem, geom0, u1, Rfield1, v1, dchi);
            gathersysvec!(stepdchi, stepdchiv)
            @. TMPv = ((-1/(nb*dt^2))*stepdchiv+(1/(nb*dt^2))*dchipv)
            rhs .+= M*TMPv
            @. TMPv = ((-ng/nb/dt)*stepdchiv+(ng/nb/dt)*dchipv - vpv)
            rhs .+= G*TMPv;
            dchi = scattersysvec!(dchi, (K+(ng/nb/dt)*G+(1/(nb*dt^2))*M)\rhs); # Disp. incr
            u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
            stepdchi.values[:] += dchi.values[:]
            v1.values[:] += (ng/nb/dt)*dchi.values[:];
            a1.values[:] += (1/nb/dt^2)*dchi.values[:];
            update_rotation_field!(Rfield1, dchi)
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
        v0.values[:] = v1.values[:];       # update the velocities
        a0.values[:] = a1.values[:];       # update the accelerations

        if mod(step, 4) == 0
            @show step
            tshape1 = plot_solid(fens, fes; x = geom0.values, u = u1.values, R = Rfield1.values, facecolor = "rgb(125, 15, 15)");
            plots = cat(tbox,  tshape0,  tshape1; dims = 1)
            pl.plot.layout[:title] = "t=$t"
            react!(pl, plots, pl.plot.layout)
            sleep(0.5)
        end
        
        step=step+1;
    end

    return true
end # twisting_circle

function allrun()
    println("#####################################################")
    println("# twisting_circle ")
    twisting_circle()
    return true
end # function allrun


@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
