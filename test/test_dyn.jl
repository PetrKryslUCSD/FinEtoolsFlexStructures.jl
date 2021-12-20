module mfasttoptest1
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
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
function test()
    # Parameters:
    E = 71240.0;#MPa
    nu = 0.31;# Poisson ratio
    rho = 2.7e-9;
    Width = 60; Height = 60; Length = 4*Width;; # cross-sectional dimensions 
    Mass = Width*Height*Length*rho;
    IL =  (1/12)*Mass*(Width^2+Height^2);
    Ib =  (1/12)*Mass*(Length^2+Height^2)+Mass*(Length/2)^2;
    Ih =  (1/12)*Mass*(Width^2+Length^2)+Mass*(Length/2)^2;
    Omega0 = 313*pi;
    R0 = rotmat3([0.05 0 0]);
    g = 9.81e3; # 9.81 meters per second squared
    q = [0,0,-g*Width*Height*rho];
    utol = 1e-3;
    maxit = 12
    dt = min(2*pi/norm(Omega0)/10, 0.005);
    tend = 0.8;
    ng = 1/2; nb = 1/4*(1/2+ng)^2;
    # Choose the mass formulation:
    mass_type=1;

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> Width, s -> Width, s -> [1.0, 0.0, 0.0])

    # Select the number of elements per leg.
    spin_vector = R0*[0, 0, Omega0];
    X=[0 0 0;    reshape(R0*[0,0,Length], 1, 3)];
    n=8;
    tolerance=Length/n/100;
    members = []
    push!(members, frame_member(X, n, cs))
    fens, fes = merge_members(members; tolerance = tolerance);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    box = fill(0.0, 6)
    initbox!(box, X[1,:])
    supportn = selectnode(fens; box = box, inflate = tolerance)
    initbox!(box, X[2,:])
    tipn = selectnode(fens; box = box, inflate = tolerance)
    for i in [1, 2, 3]
        setebc!(dchi, supportn, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Additional fields
    stepdchi = deepcopy(dchi)
    u1 = deepcopy(u0)
    v0 = deepcopy(dchi) # need the numbering of the degrees of freedom
    for i in 1:size(v0.values, 1)
        v0.values[i, 4:6] = spin_vector
    end
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
    utol = 1e-13*dchi.nfreedofs;

    # tbox = plot_space_box([[-1.1*Width -1.1*Width 0]; [1.1*Width 1.1*Width 1.1*Length]])
    # tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    # plots = cat(tbox,  tshape0; dims = 1)
    # pl = render(plots)
    # sleep(3.5)

    tipx = Float64[]
    tipy = Float64[]
    push!(tipx, X[2,1])
    push!(tipy, X[2,2])
    # tbox = scatter(;x=[0.0, 0.06], y=[-0.06, 0.02], mode="markers")
    # plots = cat(tbox, scatter(;x=tipx./Length, y=tipy./Length, mode="markers+lines"); dims = 1)
    # layout = Layout(width=500, height=500)
    # pl = plot(plots, layout)
    # display(pl)

    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    fi = ForceIntensity(q);

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

        iter = 1;
        while true
            F = distribloads_global(femm, geom0, u1, Rfield1, dchi, fi)
            Fr = restoringforce(femm, geom0, u1, Rfield1, dchi);       # Internal forces
            @. rhs = F + Fr;
            K = stiffness(femm, geom0, u1, Rfield1, dchi);
            M = mass(femm, geom0, u1, Rfield1, dchi);
            G = gyroscopic(femm, geom0, u1, Rfield1, v1, dchi);
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
            if maximum(abs.(dchi.values[:])) < utol# convergence check
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

        # tshape1 = plot_solid(fens, fes; x = geom0.values, u = u1.values, R = Rfield1.values, facecolor = "rgb(125, 15, 15)");
        # plots = cat(tbox,  tshape0,  tshape1; dims = 1)
        # react!(pl, plots, pl.plot.layout)

        if (mod(step,50)==0)
            push!(tipx, X[2,1]+u1.values[tipn[1], 1])
            push!(tipy, X[2,2]+u1.values[tipn[1], 2])
            # plots = cat(tbox, scatter(;x=tipx./Length, y=tipy./Length, mode="markers+lines"); dims = 1)
            # react!(pl, plots, pl.plot.layout)
            # sleep(0.01)
        end

        step=step+1;
    end

    # # Visualize vibration modes
    # scattersysvec!(dchi, v[:, 1])
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_nodes(fens),
    #     plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    #     dims = 1)
    # render(plots; aspectratio = space_aspectratio(fens.xyz))
    reftipx = [0.0, 7.436031085824281e-7, 0.11777019275816547, 0.8203709506659435, 2.241864039281919, 3.9679895248072063, 5.3454217754301885, 5.971372979336259, 5.994274006124595, 5.998838219173172, 6.562336269473291, 7.816659235667364, 9.346281101346545, 10.499101141896105, 10.875849911112233, 10.628929261212697, 10.345322488144156, 10.602025422687818, 11.529299245601559, 12.710964257026657, 13.496913529329909, 13.492271339298975, 12.854269845420365, 12.172880227196368, 12.025322964890227, 12.540504127222404, 13.301744975038586]                                                                              
    reftipy = [-11.995000624962799, -11.99507329138783, -12.34434321937684, -13.092699071921164, -13.606709453527161, -13.384277801594491, -12.433899377302394, -11.269034682594993, -10.521805428897146, -10.467440913151172, -10.818091285108137, -10.941387824929802, -10.338997425035933, -9.024387391352876, -7.51482956095552, -6.443555337060203, -6.084606620045228, -6.148711597765479, -6.004002866512018, -5.155146068191021, -3.6197117072489498, -1.9179359415572481, -0.6834506841073384, -0.18854220120395482, -0.1419789099985156, 0.08825549466116911, 0.9953416189668296]   
    @test norm(reftipx - tipx) / norm(reftipx) <= 1.0e-3
    @test norm(reftipy - tipy) / norm(reftipy) <= 1.0e-3
    return true
end # test
end
using .mfasttoptest1
mfasttoptest1.test()

