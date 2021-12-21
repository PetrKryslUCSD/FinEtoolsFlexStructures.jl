module simo_vuquoc_examples
# Simo, J. C. and L. Vuquoc (1988). "ON THE DYNAMICS IN SPACE OF RODS 
# UNDERGOING LARGE MOTIONS - A GEOMETRICALLY EXACT APPROACH." 
# Computer Methods in Applied Mechanics and Engineering 66(2): 125-161.

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
using Examples.VisUtilModule: plot_space_box, plot_midline, plot_solid, render, react!, default_layout_3d, save_to_json
using PlotlyJS
using JSON


Aus =[0  0
   9.4612e-01   2.5692e-02
   1.5375e+00   9.2688e-01
   2.3259e+00   2.3913e+00
   3.0355e+00   3.6680e+00
   3.8633e+00   5.0949e+00
   4.5729e+00   6.3715e+00
   5.1643e+00   7.2352e+00
   5.8739e+00   8.0237e+00
   6.5834e+00   8.3241e+00
   7.4507e+00   8.0237e+00
   8.2786e+00   6.7846e+00
   8.8305e+00   5.6957e+00
   9.2641e+00   4.4190e+00
   9.8949e+00   2.6917e+00
   1.0329e+01   1.1522e+00
   1.0723e+01  -5.3755e-01
   1.1275e+01  -3.4664e+00
   1.1787e+01  -5.5692e+00
   1.2260e+01  -7.0711e+00
   1.2970e+01  -7.9723e+00
   1.3758e+01  -8.5356e+00
   1.4428e+01  -9.2115e+00
   1.5059e+01  -9.2115e+00
   1.5966e+01  -8.4229e+00
   1.6873e+01  -6.6581e+00
   1.7740e+01  -4.3300e+00
   1.8568e+01  -1.8893e+00
   1.8962e+01  -8.6957e-02
   1.9632e+01   2.5791e+00
   2.0302e+01   4.4565e+00
   2.0894e+01   5.5079e+00
   2.1879e+01   6.5217e+00
   2.2746e+01   6.7846e+00
   2.3693e+01   6.2964e+00
   2.4481e+01   6.1462e+00
   2.5191e+01   6.6719e+00
   2.5821e+01   7.4980e+00
   2.6413e+01   6.9348e+00
   2.6965e+01   5.4704e+00
   2.7438e+01   4.1561e+00
   2.7792e+01   2.8419e+00
   2.8423e+01   1.2648e+00
   2.9409e+01   8.1423e-01
   2.9961e+01   7.3913e-01];


function simo_vuquoc_animated()
    # Parameters:
    # The  parameters given in the paper  are in terms of stiffness and mass coefficients. 
    asquared = 12.0e-3
    E = 1.0e6 / asquared;
    a = sqrt(asquared)
    E * 1/12 * a^4
    nu = 0.0;# Poisson ratio
    rho = 1 / asquared;
    b = a; h = a; L = 10; # cross-sectional dimensions and length of each leg
    utol = 1e-3;
    maxit = 30
    Fmag = 50;
    dt = 0.05;
    tend = 30.0;
    ng = 1/2; nb = 1/4*(1/2+ng)^2;
    # Choose the mass formulation:
    mass_type=1;

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0])

    # Select the number of elements per leg.
    n=4;
    tolerance=L/n/100;
    members = []
    push!(members, frame_member([0 L 0;L L 0], n, cs))
    push!(members, frame_member([L L 0;L 0 0], n, cs))
    fens, fes = merge_members(members; tolerance = tolerance);
    @show E * fes.A[1], E * fes.I2[1], E/2 * fes.J[1], rho * fes.A[1], rho * fes.I2[1]

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
    initbox!(box, [0,L,0])
    pinnedn = selectnode(fens; box = box, inflate = tolerance)
    # The force is applied at the elbow
    initbox!(box, [L,L,0])
    elbown = selectnode(fens; box = box, inflate = tolerance)
    # Tip is tracked
    initbox!(box, [L,0,0])
    tipn = selectnode(fens; box = box, inflate = tolerance)


    # Apply EBC's
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, pinnedn, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    dchi.dofnums

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
    utol = 1e-13*dchi.nfreedofs;

    tbox = plot_space_box([[-0.5*L -0.5*L -2.0*L]; [1.5*L 1.5*L 2.0*L]])
    tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    plots = cat(tbox,  tshape0; dims = 1)
    pl = render(plots)
    sleep(0.5)

    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    loadbdry = FESetP1(reshape(elbown, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))

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
        fi = ForceIntensity(FFlt[0, 0, (t<2.0)*((t<1.0)*t+(t>1.0)*(2.0-t))*Fmag, 0, 0, 0]);

        iter = 1;
        while true
            F = distribloads(lfemm, geom0, dchi, fi, 3);
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
end # simo_vuquoc_animated

function simo_vuquoc_compare()
    # Parameters:
    # The  parameters given in the paper  are in terms of stiffness and mass coefficients. 
    asquared = 12.0e-3
    E = 1.0e6 / asquared;
    a = sqrt(asquared)
    E * 1/12 * a^4
    nu = 0.0;# Poisson ratio
    rho = 1 / asquared;
    b = a; h = a; L = 10; # cross-sectional dimensions and length of each leg
    utol = 1e-3;
    maxit = 30
    Fmag = 50;
    dt = 0.05;
    tend = 30.0;
    ng = 1/2; nb = 1/4*(1/2+ng)^2;
    # Choose the mass formulation:
    mass_type=1;

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0])

    # Select the number of elements per leg.
    n=4;
    tolerance=L/n/100;
    members = []
    push!(members, frame_member([0 L 0;L L 0], n, cs))
    push!(members, frame_member([L L 0;L 0 0], n, cs))
    fens, fes = merge_members(members; tolerance = tolerance);
    @show E * fes.A[1], E * fes.I2[1], E/2 * fes.J[1], rho * fes.A[1], rho * fes.I2[1]

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
    initbox!(box, [0,L,0])
    pinnedn = selectnode(fens; box = box, inflate = tolerance)
    # The force is applied at the elbow
    initbox!(box, [L,L,0])
    elbown = selectnode(fens; box = box, inflate = tolerance)
    # Tip is tracked
    initbox!(box, [L,0,0])
    tipn = selectnode(fens; box = box, inflate = tolerance)


    # Apply EBC's
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, pinnedn, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    dchi.dofnums

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
    utol = 1e-13*dchi.nfreedofs;

    plots = cat(
            scatter(;x=Aus[:, 1], y=Aus[:, 2], mode="markers+lines", line_color = "rgb(255, 0, 0)", name="A"); 
            dims = 1)
        layout = Layout(width=700, height=500)
        pl = plot(plots, layout)
        display(pl)

    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    loadbdry = FESetP1(reshape(elbown, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))

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
        fi = ForceIntensity(FFlt[0, 0, (t<2.0)*((t<1.0)*t+(t>1.0)*(2.0-t))*Fmag, 0, 0, 0]);

        iter = 1;
        while true
            F = distribloads(lfemm, geom0, dchi, fi, 3);
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

        if (mod(step,10)==0)
            plots = cat(pl.plot.data, 
                scatter(;x=[t], y=[u1.values[elbown[1], 3]], mode="markers", line_color = "rgb(0, 0, 255)", showlegend=false); 
                dims = 1)
            react!(pl, plots, pl.plot.layout)
            sleep(0.1)
        end
        
        step=step+1;
    end

    return true
end # simo_vuquoc_compare

# function argyr_swing_compare()
#     # Parameters:
#     E = 71240.0;#MPa
#     nu = 0.31;# Poisson ratio
#     rho = 5.0e-9;
#     b=10; h=30; L=240; #cross-sectional dimensions and length of each leg in millimeters
#     utol = 1e-3;
#     maxit = 30
#     Fmag = 500;
#     dt = 0.001;
#     tend = 0.5;
#     ng = 1/2; nb = 1/4*(1/2+ng)^2;
#     # Choose the mass formulation:
#     mass_type=1;

#     # Cross-sectional properties
#     cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0])

#     # Select the number of elements per leg.
#     n=4;
#     tolerance=L/n/100;
#     members = []
#     push!(members, frame_member([0 L 0;L L 0], n, cs))
#     push!(members, frame_member([L L 0;L 0 0], n, cs))
#     fens, fes = merge_members(members; tolerance = tolerance);

#     # Material properties
#     material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

#     # Construct the requisite fields, geometry and displacement
#     # Initialize configuration variables
#     geom0 = NodalField(fens.xyz)
#     u0 = NodalField(zeros(size(fens.xyz,1), 3))
#     Rfield0 = initial_Rfield(fens)
#     dchi = NodalField(zeros(size(fens.xyz,1), 6))

#     # Locate some nodes
#     box = fill(0.0, 6)
#     # The pinned end is
#     initbox!(box, [0,L,0])
#     pinnedn = selectnode(fens; box = box, inflate = tolerance)
#     # The force is applied at the elbow
#     initbox!(box, [L,L,0])
#     elbown = selectnode(fens; box = box, inflate = tolerance)
#     # Tip is tracked
#     initbox!(box, [L,0,0])
#     tipn = selectnode(fens; box = box, inflate = tolerance)


#     # Apply EBC's
#     for i in [1, 2, 3]
#         setebc!(dchi, pinnedn, true, i)
#     end
#     applyebc!(dchi)
#     numberdofs!(dchi);
#     dchi.dofnums

#     # Additional fields
#     stepdchi = deepcopy(dchi)
#     u1 = deepcopy(u0)
#     v0 = deepcopy(dchi) # need the numbering of the degrees of freedom
#     v1 = deepcopy(dchi)
#     a1 = deepcopy(dchi) # zero out the Acceleration
#     a0 = deepcopy(dchi) # zero out the Acceleration
#     Rfield1 = deepcopy(Rfield0)
#     v0v = gathersysvec(v0)
#     a0v = gathersysvec(a0)
#     vpv = gathersysvec(v0);
#     dchipv = gathersysvec(dchi);
#     stepdchiv = gathersysvec(dchi);
#     rhs = gathersysvec(dchi);
#     TMPv = deepcopy(rhs)
#     utol = 1e-13*dchi.nfreedofs;

#     reft = range(0, 0.5, length = size(Aus,1))
#     plots = cat(
#         scatter(;x=reft, y=Aus[:,1], mode="lines", line_color = "rgb(255, 0, 0)", name="Tip-x"),
#         scatter(;x=reft, y=Aus[:,2], mode="lines", line_color = "rgb(0, 255, 0)", name="Tip-y"),
#         # scatter(;x=reft, y=Aus[:,3], mode="lines", line_color = "rgb(0, 0, 255)", name="Tip-z"),
#         dims = 1)
#     layout = Layout(width=700, height=500)
#     pl = plot(plots, layout)
#     display(pl)

#     femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
#     loadbdry = FESetP1(reshape(elbown, 1, 1))
#     lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))

#     t = 0.0; #
#     step = 0;
#     while (t <= tend)
#         t = t + dt;
#         # println("Time $(t)"); # pause
#         # Initialization
#         applyebc!(dchi) # Apply boundary conditions
#         u1.values[:] = u0.values[:]; # guess
#         Rfield1.values[:] = Rfield0.values[:]; # guess
#         stepdchi.values[:] .= 0.0;# Total increment in current step
#         a1.values[:] = -(1/nb/dt)*v0.values[:] -(1/2-nb)/nb*a0.values[:];
#         v1.values[:] = v0.values[:] + dt*((1-ng)*a0.values[:] + ng*a1.values[:]);
#         gathersysvec!(v0, v0v)
#         gathersysvec!(a0, a0v)
#         dchipv = dt*v0v + (dt^2/2*(1-2*nb))*a0v
#         vpv = v0v +(dt*(1-ng))*a0v;
#         fi = ForceIntensity(FFlt[0, 0, (t<0.2)*((t<0.1)*t+(t>0.1)*(0.2-t))*Fmag, 0, 0, 0]);

#         iter = 1;
#         while true
#             F = distribloads(lfemm, geom0, dchi, fi, 3);
#             Fr = restoringforce(femm, geom0, u1, Rfield1, dchi);       # Internal forces
#             @. rhs = F + Fr;
#             K = stiffness(femm, geom0, u1, Rfield1, dchi);
#             M = mass(femm, geom0, u1, Rfield1, dchi);
#             G = gyroscopic(femm, geom0, u1, Rfield1, v1, dchi);
#             gathersysvec!(stepdchi, stepdchiv)
#             @. TMPv = ((-1/(nb*dt^2))*stepdchiv+(1/(nb*dt^2))*dchipv)
#             rhs .+= M*TMPv
#             @. TMPv = ((-ng/nb/dt)*stepdchiv+(ng/nb/dt)*dchipv - vpv)
#             rhs .+= G*TMPv;
#             dchi = scattersysvec!(dchi, (K+(ng/nb/dt)*G+(1/(nb*dt^2))*M)\rhs); # Disp. incr
#             u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
#             stepdchi.values[:] += dchi.values[:]
#             v1.values[:] += (ng/nb/dt)*dchi.values[:];
#             a1.values[:] += (1/nb/dt^2)*dchi.values[:];
#             update_rotation_field!(Rfield1, dchi)
#             if maximum(abs.(dchi.values[:])) < utol # convergence check
#                 break;
#             end
#             if (iter > maxit)# bailout for failed convergence
#                 error("Possible failed convergence");
#             end
#             iter += 1;
#         end
#         u0.values[:] = u1.values[:];       # update the displacement
#         Rfield0.values[:] = Rfield1.values[:]; # update the rotations
#         v0.values[:] = v1.values[:];       # update the velocities
#         a0.values[:] = a1.values[:];       # update the accelerations

#         if (mod(step,10)==0)
#             plots = cat(pl.plot.data, 
#                 scatter(;x=[t], y=[u1.values[tipn[1], 1]], mode="markers", line_color = "rgb(255, 0, 0)", showlegend=false),
#                 scatter(;x=[t], y=[u1.values[tipn[1], 2]], mode="markers", line_color = "rgb(0, 255, 0)", showlegend=false),
#                 scatter(;x=[t], y=[u1.values[tipn[1], 3]], mode="markers", line_color = "rgb(0, 0, 255)", showlegend=false),
#                 scatter(;x=[t], y=[u1.values[elbown[1], 1]], mode="markers", line_color = "rgb(255, 0, 0)", showlegend=false),
#                 scatter(;x=[t], y=[u1.values[elbown[1], 2]], mode="markers", line_color = "rgb(0, 255, 0)", showlegend=false),
#                 scatter(;x=[t], y=[u1.values[elbown[1], 3]], mode="markers", line_color = "rgb(0, 0, 255)", showlegend=false); 
#                 dims = 1)
#             react!(pl, plots, pl.plot.layout)
#             sleep(0.1)
#         end

#         step=step+1;
#     end

#     return true
# end # argyr_swing_animated

function allrun()
    println("#####################################################")
    println("# simo_vuquoc_animated ")
    simo_vuquoc_animated()
    println("#####################################################")
    println("# simo_vuquoc_compare ")
    simo_vuquoc_compare()
    # println("#####################################################")
    # println("# argyr_swing_compare ")
    # argyr_swing_compare()
    return true
end # function allrun


@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
