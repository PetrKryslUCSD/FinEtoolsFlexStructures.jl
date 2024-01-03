module curved_cantilever_examples
##  Curved (45 deg circular arc) cantilever, Transverse force at the tip
# 
# Large-deflection problem. Solved many times:
#
# Bathe, Boulourchi 1979
# Simo, Vu-Quoc 1986
# Cardona, Geradin 1988
# Krysl 1993
# Krysl, FAESOR script curved_cantilever
# 
# Present calculation refers to the data of
# Reference: EULERIAN FORMULATION FOR LARGE-DISPLACEMENT ANALYSIS OF SPACE
# FRAMES, by B.  A.  lzzuddin I and A.  S.  Elnashai,
# Journal of Engineering Mechanics, Vol.  119, No.  3,  March, 1993.


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
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using VisualStructures: plot_points, plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio
using PlotlyJS
using JSON


# Comparison with reference solution.

izzelna_ux =[    0.0317    0.6048
-0.0628    0.6048
8.7273   34.7095
16.6668   85.8664
23.7557  134.8918
32.6404  213.7587
37.8389  264.9157
41.7141  318.2041
44.9277  382.1503
47.9523  446.0965
50.4098  510.0426
51.9220  578.2519
53.9069  657.1188
53.9069  657.1188];

izzelna_uy =[-0.1574    0.6048
0.0317    2.7364
1.1659   73.0772
2.7727  132.7603
4.9466  211.6272
7.1205  284.0995
8.9164  350.1772
10.2396  414.1234
11.3738  480.2011
13.1697  548.4104
14.0203  618.7511
14.6820  657.1188
14.5875  657.1188];

izzelna_uz =[   -0.5354    2.7364
-0.2519    2.7364
-1.1970   51.7618
-2.2367   90.1295
-4.3161  124.2341
-6.8681  183.9172
-9.1366  235.0741
-11.8776  288.3626
-14.5241  354.4403
-16.9815  411.9919
-19.3445  475.9380
-21.7074  550.5419
-23.9759  618.7511
-25.3936  661.3819];


function curved_cantilever()
    # Parameters:
    E=1e7;# lb.in^-2
    nu=0.3;
    h= 1; b =1;# in
    Fmag=600;#  lb
    radius=100;# in
    ang=45;
    nel=8;# number of elements
    utol=1e-6;
    maxit = 150
    
    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0])
    
    # Select the number of elements per leg.
    fens, fes = frame_member([0 0 0; nel 0 0], nel, cs)
    for i in 1:count(fens)
        k = fens.xyz[i, 1]
        fens.xyz[i, :] .= (radius*cos(k/nel*ang*2*pi/360), radius*sin(k/nel*ang*2*pi/360), 0)
    end
    clampedn = [1]
    tipn = [nel+1];
    
    tolerance = radius/nel/100;
    
    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, 0.0, E, nu, 0.0)
    
    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))
    
    # Apply EBC's
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, clampedn, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    
    # Additional fields
    u1 = deepcopy(u0)
    Rfield1 = deepcopy(Rfield0)
    rhs = gathersysvec(dchi);
    TMPv = deepcopy(rhs)
    utol = 1e-13*nfreedofs(dchi);
    
    # tbox = plot_space_box([[0.0*radius 0.0*radius 0.0*radius]; [1.15*radius 0.75*radius 2.0*radius]])
    # tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    # plots = cat(tbox,  tshape0; dims = 1)
    # pl = render(plots)
    # sleep(0.5)
    
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    loadbdry = FESetP1(reshape(tipn, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    load_parameters = 0.0:0.05:1;
    tip_displacement = fill(0.0, length(load_parameters), 3);
    step = 0
    for load_parameter in load_parameters
        applyebc!(dchi) # Apply boundary conditions
        u1.values[:] = u0.values[:]; # guess
        Rfield1.values[:] = Rfield0.values[:]; # guess
        fi = ForceIntensity(Float64[0, 0, load_parameter*Fmag, 0, 0, 0]);
        
        println("Load: $load_parameter")
        iter = 1;
        while true
            F = distribloads(lfemm, vassem, geom0, dchi, fi, 3);
            Fr = restoringforce(femm, vassem, geom0, u1, Rfield1, dchi);       # Internal forces
            @. rhs = F + Fr;
            K = stiffness(femm, massem, geom0, u1, Rfield1, dchi) + geostiffness(femm, massem, geom0, u1, Rfield1, dchi);
            dchi = scattersysvec!(dchi, (K)\rhs); # Disp. incr
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
        tip_displacement[step, :] .= u1.values[tipn[1], :]

    end
    step = step + 1
    
    plots = cat(
    scatter(;x=tip_displacement[:,1], y=load_parameters.*Fmag, mode="lines", line_color = "rgb(0, 0, 0)", name="Tip-x"),
    scatter(;x=tip_displacement[:,2], y=load_parameters.*Fmag, mode="lines", line_color = "rgb(0, 0, 0)", name="Tip-y"),
    scatter(;x=tip_displacement[:,3], y=load_parameters.*Fmag, mode="lines", line_color = "rgb(0, 0, 0)", name="Tip-z"),
    scatter(;x=izzelna_ux[:,1], y=izzelna_ux[:,2], mode="markers", line_color = "rgb(255, 0, 0)", name="Ref-x"),
    scatter(;x=izzelna_uy[:,1], y=izzelna_uy[:,2], mode="markers", line_color = "rgb(0, 255, 0)", name="Ref-y"),
    scatter(;x=izzelna_uz[:,1], y=izzelna_uz[:,2], mode="markers", line_color = "rgb(0, 0, 255)", name="Ref-z"); 
    dims = 1)
    layout = Layout(width=700, height=500, xaxis=attr(title="Tip displacement [in]"),
    yaxis=attr(title="Load [lb]"))
    pl = plot(plots, layout)
    display(pl)
    
    return true
end # curved_cantilever

function curved_cantilever_thin()
    # Parameters:
    E=1e7;# lb.in^-2
    nu=0.3;
    thinness_factor = 100
    h= 1.0/thinness_factor; b = 1.0/thinness_factor;# in
    Fmag=600.0/(thinness_factor^4);#  lb
    radius=100;# in
    ang=45;
    nel=8;# number of elements
    utol=1e-6;
    maxit = 150
    
    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0])
    
    # Select the number of elements per leg.
    fens, fes = frame_member([0 0 0; nel 0 0], nel, cs)
    for i in 1:count(fens)
        k = fens.xyz[i, 1]
        fens.xyz[i, :] .= (radius*cos(k/nel*ang*2*pi/360), radius*sin(k/nel*ang*2*pi/360), 0)
    end
    clampedn = [1]
    tipn = [nel+1];
    
    tolerance = radius/nel/100;
    
    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, 0.0, E, nu, 0.0)
    
    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))
    
    # Apply EBC's
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, clampedn, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    
    # Additional fields
    u1 = deepcopy(u0)
    Rfield1 = deepcopy(Rfield0)
    rhs = gathersysvec(dchi);
    TMPv = deepcopy(rhs)
    utol = 1e-13*nfreedofs(dchi);
    
    # tbox = plot_space_box([[0.0*radius 0.0*radius 0.0*radius]; [1.15*radius 0.75*radius 2.0*radius]])
    # tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    # plots = cat(tbox,  tshape0; dims = 1)
    # pl = render(plots)
    # sleep(0.5)
    
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    loadbdry = FESetP1(reshape(tipn, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    
    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    load_parameters = 0.0:0.05:1;
    tip_displacement = fill(0.0, length(load_parameters), 3);
    step = 0
    for load_parameter in load_parameters
        applyebc!(dchi) # Apply boundary conditions
        u1.values[:] = u0.values[:]; # guess
        Rfield1.values[:] = Rfield0.values[:]; # guess
        fi = ForceIntensity(Float64[0, 0, load_parameter*Fmag, 0, 0, 0]);
        
        println("Load: $load_parameter")
        iter = 1;
        while true
            F = distribloads(lfemm, vassem, geom0, dchi, fi, 3);
            Fr = restoringforce(femm, vassem, geom0, u1, Rfield1, dchi);       # Internal forces
            @. rhs = F + Fr;
            K = stiffness(femm, massem, geom0, u1, Rfield1, dchi) + geostiffness(femm, massem, geom0, u1, Rfield1, dchi);
            dchi = scattersysvec!(dchi, (K)\rhs); # Disp. incr
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
        tip_displacement[step, :] .= u1.values[tipn[1], :]

    end
    step = step + 1
    
    plots = cat(
    scatter(;x=tip_displacement[:,1], y=load_parameters.*Fmag, mode="lines", line_color = "rgb(0, 0, 0)", name="Tip-x"),
    scatter(;x=tip_displacement[:,2], y=load_parameters.*Fmag, mode="lines", line_color = "rgb(0, 0, 0)", name="Tip-y"),
    scatter(;x=tip_displacement[:,3], y=load_parameters.*Fmag, mode="lines", line_color = "rgb(0, 0, 0)", name="Tip-z"),
    scatter(;x=izzelna_ux[:,1], y=izzelna_ux[:,2]/(thinness_factor^4), mode="markers", line_color = "rgb(255, 0, 0)", name="Ref-x"),
    scatter(;x=izzelna_uy[:,1], y=izzelna_uy[:,2]/(thinness_factor^4), mode="markers", line_color = "rgb(0, 255, 0)", name="Ref-y"),
    scatter(;x=izzelna_uz[:,1], y=izzelna_uz[:,2]/(thinness_factor^4), mode="markers", line_color = "rgb(0, 0, 255)", name="Ref-z"); 
    dims = 1)
    layout = Layout(width=700, height=500, xaxis=attr(title="Tip displacement [in]"),
    yaxis=attr(title="Load [lb]"))
    pl = plot(plots, layout)
    display(pl)
    
    return true
end # curved_cantilever_thin

function allrun()
    println("#####################################################")
    println("# curved_cantilever ")
    curved_cantilever()
    println("#####################################################")
    println("# curved_cantilever_thin ")
    curved_cantilever_thin()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
