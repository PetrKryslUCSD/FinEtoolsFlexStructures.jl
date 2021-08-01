"""
The initially twisted cantilever beam is one of the standard test
problems for verifying the finite-element accuracy [1]. The beam is
clamped at one end and loaded either with unit in-plane or 
unit out-of-plane force at the other. The centroidal axis of the beam is
straight at the undeformed  configuration, while its cross-sections are
twisted about the centroidal axis from 0 at the clamped end to pi/2 at
the free end. 

Reference:
Zupan D, Saje M (2004) On "A proposed standard set of problems to test
finite element accuracy": the twisted beam. Finite Elements in Analysis
and Design 40: 1445-1451.  

#     Loading in the Z direction
dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
#     Loading in the Y direction
dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
"""
module twisted_beam_examples

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellDSG3Module
using FinEtoolsFlexStructures.FEMMShellCSDSG3Module
using FinEtoolsFlexStructures.FEMMShellIsoPModule
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json, plot_triads

params_thicker_dir_3 = (t =  0.32, force = 1.0, dir = 3, uex = 0.005424534868469); 
params_thicker_dir_2 = (t =  0.32, force = 1.0, dir = 2, uex = 0.001753248285256); 

params_thinner_dir_3 = (t =  0.0032, force = 1.0e-6, dir = 3, uex = 0.005256); 
params_thinner_dir_2 = (t =  0.0032, force = 1.0e-6, dir = 2, uex = 0.001294); 


function test_dsg3(t = 0.32, force = 1.0, dir = 3, uex = 0.005424534868469, nL = 2, nW = 1, visualize = true)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.0;
    
    tolerance = W/nW/10
    fens, fes = T3block(L,W,nL,nW,:a);
    # conn = connasarray(fes)
    # for e in 1:size(conn, 1)
    #     c = rand(0:1:2)
    #     if c == 1
    #         conn[e, :] = conn[e, [2, 3, 1]]
    #     elseif c == 2
    #         conn[e, :] = conn[e, [3, 1, 2]]
    #     end
    # end
    # fes = fromarray!(fes, conn)
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i,1]/L*(pi/2); y=fens.xyz[i,2]-(W/2); z=fens.xyz[i,3];
        fens.xyz[i,:]=[fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    formul = FEMMShellDSG3Module
    femm = formul.FEMMShellDSG3(IntegDomain(fes, TriRule(1), t), mater)
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped end]

    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in 1:6
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[L L 0 0 0 0], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    v = FFlt[0, 0, 0, 0, 0, 0]
    v[dir] = force
    fi = ForceIntensity(v);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    @show dchi.values[nl, dir][1]/uex*100

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_csdsg3(t = 0.32, force = 1.0, dir = 3, uex = 0.005424534868469, nL = 32, nW = 4, visualize = true)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.0;
    
    tolerance = W/nW/10
    fens, fes = T3block(L,W,nL,2*nW,:a);
    # conn = connasarray(fes)
    # for e in 1:size(conn, 1)
    #     c = rand(0:1:2)
    #     if c == 1
    #         conn[e, :] = conn[e, [2, 3, 1]]
    #     elseif c == 2
    #         conn[e, :] = conn[e, [3, 1, 2]]
    #     end
    # end
    # fes = fromarray!(fes, conn)
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i,1]/L*(pi/2); y=fens.xyz[i,2]-(W/2); z=fens.xyz[i,3];
        fens.xyz[i,:]=[fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    formul = FEMMShellCSDSG3Module
    femm = formul.FEMMShellCSDSG3(IntegDomain(fes, TriRule(1), t), mater)
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped end]

    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in 1:6
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[L L 0 0 0 0], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    v = FFlt[0, 0, 0, 0, 0, 0]
    v[dir] = force
    fi = ForceIntensity(v);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    @show dchi.values[nl, dir][1]/uex*100

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        plot_nodes(fens),
        plot_triads(fens; triad_length = 0.14, x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_t3(t = 0.32, force = 1.0, dir = 3, uex = 0.005424534868469, n = 2, visualize = true)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.0;
    
    tolerance = W/n/100
    fens, fes = T6block(L,W,6*n,2*n,:a);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i,1]/L*(pi/2); y=fens.xyz[i,2]-(W/2); z=fens.xyz[i,3];
        fens.xyz[i,:]=[fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    formul = FEMMShellIsoPModule
    femm = formul.make(IntegDomain(fes, TriRule(3), t), mater)
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped end
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in 1:6
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[L L 0 0 0 0], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    v = FFlt[0, 0, 0, 0, 0, 0]
    v[dir] = force
    fi = ForceIntensity(v);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    @show dchi.values[nl, dir][1]/uex*100

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_t6(t = 0.32, force = 1.0, dir = 3, uex = 0.005424534868469, n = 2, visualize = true)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.0;
    
    tolerance = W/n/100
    fens, fes = T6block(L,W,6*n,2*n,:a);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i,1]/L*(pi/2); y=fens.xyz[i,2]-(W/2); z=fens.xyz[i,3];
        fens.xyz[i,:]=[fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    formul = FEMMShellIsoPModule
    femm = formul.make(IntegDomain(fes, TriRule(3), t), mater)
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped end
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in 1:6
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[L L 0 0 0 0], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    v = FFlt[0, 0, 0, 0, 0, 0]
    v[dir] = force
    fi = ForceIntensity(v);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    @show dchi.values[nl, dir][1]/uex*100

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_q4sri(n = 2, visualize = true)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.0;
    t =  0.32;
    force = 1.0; 9
    #     Loading in the Z direction
    dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
    #     Loading in the Y direction
    # dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;

    tolerance = W/n/100
    fens, fes = Q4block(L,W,2*n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i,1]/L*(pi/2); y=fens.xyz[i,2]-(W/2); z=fens.xyz[i,3];
        fens.xyz[i,:]=[fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    formul = FEMMShellQ4SRIModule
    femm = formul.FEMMShellQ4SRI(IntegDomain(fes, GaussRule(2, 2), t), mater)
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped end
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in 1:6
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[L L 0 0 0 0], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    v = FFlt[0, 0, 0, 0, 0, 0]
    v[dir] = force
    fi = ForceIntensity(v);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    @show dchi.values[nl, dir][1]/uex*100

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 -L/4 -L/4]; [L L/2 L/2]]),
        plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_dsg3_convergence()
    for n in [2, 4, 8, 16, ]
        test_dsg3(params_thicker_dir_3..., 6*n, n, false)
    end
    return true
end

function test_csdsg3_convergence()
    for n in [2, 4, 8, 16]
        test_csdsg3(params_thicker_dir_3..., 6*n, n, false)
    end
    return true
end

function test_t3_convergence()
    for n in [2, 4, 8, 16, ]
        test_t3(params_thicker_dir_3..., n, false)
    end
    return true
end

function test_t6_convergence()
    for n in [2, 4, 8, 16, ]
        test_t6(params_thicker_dir_2..., n, false)
    end
    return true
end

function test_q4sri_convergence()
    for n in [2, 4, 8, ]
        test_q4sri(n, false)
    end
    return true
end

function allrun()
    println("#####################################################")
    println("# test_dsg3 ")
    test_dsg3()
    println("#####################################################")
    println("# test_q4sri ")
    test_q4sri()
    println("#####################################################")
    println("# test_dsg3_convergence  ")
    test_dsg3_convergence()
    println("#####################################################")
    println("# test_csdsg3_convergence  ")
    test_csdsg3_convergence()
    println("#####################################################")
    println("# test_t3_convergence  ")
    test_t3_convergence()
    println("#####################################################")
    println("# test_q4sri_convergence  ")
    test_q4sri_convergence()
    return true
end # function allrun

end # module

using .twisted_beam_examples
# twisted_beam_examples.test_dsg3_convergence()
# twisted_beam_examples.test_csdsg3()
# twisted_beam_examples.test_csdsg3_convergence()
twisted_beam_examples.test_t6_convergence()
# twisted_beam_examples.test_q4sri_convergence()