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
using FinEtoolsFlexStructures.FEMMShellDSG3IModule
using FinEtoolsFlexStructures.FEMMShellDSG3IFModule
using FinEtoolsFlexStructures.FEMMShellCSDSG3Module
using FinEtoolsFlexStructures.FEMMShellIsoPModule
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json, plot_triads

params_thicker_dir_3 = (t =  0.32, force = 1.0, dir = 3, uex = 0.005424534868469); 
params_thicker_dir_2 = (t =  0.32, force = 1.0, dir = 2, uex = 0.001753248285256); 

params_thinner_dir_3 = (t =  0.0032, force = 1.0e-6, dir = 3, uex = 0.005256); 
params_thinner_dir_2 = (t =  0.0032, force = 1.0e-6, dir = 2, uex = 0.001294); 


function test_dsg3if(args...)
    return _execute_dsg_model(FEMMShellDSG3IFModule, args...)
end

function test_dsg3i(args...)
  return _execute_dsg_model(FEMMShellDSG3IModule, args...)
end

function test_dsg3(args...)
  return _execute_dsg_model(FEMMShellDSG3Module, args...)
end

function _execute_dsg_model(formul, t = 0.32, force = 1.0, dir = 3, uex = 0.005424534868469, nL = 8, nW = 2, visualize = true)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.0;
    
    tolerance = W/nW/100
    fens, fes = T3block(L,W,nL,nW,:a);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i,1]/L*(pi/2); y=fens.xyz[i,2]-(W/2); z=fens.xyz[i,3];
        fens.xyz[i,:]=[fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), t), mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

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
    associategeometry!(femm, geom0)
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
    scattersysvec!(dchi, (L/4)/dchi.values[nl, dir][1].*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_t3(t = 0.32, force = 1.0, dir = 3, uex = 0.005424534868469, n = 8, visualize = true)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.0;
    
    tolerance = W/n/100
    fens, fes = T3block(L,W,6*n,2*n,:a);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i,1]/L*(pi/2); y=fens.xyz[i,2]-(W/2); z=fens.xyz[i,3];
        fens.xyz[i,:]=[fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    formul = FEMMShellIsoPModule
    femm = formul.make(IntegDomain(fes, TriRule(1), t), mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

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
    associategeometry!(femm, geom0)
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
    v = FinEtools.MeshExportModule.VTK
    v.vtkexportmesh("twisted-$n" * ".vtk", fens, fes;   vectors=[("u", deepcopy(dchi.values[:, 1:3]))], scalars=[("ux", deepcopy(dchi.values[:, 1])), ("uy", deepcopy(dchi.values[:, 2])), ("uz", deepcopy(dchi.values[:, 3])), ("rx", deepcopy(dchi.values[:, 4])), ("ry", deepcopy(dchi.values[:, 5])), ("rz", deepcopy(dchi.values[:, 6]))])
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

function test_convergence(t)
    for n in [2, 4, 8, 16, ]
        t(params_thicker_dir_2..., 2*n, n, false)
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
    
    return true
end # function allrun

end # module

using .twisted_beam_examples
m = twisted_beam_examples
m.test_dsg3if()
m.test_convergence(m.test_dsg3if)