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

Thicker cross section (t = 0.32)
#     Loading in the Z direction
dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
#     Loading in the Y direction
dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;

Thinner cross section (t = 0.0032)
TO DO What is the reference for the thinner section?
#     Loading in the Z direction
dir = 3; uex = 0.005256;
#     Loading in the Y direction
dir = 2; uex = 0.001294;

MacNeal,  R. H., and R. L. Harder, “A Proposed Standard Set of Problems to Test Finite Element Accuracy,” Finite Elements in Analysis Design, vol. 11, pp. 3–20, 1985.

Simo,  J. C., D. D. Fox, and M. S. Rifai, “On a Stress Resultant Geometrically Exact Shell Model. Part II: The Linear Theory; Computational Aspects,” Computational Methods in Applied Mechanical Engineering, vol. 73, pp. 53–92, 1989.
"""
module twisted_beam_examples

using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, DOF_KIND_FREE, DOF_KIND_ALL
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.FEMMShellQ4RSModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json, plot_triads
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

params_thicker_dir_3 = (t =  0.32, force = 1.0, dir = 3, uex = 0.005424534868469); 
params_thicker_dir_2 = (t =  0.32, force = 1.0, dir = 2, uex = 0.001753248285256); 

params_thinner_dir_3 = (t =  0.0032, force = 1.0e-6, dir = 3, uex = 0.005256); 
params_thinner_dir_2 = (t =  0.0032, force = 1.0e-6, dir = 2, uex = 0.001294); 


function _execute_t3ff(t = 0.32, force = 1.0, dir = 3, uex = 0.005424534868469, nL = 24, nW = 2, visualize = true)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.0;
    formul = FEMMShellT3FFModule
    
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
    femm.drilling_stiffness_scale = 0.1
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
    v = Float64[0, 0, 0, 0, 0, 0]
    v[dir] = force
    fi = ForceIntensity(v);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve
    solve_blocked!(dchi, K, F)
    result = dchi.values[nl, dir][1]
    @show result/uex*100

    # Visualization
    if visualize
        vtkwrite("twisted-geometry.vtu", fens, fes)
        scattersysvec!(dchi, (L/4)/dchi.values[nl, dir][1].*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
        dims = 1)
        pl = render(plots)
    end
    
    return result
end

function _execute_Q4RS(t = 0.32, force = 1.0, dir = 3, uex = 0.005424534868469, nL = 24, nW = 2, visualize = true)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.0;
    formul = FEMMShellQ4RSModule
    
    tolerance = W/nW/100
    fens, fes = Q4block(L,W,nL,nW);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i,1]/L*(pi/2); y=fens.xyz[i,2]-(W/2); z=fens.xyz[i,3];
        fens.xyz[i,:]=[fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, GaussRule(2, 2), t), mater)
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
    v = Float64[0, 0, 0, 0, 0, 0]
    v[dir] = force
    fi = ForceIntensity(v);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve
    solve_blocked!(dchi, K, F)
    result = dchi.values[nl, dir][1]
    @show result/uex*100

    # Visualization
    if visualize
        vtkwrite("twisted-geometry.vtu", fens, fes)
        vtkwrite("twisted-$(n)-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        # scattersysvec!(dchi, (L/4)/dchi.values[nl, dir][1].*U)
        # update_rotation_field!(Rfield0, dchi)
        # plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        # plot_nodes(fens),
        # plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
        # dims = 1)
        # pl = render(plots)
    end
    
    return result
end

function test_convergence_t3ff(ns = [2, 4, 8, 16, 32, ])
    @info "Twisted, thicker, T3FF elements"
    results = []
    for n in ns
        v = _execute_t3ff(params_thicker_dir_2..., 12*n, n, false)
        push!(results, v)
    end
    @show results
    results = []
    for n in ns
        v = _execute_t3ff(params_thicker_dir_3..., 12*n, n, false)
        push!(results, v)
    end
    @show results
    @info "Twisted, thinner, T3FF elements"
    results = []
    for n in ns
        v = _execute_t3ff(params_thinner_dir_2..., 12*n, n, false)
        push!(results, v)
    end
    @show results
    results = []
    for n in ns
        v = _execute_t3ff(params_thinner_dir_3..., 12*n, n, false)
        push!(results, v)
    end
    @show results
    return true
end

function test_convergence_Q4RS(ns = [2, 4, 8, 16, 32, ])
    @info "Twisted, thicker, Q4RS elements"
    results = []
    for n in ns
        v = _execute_Q4RS(params_thicker_dir_2..., 12*n, n, false)
        push!(results, v)
    end
    @show results
    results = []
    for n in ns
        v = _execute_Q4RS(params_thicker_dir_3..., 12*n, n, false)
        push!(results, v)
    end
    @show results
    @info "Twisted, thinner, Q4RS elements"
    results = []
    for n in ns
        v = _execute_Q4RS(params_thinner_dir_2..., 12*n, n, false)
        push!(results, v)
    end
    @show results
    results = []
    for n in ns
        v = _execute_Q4RS(params_thinner_dir_3..., 12*n, n, false)
        push!(results, v)
    end
    @show results
    return true
end

function allrun()
    println("#####################################################")
    println("# test_convergence_t3ff ")
    test_convergence_t3ff()
    println("#####################################################")
    println("# test_convergence_Q4RS ")
    test_convergence_Q4RS()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
