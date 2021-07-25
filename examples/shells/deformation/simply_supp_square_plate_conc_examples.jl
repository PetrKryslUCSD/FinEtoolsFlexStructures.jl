# Simply supported square plate with center concentrated load
module simply_supp_square_plate_conc_examples

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellDSG3Module
# using FinEtoolsFlexStructures.FEMMShellT3Module
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function test_dsg3(n = 8, visualize = true)
    E = 30e6;
    nu = 0.3;
    force = 40;
    thickness = 0.1;
    L= 10;
    # analytical solution for the vertical deflection under the load
    analyt_sol=-0.0168;

    tolerance = L/n/1000
    fens, fes = T3block(L/2,L/2,n,n);
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellDSG3Module.FEMMShellDSG3(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = FEMMShellDSG3Module.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 0 L/2 -Inf Inf], inflate = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[0 L/2 0 0 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    # simple support
    l1 = selectnode(fens; box = Float64[L/2 L/2 0 L/2 -Inf Inf], inflate = tolerance)
    for i in [3,4,6]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 L/2 L/2 L/2 -Inf Inf], inflate = tolerance)
    for i in [3,5,6]
        setebc!(dchi, l1, true, i)
    end
    # in-plane, rotations
    l1 = selectnode(fens; box = Float64[0 L/2 0 L/2 -Inf Inf], inflate = tolerance)
    for i in [1, 2, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, 0, -force/4, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    @show dchi.values[nl, 3]/analyt_sol*100

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_q4sri(n = 8, visualize = true)
    E = 30e6;
    nu = 0.3;
    force = 40;
    thickness = 0.1;
    L= 10;
    # analytical solution for the vertical deflection under the load
    analyt_sol=-0.0168;

    tolerance = L/n/1000
    fens, fes = Q4block(L/2,L/2,n,n);
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = FEMMShellQ4SRIModule.FEMMShellQ4SRI(IntegDomain(fes, GaussRule(2, 2), thickness), mater)
    stiffness = FEMMShellQ4SRIModule.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 0 L/2 -Inf Inf], inflate = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[0 L/2 0 0 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    # simple support
    l1 = selectnode(fens; box = Float64[L/2 L/2 0 L/2 -Inf Inf], inflate = tolerance)
    for i in [3,4,6]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 L/2 L/2 L/2 -Inf Inf], inflate = tolerance)
    for i in [3,5,6]
        setebc!(dchi, l1, true, i)
    end
    # in-plane, rotations
    l1 = selectnode(fens; box = Float64[0 L/2 0 L/2 -Inf Inf], inflate = tolerance)
    for i in [1, 2, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, 0, -force/4, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    @show dchi.values[nl, 3]/analyt_sol*100

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_dsg3_convergence()
    for n in [2, 4, 8, 16, 32, 64]
        test_dsg3(n, false)
    end
    return true
end

function test_q4sri_convergence()
    for n in [2, 4, 8, 16, 32, 64]
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
    println("# test_q4sri_convergence  ")
    test_q4sri_convergence()
    return true
end # function allrun

end # module

using .simply_supp_square_plate_conc_examples
simply_supp_square_plate_conc_examples.test_dsg3_convergence()
simply_supp_square_plate_conc_examples.test_q4sri_convergence()