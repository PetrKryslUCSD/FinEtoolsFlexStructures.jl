# Simply supported square plate with uniform distributed load
module LE5_Z_cantilever_examples

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellDSG3Module
using FinEtoolsFlexStructures.FEMMShellDSG3IModule
using FinEtoolsFlexStructures.FEMMShellCSDSG3Module
using FinEtoolsFlexStructures.FEMMShellT3Module
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

using Infiltrator

function test_dsg3(n = 8, tL_ratio = 0.01, visualize = true)
    E = 30e6;
    nu = 0.3;
    L = 10.0;
    thickness = L * tL_ratio;
    D = E*thickness^3/12/(1-nu^2)
    p = 1.0*tL_ratio
    # analytical solution for the vertical deflection under the load
    analyt_sol=-4.062e-3*p*L^4/D;

    tolerance = L/n/1000
    fens, fes = T3block(L/2,L/2,n,n);
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    formul = FEMMShellDSG3Module
    # Report
    @info "SS sq. plate with unif. distr. load, t/L=$(tL_ratio), formulation=$(formul)"
    @info "Mesh: $n elements per side"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = formul.stiffness

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
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt[0, 0, -p, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);

    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    targetu =  dchi.values[nl, 3][1]
    @info "Target: $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

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

function test_dsg3i(input = "nle5xf3c.inp", visualize = true)
    E = 210e9;
    nu = 0.3;
    L = 10.0;
    thickness = 0.1

    tolerance = thickness/1000
    output = import_ABAQUS(input)
    fens = output["fens"]
    fes = output["fesets"][1]

    connected = findunconnnodes(fens, fes);
    fens, new_numbering = compactnodes(fens, connected);
    fes = renumberconn!(fes, new_numbering);

    # fens, fes = T3refine(fens, fes)
    # fens, fes = T3refine(fens, fes)
    # fens, fes = T3refine(fens, fes)
    # fens, fes = T3refine(fens, fes)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    formul = FEMMShellT3Module
    # Report
    @info "LE5 Z-cantilever, formulation=$(formul)"
    @info "Mesh: $input"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        plot_nodes(fens),
        plot_midsurface(fens, fes);
    dims = 1)
    pl = render(plots)

    # Apply EBC's
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,]
        setebc!(dchi, l1, true, i)
    end
    
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    
    # Load
    nl = selectnode(fens; box = Float64[10.0 10.0 1.0 1.0 0 0], tolerance = tolerance)
    loadbdry1 = FESetP1(reshape(nl, 1, 1))
    lfemm1 = FEMMBase(IntegDomain(loadbdry1, PointRule()))
    fi1 = ForceIntensity(FFlt[0, 0, +0.6e6, 0, 0, 0]);
    nl = selectnode(fens; box = Float64[10.0 10.0 -1.0 -1.0 0 0], tolerance = tolerance)
    loadbdry2 = FESetP1(reshape(nl, 1, 1))
    lfemm2 = FEMMBase(IntegDomain(loadbdry2, PointRule()))
    fi2 = ForceIntensity(FFlt[0, 0, -0.6e6, 0, 0, 0]);
    F = distribloads(lfemm1, geom0, dchi, fi1, 3) + distribloads(lfemm2, geom0, dchi, fi2, 3);

    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    targetu =  minimum(dchi.values[:, 3]), maximum(dchi.values[:, 3])
    @info "Target: $(round.(targetu, digits=8))"

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

function test_csdsg3(n = 8, tL_ratio = 0.01, visualize = true)
    E = 30e6;
    nu = 0.3;
    L = 10.0;
    thickness = L * tL_ratio;
    D = E*thickness^3/12/(1-nu^2)
    p = 1.0*tL_ratio
    # analytical solution for the vertical deflection under the load
    analyt_sol=-4.062e-3*p*L^4/D;

    tolerance = L/n/1000
    fens, fes = T3block(L/2,L/2,n,n);
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    formul = FEMMShellCSDSG3Module
    # Report
    @info "SS sq. plate with unif. distr. load, t/L=$(tL_ratio), formulation=$(formul)"
    @info "Mesh: $n elements per side"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.FEMMShellCSDSG3(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = formul.stiffness

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
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt[0, 0, -p, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);

    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    targetu =  dchi.values[nl, 3][1]
    @info "Target: $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

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

function test_t3(n = 8, tL_ratio = 0.01, visualize = true)
    E = 30e6;
    nu = 0.3;
    L = 10.0;
    thickness = L * tL_ratio;
    D = E*thickness^3/12/(1-nu^2)
    p = 1.0*tL_ratio
    # analytical solution for the vertical deflection under the load
    analyt_sol=-4.062e-3*p*L^4/D;

    tolerance = L/n/1000
    fens, fes = T3block(L/2,L/2,n,n);
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    formul = FEMMShellT3Module
    # Report
    @info "SS sq. plate with unif. distr. load, t/L=$(tL_ratio), formulation=$(formul)"
    @info "Mesh: $n elements per side"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.FEMMShellT3(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = formul.stiffness

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
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt[0, 0, -p, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);

    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    targetu =  dchi.values[nl, 3][1]
    @info "Target: $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

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
        test_dsg3(n, 0.001, false)
    end
    return true
end

function test_dsg3i_convergence()
    for n in [2, 4, 8, 16, 32, 64]
        test_dsg3i(n, 0.001, false)
    end
    return true
end

function test_csdsg3_convergence()
    for n in [2, 4, 8, 16, 32, 64]
        test_csdsg3(n, 0.00001, false)
    end
    return true
end

function test_t3_convergence()
    for n in [2, 4, 8, 16, 32, 64]
        test_t3(n, 0.0001, false)
    end
    return true
end

function test_q4sri_convergence()
    for n in [2, 4, 8, 16, 32, 64]
        test_q4sri(n, 0.001, false)
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

using .LE5_Z_cantilever_examples
LE5_Z_cantilever_examples.test_dsg3i()
# LE5_Z_cantilever_examples.test_dsg3_convergence()
# LE5_Z_cantilever_examples.test_dsg3i_convergence()
# LE5_Z_cantilever_examples.test_csdsg3_convergence()
# LE5_Z_cantilever_examples.test_t3_convergence()
# LE5_Z_cantilever_examples.test_q4sri_convergence()