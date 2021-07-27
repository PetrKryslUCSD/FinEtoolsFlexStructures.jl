# Pinched cylinder with diagphram supports and concentrated force
module pinched_cylinder_examples

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3, local_frame!
using FinEtoolsFlexStructures.FEMMShellDSG3Module
using FinEtoolsFlexStructures.FEMMShellCSDSG3Module
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function test_dsg3(n = 4, visualize = true)
    E = 3e6;
    nu = 0.3;
    thickness = 3.0;
    # analytical solution for the vertical deflection under the load
    analyt_sol=-1.82488e-5;

    # Mesh
    R = 300.0;
    L = 600.0;
    
    tolerance = R/n/1000
    fens, fes = T3block(90/360*2*pi,L/2,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*cos(a))
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    formul = FEMMShellDSG3Module
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.FEMMShellDSG3(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # rigid diaphragm
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,3]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Z
    l1 = selectnode(fens; box = Float64[-Inf Inf -Inf Inf 0.0 0.0], inflate = tolerance)
    for i in [3,4,5]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 L/2 L/2 -Inf Inf], inflate = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, 0, -1/4, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    @show n, dchi.values[nl, 3],  dchi.values[nl, 3]/analyt_sol*100

    # Visualization
    scattersysvec!(dchi, (L/8)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
end


function test_csdsg3(n = 4, visualize = true)
    E = 3e6;
    nu = 0.3;
    thickness = 3.0;
    # analytical solution for the vertical deflection under the load
    analyt_sol=-1.82488e-5;

    # Mesh
    R = 300.0;
    L = 600.0;
    
    tolerance = R/n/1000
    fens, fes = T3block(90/360*2*pi,L/2,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*cos(a))
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    formul = FEMMShellCSDSG3Module
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
    # rigid diaphragm
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,3]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Z
    l1 = selectnode(fens; box = Float64[-Inf Inf -Inf Inf 0.0 0.0], inflate = tolerance)
    for i in [3,4,5]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 L/2 L/2 -Inf Inf], inflate = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, 0, -1/4, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    @show n, dchi.values[nl, 3],  dchi.values[nl, 3]/analyt_sol*100

    # Visualization
    scattersysvec!(dchi, (L/8)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
end

function test_dsg3_convergence()
    for n in [2, 4, 8, 16, 32, 64, 128]
        test_dsg3(n, false)
    end
    return true
end

function test_csdsg3_convergence()
    for n in [2, 4, 8, 16, 32, 64, 128]
        test_csdsg3(n, false)
    end
    return true
end

function test_q4sri_convergence()
    for n in [2, 4, 8, 16, 32, 64, 128]
        test_q4sri(n, false)
    end
    return true
end

function allrun()
    println("#####################################################")
    println("# test_dsg3 ")
    test_dsg3()
    println("#####################################################")
    println("# test_dsg3_convergence  ")
    test_dsg3_convergence()
    println("#####################################################")
    println("# test_q4sri_convergence  ")
    # test_q4sri_convergence()
    return true
end # function allrun

end # module

using .pinched_cylinder_examples
pinched_cylinder_examples.test_dsg3_convergence()
pinched_cylinder_examples.test_csdsg3_convergence()