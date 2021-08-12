# Pinched cylinder with diagphram supports and concentrated force
module pinched_cylinder_examples

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3, local_frame!
using FinEtoolsFlexStructures.FEMMShellDSG3Module
using FinEtoolsFlexStructures.FEMMShellDSG3IModule
using FinEtoolsFlexStructures.FEMMShellDSG3IFModule
using FinEtoolsFlexStructures.FEMMShellCSDSG3Module
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function test_dsg3if(args...)
    return _execute_dsg_model(FEMMShellDSG3IFModule, args...)
end

function test_dsg3i(args...)
  return _execute_dsg_model(FEMMShellDSG3IModule, args...)
end

function test_dsg3(args...)
  return _execute_dsg_model(FEMMShellDSG3Module, args...)
end

function test_csdsg3(args...)
  return _execute_dsg_model(FEMMShellCSDSG3Module, args...)
end

function _execute_dsg_model(formul, n = 2, visualize = true)
    E = 3e6;
    nu = 0.3;
    thickness = 3.0;
    # analytical solution for the vertical deflection under the load
    analyt_sol=-1.82488e-5;

    @info "Mesh: $n elements per side"

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

    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    associategeometry! = formul.associategeometry!
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
    associategeometry!(femm, geom0)
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
    @show n,  dchi.values[nl, 3]/analyt_sol*100

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (L/8)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_convergence(t)
    @info "Pinched cylinder, formulation=$(t)"
    for n in [8, 12, 16, 24]
        t(n, false)
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
m = pinched_cylinder_examples
m.test_convergence(m.test_dsg3if)
m.test_convergence(m.test_dsg3i)
m.test_convergence(m.test_dsg3)
