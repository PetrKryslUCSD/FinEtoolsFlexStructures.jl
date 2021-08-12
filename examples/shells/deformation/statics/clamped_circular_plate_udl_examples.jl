# Clamped circular plate with uniform distributed load

module clamped_circular_plate_udl_examples

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
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

using Infiltrator

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

function _execute_dsg_model(formul, mesh_procedure = :q4_t3, n = 2, t_radius_ratio = 0.01, visualize = true)
    E = 200*phun("GPa");
    nu = 0.3;
    a = 1.0*phun("m");
    thickness = a * t_radius_ratio;
    D = E*thickness^3/12/(1-nu^2)
    q = 1.0e10*t_radius_ratio^3
    # analytical solution for the vertical deflection under the load
    analyt_sol = -q*a^4/64/D*(1+16/5/(1-nu)*thickness^2/a^2);

    tolerance = a/n/1000
    if mesh_procedure == :q4_t3
        fens, fes = Q4circlen(a, n);
        fens, fes = Q4toT3(fens, fes)
    elseif mesh_procedure == :t3_nice
        fens, fes = T3circlen(a, n);
    elseif mesh_procedure == :t3_poor
        fens, fes = T3circleseg(pi/2, a, n, n);
    else
        @error "Unknown mesh procedure"
    end
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    
    @info "Mesh: $n elements per side"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    bfes = meshboundary(fes)
    lx = selectelem(fens, bfes; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    ly = selectelem(fens, bfes; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    lc = setdiff(1:count(bfes), vcat(lx, ly))
    # plane of symmetry perpendicular to X
    l1 = connectednodes(subset(bfes, lx))
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = connectednodes(subset(bfes, ly))
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    # clamped support
    l1 = connectednodes(subset(bfes, lc))
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], tolerance = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt[0, 0, -q, 0, 0, 0]);
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
    scattersysvec!(dchi, (a/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -a/2]; [a a a/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_convergence(t)
    t_radius_ratio = 0.00001
    @info "Simply supported square plate with uniform load,"
    @info "thickness/length = $t_radius_ratio formulation=$(t)"
    for n in [2, 4, 8, 16, 32, 64]
        t(:q4_t3, n, t_radius_ratio, false)
    end
    return true
end

end # module

using .clamped_circular_plate_udl_examples
m = clamped_circular_plate_udl_examples
m.test_convergence(m.test_dsg3)
m.test_convergence(m.test_dsg3i)
m.test_convergence(m.test_dsg3if)
# m.test_dsg3if()