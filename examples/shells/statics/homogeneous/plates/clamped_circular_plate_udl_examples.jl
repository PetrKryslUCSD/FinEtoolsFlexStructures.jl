"""
Clamped circular plate with uniform distributed load
"""
module clamped_circular_plate_udl_examples

using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function _execute(mesh_procedure = :q4_t3, n = 2, t_radius_ratio = 0.01, visualize = true)
    E = 200*phun("GPa");
    nu = 0.3;
    a = 1.0*phun("m");
    thickness = a * t_radius_ratio;
    D = E*thickness^3/12/(1-nu^2)
    q = 1.0e10*t_radius_ratio^3
    # analytical solution for the vertical deflection under the load
    analyt_sol = -q*a^4/64/D*(1+16/5/(1-nu)*thickness^2/a^2);
    # analyt_sol = -q*a^4/64/D
    formul = FEMMShellT3FFModule

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

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    Kff = stiffness(femm, massem, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], tolerance = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64[0, 0, -q, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 2);

    # @infiltrate
    # Solve
     Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)
    
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

function test_convergence()
    t_radius_ratio = 0.1
    @info "Simply supported square plate with uniform load,"
    @info "thickness/length = $t_radius_ratio "
    for n in [2, 4, 8, 16, 32, 64]
        _execute(:q4_t3, n, t_radius_ratio, false)
    end
    return true
end

function allrun()
    println("#####################################################")
    println("# test_convergence ")
    test_convergence()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
