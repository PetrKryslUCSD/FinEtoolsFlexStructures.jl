# Simply supported equilateral triangle plate with uniform distributed load
module simply_supp_tri_plate_udl_examples

using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite



function   cartesian!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, feid::FInt, qpid::FInt)
    csmatout[:, 1] .= (1.0, 0.0, 0.0)
    csmatout[:, 2] .= (0.0, 1.0, 0.0)
    csmatout[:, 3] .= (0.0, 0.0, 1.0)
    return csmatout
end

function _execute_full_model(formul, n = 2, tL_ratio = 0.1, visualize = true)
    E = 10.0e9
    nu = 0.2
    L = 50.0;
    thickness = L * tL_ratio;
    D = E*thickness^3/12/(1-nu^2)
    G = E / 2 / (1+nu)
    ks = 5/6
    q0 = 1000*sqrt(tL_ratio)
    a = sqrt(3.0) * L / 2
    # Teik-Cheng Lim (2016) Higher-order shear deformation of very thick simply
    # supported equilateral triangular plates under uniform load, Mechanics Based Design of
    # Structures and Machines, 44:4, 514-522, DOI: 10.1080/15397734.2015.1124784
    cpt_solution = q0 * a^4 / 4 / D * (4.0/27) * (4.0/9/16)
    fsdt_solution = q0 * a^4 / 64 / D * (4.0/27) * (4.0/9 + 68 * thickness^2 / 21 / (1 - nu) / a^2)
    # REDDY
    fsdt_solution = q0 * a^4 / 64 / D * (4.0/27) * (4.0/9 + 16 * D / (ks * G * thickness * a^2))

    tolerance = L/n/1000
    fens = FENodeSet([
        -a/3 -a/sqrt(3)
        2*a/3 0
        -a/3 +a/sqrt(3)
        ])
    fes = FESetT3(reshape([1 2 3], 1, 3))

    for _ in 1:n
        fens, fes = T3refine(fens, fes);
    end

    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    ocsys = CSys(3, 3, cartesian!)

    @info "Mesh: $n refinements, number of elements: $(count(fes))"

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
    # plane of symmetry perpendicular to X
    l1 = connectednodes(meshboundary(fes))
    for i in [1, 2, 3]
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
    fi = ForceIntensity(FFlt[0, 0, q0, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);

    # @infiltrate
    # Solve
    solve!(dchi, K, F)
    targetu =  maximum(dchi.values[:, 3])
    @info "Target: $(round(targetu, digits=8))"
    @info "CPT: $(round(targetu/cpt_solution, digits = 4)*100)%"
    @info "FSDT: $(round(targetu/fsdt_solution, digits = 4)*100)%"

    # Visualization
    # if visualize

    #     # Generate a graphical display of resultants
    #     # scalars = []
    #     # for nc in 1:3
    #     #     fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
    #     #     push!(scalars, ("m$nc", fld.values))
    #     #     fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
    #     #     push!(scalars, ("em$nc", fld.values))
    #     # end
    #     # vtkwrite("scordelis_lo_examples-$(n)-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    #     # scalars = []
    #     # for nc in 1:3
    #     #     fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
    #     #     push!(scalars, ("n$nc", fld.values))
    #     #     fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
    #     #     push!(scalars, ("en$nc", fld.values))
    #     # end
    #     # vtkwrite("scordelis_lo_examples-$(n)-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    #     scalars = []
    #     for nc in 1:2
    #         fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
    #         push!(scalars, ("q$nc", fld.values))
    #         fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
    #         push!(scalars, ("eq$nc", fld.values))
    #     end
    #     vtkwrite("simply_supp_square_plate_udl-quarter-$(n)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

    #     scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
    #     update_rotation_field!(Rfield0, dchi)
    #     plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
    #     #plot_nodes(fens),
    #         plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    #         dims = 1)
    #     pl = render(plots)
    # end
    return true
end

function test_convergence_full()
    formul = FEMMShellT3FFModule
    tL_ratio = 0.001
    @info "Simply supported square plate with uniform load,"
    @info "Thickness/length = $tL_ratio"
    for n in 2:7
    # for n in [2, 4, 8, 16, 32, 64]
        _execute_full_model(formul, n, tL_ratio, true)
    end
    return true
end

function test_convergence_full_thickness(tL_ratios = [0.1, 0.01, 0.001, ])
    formul = FEMMShellT3FFModule
    @info "Simply supported square plate with uniform load,"
    for n in 2:8
        for tL_ratio in tL_ratios
            @info "Thickness/length = $tL_ratio"
            _execute_full_model(formul, n, tL_ratio, false)
        end
    end
    return true
end

function allrun()
    println("#####################################################")
    println("# test_convergence_full ")
    test_convergence_full()
    println("#####################################################")
    println("# test_convergence_full_thickness ")
    test_convergence_full_thickness()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
