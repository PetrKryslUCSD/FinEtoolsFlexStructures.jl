"""
Square plate with uniform distributed load.
Two opposite sides simply supported, the other ones free.

FEMOL 
    w(nmidfreeedge) * (p*L^4/D/100) = -1.507
    q_2(A) * (p * L) = 4.311
"""
module sfsf_square_plate_udl_examples

using Arpack
using LinearAlgebra: norm
using Random: rand, Xoshiro
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.FEMMShellQ4RSModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule: VTKWrite
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using PGFPlotsX

const E = 30e6
const nu = 0.3
const L = 10.0

loading(tL_ratio) = 1.0e10 * (tL_ratio)^3

function _execute_q4rs_quarter_model(
    n=2,
    tL_ratio=0.01,
    simple_support=:hard,
    stab_alpha=0.1,
    mesh=:none,
    visualize=true
    )
    formul = FEMMShellQ4RSModule
    thickness = L * tL_ratio
    D = E / 12 / (1 - nu^2) * thickness^3
    p = loading(tL_ratio)
    
    if mesh == :uniform
        fens, fes = Q4block(L/2, L/2, n, n);
    elseif mesh == :biased
        xs = biasedspace(0.0, L/2, n+1, 100)
        ys = biasedspace(0.0, L/2, n+1, 100)
        fens, fes = Q4blockx(xs, ys);
    else
        xs = L/2 .* vcat(linearspace(0.0, tL_ratio, n), linearspace(tL_ratio, 1.0, n)[2:end])
        ys = L/2 .* vcat(linearspace(0.0, tL_ratio, n), linearspace(tL_ratio, 1.0, n)[2:end])
        fens, fes = Q4blockx(xs, ys);
    end
    bfes = meshboundary(fes)
    elleft = selectelem(fens, bfes; facing = true, direction = Float64[-1, 0])
    elright = selectelem(fens, bfes; facing = true, direction = Float64[+1, 0])
    eldown = selectelem(fens, bfes; facing = true, direction = Float64[0, -1])
    elup = selectelem(fens, bfes; facing = true, direction = Float64[0, +1])
    ncenter = selectnode(fens; box=Float64[(L/2) (L/2) (L/2) (L/2)], tolerance=eps(1.0))
    nmidfreeedge = selectnode(fens; box=Float64[(0.0) (0.0) (L/2) (L/2)], tolerance=eps(1.0))
    fens.xyz = xyz3(fens)
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    function _cartesian!(csmatout, XYZ, tangents, feid, qpid)
        csmatout[:, 1] .= (1.0, 0.0, 0.0)
        csmatout[:, 2] .= (0.0, 1.0, 0.0)
        csmatout[:, 3] .= (0.0, 0.0, 1.0)
        return csmatout
    end

    ocsys = CSys(3, 3, _cartesian!)

    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    ir = Simpson13Rule(2) # 
    ir = GaussRule(2, 2)
    femm = formul.make(IntegDomain(fes, ir, thickness),
        mater,
        (t, h) -> t^2 / (t^2 + stab_alpha * h^2),
    )
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # left
    dof = [1, 2, 6]
    for i in dof
        setebc!(dchi, connectednodes(subset(bfes, elleft)), true, i)
    end
    # right -- symmetry
    dof = [1, 2, 5, 6]
    for i in dof
        setebc!(dchi, connectednodes(subset(bfes, elright)), true, i)
    end
    # down
    dof = [1, 2, 3, 6]
    simple_support == :hard && (dof = [1, 2, 3, 5, 6])
    for i in dof
        setebc!(dchi, connectednodes(subset(bfes, eldown)), true, i)
    end
    # up -- symmetry
    dof = [1, 2, 4, 6]
    for i in dof
        setebc!(dchi, connectednodes(subset(bfes, elup)), true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # AE = AbaqusExporter("ss-qua-$(simple_support)-$(mesh)-tL=$(tL_ratio)-n=$(n).inp")
    # HEADING(AE, "Simply supported square plate with uniform distributed load")
    # PART(AE, "part1")
    # NODE(AE, fens.xyz)
    # @assert nodesperelem(fes) == 4
    # ELEMENT(
    #     AE,
    #     "S4",
    #     "AllElements",
    #     connasarray(fes),
    # )
    # NSET_NSET(AE, "edges_parallel_to_y", edges_parallel_to_y[1])
    # NSET_NSET(AE, "edges_parallel_to_x", edges_parallel_to_x[1])
    # NSET_NSET(AE, "symmetry_orthogonal_to_y", symmetry_orthogonal_to_y[1])
    # NSET_NSET(AE, "symmetry_orthogonal_to_x", symmetry_orthogonal_to_x[1])
    # SHELL_SECTION(AE, "elasticity", "AllElements", thickness)
    # END_PART(AE)
    # ASSEMBLY(AE, "ASSEM1")
    # ORIENTATION(AE, "GlobalOrientation", vec([1.0 0 0]), vec([0 1.0 0]))
    # INSTANCE(AE, "INSTNC1", "PART1")
    # END_INSTANCE(AE)
    # END_ASSEMBLY(AE)
    # MATERIAL(AE, "elasticity")
    # ELASTIC(AE, E, nu)
    # STEP_PERTURBATION_STATIC(AE)
    # DLOAD(AE, "ASSEM1.INSTNC1.AllElements", [0, 0, -p])   
    # for d in edges_parallel_to_y[2]
    #     BOUNDARY(AE, "ASSEM1.INSTNC1.edges_parallel_to_y", d)
    # end
    # for d in edges_parallel_to_x[2]
    #     BOUNDARY(AE, "ASSEM1.INSTNC1.edges_parallel_to_x", d)
    # end
    # for d in symmetry_orthogonal_to_y[2]
    #     BOUNDARY(AE, "ASSEM1.INSTNC1.symmetry_orthogonal_to_y", d)
    # end
    # for d in symmetry_orthogonal_to_x[2]
    #     BOUNDARY(AE, "ASSEM1.INSTNC1.symmetry_orthogonal_to_x", d)
    # end
    # FIELD_OUTPUT(AE, "PRESELECT", "LE, S, SF")
    # END_STEP(AE)
    # close(AE)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Load
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(Float64[0, 0, -p, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 2)

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    targetu = dchi.values[nmidfreeedge, 3][1]
    @info "w(nmidfreeedge): $(round(targetu, digits=8) / (p*L^4/D/100))"

    # Visualization
    if visualize
        vtkwrite("sfsfsqpludl-q4rs-$(simple_support)-$(mesh)-tL=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-uur.vtu", fens, fes; vectors=[("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        # ocsys = CSys(3)
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys=ocsys)
            push!(scalars, ("m$nc", fld.values))
        end
        vtkwrite("sfsfsqpludl-q4rs-$(simple_support)-$(mesh)-tL=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-m.vtu", fens, fes; 
            scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys=ocsys)
            push!(scalars, ("q$nc", fld.values))
            @info "q$nc Range: $(minimum(fld.values) / (p * L )) to $(maximum(fld.values) / (p * L ))"
        end
        vtkwrite("sfsfsqpludl-q4rs-$(simple_support)-$(mesh)-tL=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-q.vtu", fens, fes; scalars=scalars,
            vectors=[("u", dchi.values[:, 1:2])])

    end
    return targetu
end

function test_convergence_quarter_thickness(tL_ratios=[1/50])
    visualize = true
    ns = [16, 32, 64, 128]
    stab_alpha = 0.1
    for support in [:hard, :soft]
        @info "Support $support --------------------------------------------------"
        for mesh in [:biased]
            @info "Simply supported square plate with uniform load, Q4RS, stab_alpha=$stab_alpha  "
            @info "Mesh distortion: $mesh"
            for tL_ratio in tL_ratios
                @info "thickness/length = $tL_ratio"
                for n in ns
                    @info "n = $n"
                    _execute_q4rs_quarter_model(n, tL_ratio, support, stab_alpha, mesh, visualize)
                end

            end
        end
    end
    return true
end



function allrun()
    # println("#####################################################")
    # println("# test_stabilization_t3ff ")
    # test_stabilization_t3ff()
    # println("#####################################################")
    # println("# test_stabilization_q4rs ")
    # test_stabilization_q4rs()
     println("#####################################################")
    println("# test_convergence_quarter_thickness ")
    test_convergence_quarter_thickness()
    # println("#####################################################")
    # println("# test_convergence_full ")
    # test_convergence_full()
    # println("#####################################################")
    # println("# test_convergence_full_thickness ")
    # test_convergence_full_thickness()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
