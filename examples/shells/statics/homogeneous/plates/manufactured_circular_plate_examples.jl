"""
Manufactured solution for a circular plate with uniform distributed load.

Resultants are presented graphically.
"""
module manufactured_circular_plate_examples

using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtools.RotationUtilModule: cross3!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Targe2

C = -1/48
manufactured_solution_w(x, y) = C * (x^4 + y^4) 
manufactured_solution_dwdx(x, y) = 4 * C * x^3 
manufactured_solution_dwdy(x, y) = 4 * C * y^3

function _execute(n = 2, t_radius_ratio = 0.01, mesh_procedure = :t3_nice, visualize = true)
    E = 200*phun("GPa");
    nu = 0.3;
    a = 1.0*phun("m");
    thickness = a * t_radius_ratio;
    D = E*thickness^3/12/(1-nu^2)
    q = 48 * C * D
    center_w = 0;
    formul = FEMMShellT3FFModule

    tolerance = a/n/1000
    if mesh_procedure == :q4_t3
        fens, fes = Q4circlen(a, n);
        fens, fes = Q4toT3(fens, fes)
    elseif mesh_procedure == :t3_nice
        fens, fes = T3circlen(a, n);
    elseif mesh_procedure == :t3_poor
        fens, fes = T3circleseg(pi/2, a, n, n);
    elseif mesh_procedure == :t3_targe2
        commands = """
   curve 1 circle center 0.0 0.0 radius $(a)
   subregion 1  property 1 boundary 1
   m-ctl-point constant $(a/n)
   """
        mesh = triangulate(commands)
        fens = FENodeSet(mesh.xy)
        fes = FESetT3(mesh.triconn)
    else
        @error "Unknown mesh procedure"
    end

    fens.xyz = xyz3(fens)

    renumb = (c) -> c[[1, 3, 2]] # make sure the mirrored copies have consistent orientation
    fens1, fes1 = mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], renumb = renumb)
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, tolerance)
    fes = cat(fesa[1], fesa[2])
    fens1, fes1 = mirrormesh(fens, fes, [0.0, -1.0, 0.0], [0.0, 0.0, 0.0], renumb = renumb)
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, tolerance)
    fes = cat(fesa[1], fesa[2])
    
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
    l1 = connectednodes(bfes)
    for i in [1,2,6]
        setebc!(dchi, l1, true, i)
    end
    for j in l1
        setebc!(dchi, j, true, 3, manufactured_solution_w(fens.xyz[j, 1:2]...))
        setebc!(dchi, j, true, 4, manufactured_solution_dwdy(fens.xyz[j, 1:2]...))
        setebc!(dchi, j, true, 5, -manufactured_solution_dwdx(fens.xyz[j, 1:2]...))
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], tolerance = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64[0, 0, q, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    
    # targetu =  dchi.values[nl, 3][1]
    # @info "Target: $(round(targetu, digits=8)),  $(round(targetu/center_w, digits = 4)*100)%"

    # Generate a graphical display of resultants
    cylindrical!(csmatout, XYZ, tangents, feid, qpid) = begin
        csmatout[:, 1] .= vec(XYZ)/norm(vec(XYZ))
        csmatout[:, 3] .= (0.0, 0.0, 1.0)
        cross3!(view(csmatout, :, 2), view(csmatout, :, 3), view(csmatout, :, 1))
        return csmatout
    end
    ocsys = CSys(3)
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
        push!(scalars, ("m$nc", fld.values))
    end
    vtkwrite("manufactured_circular_plate_examples-m.vtu", fens, fes; scalars = scalars)
    # Generate a graphical display of resultants
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
        push!(scalars, ("q$nc", fld.values))
    end
    vtkwrite("manufactured_circular_plate_examples-q.vtu", fens, fes; scalars = scalars)
    vtkwrite("manufactured_circular_plate_examples-uur.vtu", fens, fes; 
        vectors = [
            ("u", deepcopy(dchi.values[:, 1:3])),
            ("ur", deepcopy(dchi.values[:, 4:6])),
        ])
    # Visualization
    if !visualize
        return true
    end
    
    return true
end

function test_convergence()
    formul = FEMMShellT3FFModule
    t_radius_ratio = 0.01
    @info "Simply supported square plate with uniform load,"
    @info "thickness/length = $t_radius_ratio "
    for n in [18, ]
        _execute(n, t_radius_ratio, :t3_targe2, false)
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
