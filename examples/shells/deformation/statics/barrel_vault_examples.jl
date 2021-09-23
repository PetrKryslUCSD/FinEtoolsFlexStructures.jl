"""

"""
module barrel_vault_examples

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3DSGAModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

using Infiltrator

function _execute_dsg_model(formul, input = "barrelvault_s3r_fineirreg.inp", visualize = true)
    E = 3.0e6;
    nu = 0.0;
    thickness = 3.0;
    tolerance = thickness/20
    # analytical solution for the vertical deflection under the load
    analyt_sol = -0.3024*12;
    R = 25.0*12;
    L = 50.0*12;
    
    output = import_ABAQUS(joinpath(dirname(@__FILE__()), input))
    fens = output["fens"]
    fes = output["fesets"][1]

     # boundingbox(fens.xyz)
    connected = findunconnnodes(fens, fes);
    fens, new_numbering = compactnodes(fens, connected);
    fes = renumberconn!(fes, new_numbering);

    fens, fes = mergenodes(fens, fes, thickness/10)

    # plots = cat(plot_space_box([[0 0 -R/2]; [R/2 R/2 R/2]]),
    #     plot_nodes(fens),
    #     plot_midsurface(fens, fes);
    # dims = 1)
    # pl = render(plots)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    @info "Mesh: $input"

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

    # vtkwrite("barrel_vault_examples-mesh.vtu", fens, fes)
    # bfes = meshboundary(fes)
    # vtkwrite("barrel_vault_examples-boundary-mesh.vtu", fens, bfes)
    
    # Apply EBC's
    # Rigid diaphragm end
    l1 = selectnode(fens; box = Float64[-Inf Inf -Inf Inf 0 0], inflate = tolerance)
    for i in [1,2]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Z
    l1 = selectnode(fens; box = Float64[-Inf Inf -Inf Inf L/2 L/2], inflate = tolerance)
    for i in [3,4,5]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y (apex)
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[-Inf Inf sin(40/360*2*pi)*R sin(40/360*2*pi)*R L/2 L/2], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt[-0.625, 0, 0, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])

    targetu =  dchi.values[nl, 1][1]
    @info "Solution: $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

    # Generate a graphical display of resultants
    cylindrical!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) = begin
        r = -vec(XYZ); r[3] = 0.0
        csmatout[:, 3] .= vec(r)/norm(vec(r))
        csmatout[:, 2] .= (0.0, 0.0, 1.0)
        cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
        return csmatout
    end
    ocsys = CSys(3, 3, cylindrical!)
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("m$nc", fld.values))
    end
    vtkwrite("$(input)-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("n$nc", fld.values))
    end
    vtkwrite("$(input)-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("q$nc", fld.values))
    end
    vtkwrite("$(input)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (R/2)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -R]; [R R R]]),
        plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_convergence(formul)
    
    @info "Scordelis-Lo Abaqus model, formulation=$(formul)"

    _execute_dsg_model(formul, "barrelvault_stri3_irreg.inp", false)

    _execute_dsg_model(formul, "barrelvault_s3r_fineirreg.inp", false)

    return true
end

end # module

using .barrel_vault_examples
using FinEtoolsFlexStructures.FEMMShellT3DSGAModule
barrel_vault_examples.test_convergence(FEMMShellT3DSGAModule)