"""
MODEL DESCRIPTION

Z-section cantilever under torsional loading.

Linear elastic analysis, Young's modulus = 210 GPa, Poisson's ratio = 0.3.

All displacements are fixed at X=0.

Torque of 1.2 MN-m applied at X=10. The torque is applied by two 
uniformly distributed shear loads of 0.6 MN at each flange surface.

Objective of the analysis is to compute the axial stress at X = 2.5 from fixed end.

NAFEMS REFERENCE SOLUTION

Axial stress at X = 2.5 from fixed end (point A) at the midsurface is -108 MPa.
"""
module LE5_Z_cantilever_examples

using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

function zcant!(csmatout, XYZ, tangents, feid, qpid)
    r = vec(XYZ); 
    cross3!(r, view(tangents, :, 1), view(tangents, :, 2))
    csmatout[:, 3] .= vec(r)/norm(vec(r))
    csmatout[:, 1] .= (1.0, 0.0, 0.0)
    cross3!(view(csmatout, :, 2), view(csmatout, :, 3), view(csmatout, :, 1))
    return csmatout
end


function _execute_model(formul, input = "nle5xf3c.inp", nrefs = 0, visualize = true)
    E = 210e9;
    nu = 0.3;
    L = 10.0;
    thickness = 0.1

    tolerance = thickness/1000
    output = import_ABAQUS(joinpath(dirname(@__FILE__()), input))
    fens = output["fens"]
    fes = output["fesets"][1]

    connected = findunconnnodes(fens, fes);
    fens, new_numbering = compactnodes(fens, connected);
    fes = renumberconn!(fes, new_numbering);

    for r in 1:nrefs
        fens, fes = T3refine(fens, fes)
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    ocsys = CSys(3, 3, zcant!)
    
    # Report
    @info "Mesh: $input, nrefs = $nrefs"

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


    # plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
    #     plot_nodes(fens),
    #     plot_midsurface(fens, fes);
    # dims = 1)
    # pl = render(plots)

    # Apply EBC's
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,]
        setebc!(dchi, l1, true, i)
    end
    
    applyebc!(dchi)
    numberdofs!(dchi);
    @show nfreedofs(dchi)

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, massem, geom0, u0, Rfield0, dchi);
    
    # Load
    nl = selectnode(fens; box = Float64[10.0 10.0 1.0 1.0 0 0], tolerance = tolerance)
    loadbdry1 = FESetP1(reshape(nl, 1, 1))
    lfemm1 = FEMMBase(IntegDomain(loadbdry1, PointRule()))
    fi1 = ForceIntensity(Float64[0, 0, +0.6e6, 0, 0, 0]);
    nl = selectnode(fens; box = Float64[10.0 10.0 -1.0 -1.0 0 0], tolerance = tolerance)
    loadbdry2 = FESetP1(reshape(nl, 1, 1))
    lfemm2 = FEMMBase(IntegDomain(loadbdry2, PointRule()))
    fi2 = ForceIntensity(Float64[0, 0, -0.6e6, 0, 0, 0]);
    F = distribloads(lfemm1, vassem, geom0, dchi, fi1, 3) + distribloads(lfemm2, vassem, geom0, dchi, fi2, 3);

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    targetu =  minimum(dchi.values[:, 3]), maximum(dchi.values[:, 3])
    @info "Target: $(round.(targetu, digits=8))"

    # Generate a graphical display of displacements and rotations
    scalars = []
    for nc in 1:6
        push!(scalars, ("dchi$nc", deepcopy(dchi.values[:, nc])))
    end
    vectors = []
    push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
    push!(vectors, ("UR", deepcopy(dchi.values[:, 4:6])))
    vtkwrite("z_cant-$input-$nrefs-dchi.vtu", fens, fes; scalars = scalars, vectors = vectors)

    # Generate a graphical display of resultants
    nl = selectnode(fens; nearestto = Float64[2.5 1.0 1.0])
    
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("m$nc", fld.values))
    end
    vtkwrite("z_cant-$input-$nrefs-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    pointAstresses = Float64[]
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(pointAstresses, fld.values[nl[1], 1]/thickness)
        push!(scalars, ("n$nc", fld.values))
    end
    vtkwrite("z_cant-$input-$nrefs-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("q$nc", fld.values))
    end
    vtkwrite("z_cant-$input-$nrefs-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

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
    return pointAstresses
end

function _execute_model_alt(formul, ignore = "", nrefs = 0, visualize = true)
    E = 210e9;
    nu = 0.3;
    L = 10.0;
    thickness = 0.1

    tolerance = thickness/100
    xyz = Float64[
    0 -1 -1;
    2.5 -1 -1;
    10 -1 -1;
    0 -1 0;
    2.5 -1 0;
    10 -1 0;
    0 1 0;
    2.5 1 0;
    10 1 0;
    0 1 1;
    2.5 1 1;
    10 1 1;
    ]
    fens1, fes1 = Q4quadrilateral(xyz[[1, 2, 5, 4], :], 2, 1) 
    fens2, fes2 = Q4quadrilateral(xyz[[4, 5, 8, 7], :], 2, 1) 
    fens, fes1, fes2 = mergemeshes(fens1, fes1,  fens2, fes2, tolerance)
    fes = cat(fes1, fes2)

    fens1, fes1 = fens, fes 
    fens2, fes2 = Q4quadrilateral(xyz[[7, 8, 11, 10], :], 2, 1) 
    fens, fes1, fes2 = mergemeshes(fens1, fes1,  fens2, fes2, tolerance)
    fes = cat(fes1, fes2)

    fens1, fes1 = fens, fes 
    fens2, fes2 = Q4quadrilateral(xyz[[2, 3, 6, 5], :], 6, 1) 
    fens, fes1, fes2 = mergemeshes(fens1, fes1,  fens2, fes2, tolerance)
    fes = cat(fes1, fes2)

    fens1, fes1 = fens, fes 
    fens2, fes2 = Q4quadrilateral(xyz[[5, 6, 9, 8], :], 6, 1) 
    fens, fes1, fes2 = mergemeshes(fens1, fes1,  fens2, fes2, tolerance)
    fes = cat(fes1, fes2)

    fens1, fes1 = fens, fes 
    fens2, fes2 = Q4quadrilateral(xyz[[8, 9, 12, 11], :], 6, 1) 
    fens, fes1, fes2 = mergemeshes(fens1, fes1,  fens2, fes2, tolerance)
    fes = cat(fes1, fes2)

    fens, fes = Q4toT3(fens, fes)

    for r in 1:nrefs
        fens, fes = T3refine(fens, fes)
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    ocsys = CSys(3, 3, zcant!)
    
    # Report
    @info "Mesh: nrefs = $nrefs"

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

    # plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
    #     plot_nodes(fens),
    #     plot_midsurface(fens, fes);
    # dims = 1)
    # pl = render(plots)

    # Apply EBC's
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,]
        setebc!(dchi, l1, true, i)
    end
    
    applyebc!(dchi)
    numberdofs!(dchi);
    @show nfreedofs(dchi)

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, massem, geom0, u0, Rfield0, dchi);
    
    # Load
    # bfes = meshboundary(fes)
    # el = selectelem(fens, bfes, box = Float64[10.0 10.0 1.0 1.0 -1 1], tolerance = tolerance)
    # lfemm1 = FEMMBase(IntegDomain(subset(bfes, el), GaussRule(1, 2)))
    # fi1 = ForceIntensity(Float64[0, 0, +0.6e6, 0, 0, 0]);
    # el = selectelem(fens, bfes, box = Float64[10.0 10.0 -1.0 -1.0 -1 1], tolerance = tolerance)
    # lfemm2 = FEMMBase(IntegDomain(subset(bfes, el), GaussRule(1, 2)))
    # fi2 = ForceIntensity(Float64[0, 0, -0.6e6, 0, 0, 0]);
    # F = distribloads(lfemm1, vassem, geom0, dchi, fi1, 1) + distribloads(lfemm2, vassem, geom0, dchi, fi2, 1);

    nl = selectnode(fens; box = Float64[10.0 10.0 1.0 1.0 0 0], tolerance = tolerance)
    loadbdry1 = FESetP1(reshape(nl, 1, 1))
    lfemm1 = FEMMBase(IntegDomain(loadbdry1, PointRule()))
    fi1 = ForceIntensity(Float64[0, 0, +0.6e6, 0, 0, 0]);
    nl = selectnode(fens; box = Float64[10.0 10.0 -1.0 -1.0 0 0], tolerance = tolerance)
    loadbdry2 = FESetP1(reshape(nl, 1, 1))
    lfemm2 = FEMMBase(IntegDomain(loadbdry2, PointRule()))
    fi2 = ForceIntensity(Float64[0, 0, -0.6e6, 0, 0, 0]);
    F = distribloads(lfemm1, vassem, geom0, dchi, fi1, 3) + distribloads(lfemm2, vassem, geom0, dchi, fi2, 3);

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    targetu =  minimum(dchi.values[:, 3]), maximum(dchi.values[:, 3])
    @info "Target: $(round.(targetu, digits=8))"

    # Generate a graphical display of displacements and rotations
    scalars = []
    for nc in 1:6
        push!(scalars, ("dchi$nc", deepcopy(dchi.values[:, nc])))
    end
    vectors = []
    push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
    push!(vectors, ("UR", deepcopy(dchi.values[:, 4:6])))
    vtkwrite("z_cant-alt-$nrefs-dchi.vtu", fens, fes; scalars = scalars, vectors = vectors)

    # Generate a graphical display of resultants
    nl = selectnode(fens; nearestto = Float64[2.5 1.0 1.0])
    
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("m$nc", fld.values))
    end
    vtkwrite("z_cant-alt-$nrefs-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    pointAstresses = Float64[]
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(pointAstresses, fld.values[nl[1], 1]/thickness)
        push!(scalars, ("n$nc", fld.values))
    end
    vtkwrite("z_cant-alt-$nrefs-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("q$nc", fld.values))
    end
    vtkwrite("z_cant-alt-$nrefs-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

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
    return pointAstresses
end

function test_convergence()
    formul = FEMMShellT3FFModule
    @info "LE5 Z-cantilever, formulation=$(formul)"
    for n in [0, 1, 2, 3, 4, ]
        _execute_model(formul, "nle5xf3c.inp", n, false)
    end
    return true
end

function test_convergence_alt()
    formul = FEMMShellT3FFModule
    @info "LE5 Z-cantilever, formulation=$(formul)"
    for n in [0, 1, 2, 3, 4, 5]
        _execute_model_alt(formul, "ignore", n, false)
    end
    return true
end


function allrun()
    println("#####################################################")
    println("# test_convergence ")
    test_convergence()
    println("#####################################################")
    println("# test_convergence_alt ")
    test_convergence_alt()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
