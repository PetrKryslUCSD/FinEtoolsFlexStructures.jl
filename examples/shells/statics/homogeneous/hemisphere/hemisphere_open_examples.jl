"""
The pinched hemisphere benchmark for the configuration with a 18 deg hole at the
top 

The reference below states:

The spherical shell shown in Fig. 9 is our proposed doubly-curved shell
problem. Note that the equator is a free edge so that the problem
represents a hemisphere with four point loads alternating in sign at 90 °
intervals on the equator. The hole at the top has been introduced to
avoid the use of triangles near the axis of revolution. Convergence can
be studied by varying  mesh size. Both membrane and bending strains
contribute significantly to the radial displacement at the load point. A
theoretical value of the displacement under load has been computed for a
slightly different configuration [7] in which the hole at the axis is
closed.

Macneal RH, Harder RL (1985) A proposed standard set of problems to test
finite element accuracy. Finite Elements in Analysis and Design 1: 3-20.

Regarding the convergence: refer to the article
Performance of the MITC3+ and MITC4+ shell elements in widely-used
benchmark problems
Yeongbin Ko a, Youngyu Lee b, Phill-Seung Lee a,⇑, Klaus-Jürgen Bathe

The drilling degrees of freedom stiffness has an effect on the convergence for
coarser meshes.
"""
module hemisphere_open_examples

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

function spherical!(csmatout, XYZ, tangents, feid, qpid)
    r = vec(XYZ); 
    csmatout[:, 3] .= vec(r)/norm(vec(r))
    csmatout[:, 2] .= (0.0, 0.0, 1.0)
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    csmatout[:, 1] .= vec(view(csmatout, :, 1))/norm(vec(view(csmatout, :, 1)))
    cross3!(view(csmatout, :, 2), view(csmatout, :, 3), view(csmatout, :, 1))
    return csmatout
end

function _execute_w_approx_normals(n = 8, visualize = true, drilling_stiffness_multiplier = 1.0)
    E = 6.825e7;
    nu = 0.3;
    thickness  =  0.04;
    # analytical solution for the vertical deflection under the load
    analyt_sol = 0.093;
    R = 10.0;
    formul = FEMMShellT3FFModule

    tolerance = R/n/100
    fens, fes = Q4block(90.0, 70.0, n, n)
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        phi, psi = fens.xyz[i, 1:2]
        fens.xyz[i, :] .= (cos(psi/180*pi) .* cos(phi/180*pi)*R, cos(psi/180*pi) .* sin(phi/180*pi)*R, sin(psi/180*pi) .* R)
    end
    fens, fes = Q4toT3(fens, fes)
    
    vtkwrite("geom.vtu", fens, fes)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    # Report
    @info "Mesh: $n elements per side"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    
    femm.drilling_stiffness_scale = 0.1 * drilling_stiffness_multiplier
    femm.mult_el_size = 5/12/1.5 
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!
    num_normals = formul.num_normals

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    # one point from the list
    for i in [3]
        setebc!(dchi, l1[1:1], true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    total_normals, invalid_normals = num_normals(femm)
    @assert invalid_normals == 0
    K = stiffness(femm, massem, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 R R 0 0], inflate = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(Float64[0, -1, 0, 0, 0, 0]);
    F = distribloads(lfemm, vassem, geom0, dchi, fi, 3);
    nl = selectnode(fens; box = Float64[R R 0 0 0 0], inflate = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(Float64[1, 0, 0, 0, 0, 0]);
    F += distribloads(lfemm, vassem, geom0, dchi, fi, 3);


    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    resultpercent =  dchi.values[nl, 1][1]/analyt_sol*100
    @info "Solution: $(round(resultpercent, digits = 4))% ($analyt_sol)"

    # formul._resultant_check(femm, geom0, u0, Rfield0, dchi)

    # Generate a graphical display of displacements and rotations
    scalars = []
    for nc in 1:6
        push!(scalars, ("dchi$nc", deepcopy(dchi.values[:, nc])))
    end
    vectors = []
    push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
    push!(vectors, ("UR", deepcopy(dchi.values[:, 4:6])))
    vtkwrite("hemisphere_open-$n-dchi.vtu", fens, fes; scalars = scalars, vectors = vectors)

    # Generate a graphical display of resultants
    ocsys = CSys(3, 3, spherical!)
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("m$nc", fld.values))
    end
    vtkwrite("hemisphere_open-$n-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("n$nc", fld.values))
    end
    vtkwrite("hemisphere_open-$n-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("q$nc", fld.values))
    end
    vtkwrite("hemisphere_open-$n-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (R/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -R]; [R R R]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function _execute_w_exact_normals(n = 8, visualize = true, drilling_stiffness_multiplier = 1.0)
    E = 6.825e7;
    nu = 0.3;
    thickness  =  0.04;
    # analytical solution for the vertical deflection under the load
    analyt_sol = 0.093;
    R = 10.0;
    formul = FEMMShellT3FFModule

    tolerance = R/n/100
    fens, fes = Q4block(90.0, 70.0, n, n)
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        phi, psi = fens.xyz[i, 1:2]
        fens.xyz[i, :] .= (cos(psi/180*pi) .* cos(phi/180*pi)*R, cos(psi/180*pi) .* sin(phi/180*pi)*R, sin(psi/180*pi) .* R)
    end
    fens, fes = Q4toT3(fens, fes)
    
    vtkwrite("geom.vtu", fens, fes)
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    # Report
    @info "Mesh: $n elements per side"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    ocsys = CSys(3, 3, spherical!)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    femm.drilling_stiffness_scale = 0.1 * drilling_stiffness_multiplier
    femm.mult_el_size = 5/12/1.5  
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
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    # one point from the list
    for i in [3]
        setebc!(dchi, l1[1:1], true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, massem, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 R R 0 0], inflate = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(Float64[0, -1, 0, 0, 0, 0]);
    F = distribloads(lfemm, vassem, geom0, dchi, fi, 3);
    nl = selectnode(fens; box = Float64[R R 0 0 0 0], inflate = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(Float64[1, 0, 0, 0, 0, 0]);
    F += distribloads(lfemm, vassem, geom0, dchi, fi, 3);


    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    resultpercent =  dchi.values[nl, 1][1]/analyt_sol*100
        @info "Solution: $(round(resultpercent, digits = 4))% ($analyt_sol)"

    # formul._resultant_check(femm, geom0, u0, Rfield0, dchi)

    # Generate a graphical display of displacements and rotations
    scalars = []
    for nc in 1:6
        push!(scalars, ("dchi$nc", deepcopy(dchi.values[:, nc])))
    end
    vectors = []
    push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
    push!(vectors, ("UR", deepcopy(dchi.values[:, 4:6])))
    vtkwrite("hemisphere_open-$n-dchi.vtu", fens, fes; scalars = scalars, vectors = vectors)

    # Generate a graphical display of resultants
    
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("m$nc", fld.values))
    end
    vtkwrite("hemisphere_open-$n-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("n$nc", fld.values))
    end
    vtkwrite("hemisphere_open-$n-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("q$nc", fld.values))
    end
    vtkwrite("hemisphere_open-$n-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (R/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -R]; [R R R]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_convergence()
    drilling_stiffness_multiplier = 10000.0
    drilling_stiffness_multiplier = 0.0001
    @info "Hemisphere with opening: with exact normals"
    for n in [2, 4, 8, 16, 32, 64, 128, ]
        _execute_w_exact_normals(n, false, drilling_stiffness_multiplier)
    end
    @info "Hemisphere with opening: with approximate normals"
    for n in [2, 4, 8, 16, 32, 64, 128, ]
        _execute_w_approx_normals(n, false, drilling_stiffness_multiplier)
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

