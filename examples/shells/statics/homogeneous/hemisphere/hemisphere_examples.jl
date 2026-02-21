"""
The pinched hemisphere benchmark. Half over fully closed sphere.

"""
module hemisphere_examples

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.FEMMShellQ4RNTModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

function _execute_t3ff(n = 2, visualize = true)
    E = 6.825e7;
    nu = 0.3;
    thickness  =  0.04;
    # analytical solution for the vertical deflection under the load
    analyt_sol = 0.0924;
    R = 10.0;
    formul = FEMMShellT3FFModule

    tolerance = R/n/1000
    fens, fes = Q4spheren(R, n)
    fens, fes = Q4toT3(fens, fes)

    spherical!(csmatout, XYZ, tangents, feid, qpid) = begin
        r = vec(XYZ); 
        csmatout[:, 3] .= vec(r)/norm(vec(r))
        csmatout[:, 2] .= (0.0, 0.0, 1.0)
        cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
        csmatout[:, 1] .= vec(view(csmatout, :, 1))/norm(vec(view(csmatout, :, 1)))
        cross3!(view(csmatout, :, 2), view(csmatout, :, 3), view(csmatout, :, 1))
        return csmatout
    end
    ocsys = CSys(3, 3, spherical!)
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    # Report
    @info "Mesh: $n elements per side"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    femm.drilling_stiffness_scale = 0.1
    femm.threshold_angle = 45.0
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
    # top
    l1 = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], inflate = tolerance)
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
    nl = selectnode(fens; box = Float64[0 0 R R 0 0], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(Float64[0, -1, 0, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 3);
    nl = selectnode(fens; box = Float64[R R 0 0 0 0], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(Float64[1, 0, 0, 0, 0, 0]);
    Ff += distribloads(lfemm, vassem, geom0, dchi, fi, 3);


    # @infiltrate
    # Solve
     Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)
    
    resultpercent =  dchi.values[nl, 1][1]*100
    @info "Solution: $(round(resultpercent/analyt_sol, digits = 4))%"

    # formul._resultant_check(femm, geom0, u0, Rfield0, dchi)

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

function _execute_q4rnt(n = 2, visualize = true)
    E = 6.825e7;
    nu = 0.3;
    thickness  =  0.04;
    # analytical solution for the vertical deflection under the load
    analyt_sol = 0.0924;
    R = 10.0;
    formul = FEMMShellQ4RNTModule

    tolerance = R/n/1000
    fens, fes = Q4spheren(R, n)

    spherical!(csmatout, XYZ, tangents, feid, qpid) = begin
        r = vec(XYZ); 
        csmatout[:, 3] .= vec(r)/norm(vec(r))
        csmatout[:, 2] .= (0.0, 0.0, 1.0)
        cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
        csmatout[:, 1] .= vec(view(csmatout, :, 1))/norm(vec(view(csmatout, :, 1)))
        cross3!(view(csmatout, :, 2), view(csmatout, :, 3), view(csmatout, :, 1))
        return csmatout
    end
    ocsys = CSys(3, 3, spherical!)
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    # Report
    @info "Mesh: $n elements per side"

    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.make(
        IntegDomain(fes, GaussRule(2, 2), 
        thickness), ocsys, mater)
    femm.drilling_stiffness_scale = 0.1
    femm.threshold_angle = 45.0
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
    # top
    l1 = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], inflate = tolerance)
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
    nl = selectnode(fens; box = Float64[0 0 R R 0 0], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(Float64[0, -1, 0, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 3);
    nl = selectnode(fens; box = Float64[R R 0 0 0 0], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(Float64[1, 0, 0, 0, 0, 0]);
    Ff += distribloads(lfemm, vassem, geom0, dchi, fi, 3);


    # @infiltrate
    # Solve
     Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)
    
    resultpercent =  dchi.values[nl, 1][1]*100
    @info "Solution: $(round(resultpercent/analyt_sol, digits = 4))%"

    # formul._resultant_check(femm, geom0, u0, Rfield0, dchi)

    # Visualization
    if visualize
        vtkwrite("hemisphere-$(n)-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        return true
    end
    
    return true
end

function test_convergence()
    @info "Hemisphere benchmark, T3FF elements"
    # for n in [2, ]
    for n in [2, 4, 8, 16, 32]
        _execute_t3ff(n, false)
    end
    @info "Hemisphere benchmark, Q4RNT elements"
    # for n in [2, ]
    for n in [2, 4, 8, 16, 32]
        _execute_q4rnt(n, true)
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

