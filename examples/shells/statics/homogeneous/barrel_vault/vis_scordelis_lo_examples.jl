"""
The barrel vault (Scordelis-Lo) roof is one of the benchmarks for linear elastic
analysis of shells. 

The candidate element's usefulness in irregular geometries (and most practical
cases involve a high degree of geometric irregularity) is tested. As would be
expected,the irregular mesh results are not as good as those provided by a
regular mesh with the same number of variables. 

Problem description

The physical basis of the problem is a deeply arched roof supported only
by diaphragms at its curved edges (an aircraft hanger), deforming under its own
weight. It is interesting to observe that the geometry is such that the
centerpoint of the roof moves upward under the self-weight(downwardly directed)
load. Perhaps this is one reason why the problem is not straightforward
numerically.

Analytical solution for the vertical deflection and the midpoint of the
free edge is often cited as 0.3024. However, this number is suspect. 
"""
module vis_scordelis_lo_examples

using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.FEMMShellQ4RNTModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

# Analytical solution?
analyt_sol=-0.3024;
# Parameters:
E=4.32e8;
nu=0.0;
thickness = 0.25; # geometrical dimensions are in feet
R = 25.0;
L = 50.0;

cylindrical!(csmatout, XYZ, tangents, feid, qpid) = begin
    r = vec(XYZ); r[2] = 0.0; r[3] += R
    csmatout[:, 3] .= vec(r)/norm(vec(r))
    csmatout[:, 2] .= (0.0, 1.0, 0.0) #  this is along the axis
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

function _execute_t3ff(formul, n = 8, visualize = true)
    tolerance = R/n/10
    # fens, fes = T3blockrand(40/360*2*pi,L/2,n,n);
    fens, fes = T3block(40/360*2*pi,L/2,n,n,:b);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*(cos(a)-1))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    ocsys = CSys(3, 3, cylindrical!)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness),  mater)
    # femm.drilling_stiffness_scale = 1.0e-4
    # femm.mult_el_size = 5/12
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # rigid diaphragm
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,3,5]
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
    applyebc!(dchi)
    numberdofs!(dchi);

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    Kff = stiffness(femm, massem, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[sin(40/360*2*pi)*25 sin(40/360*2*pi)*25 L/2 L/2 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64[0, 0, -90, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 3);
    
    # Solve
     Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)
    
    result =   dchi.values[nl, 3][1]
    @info "Solution: $(result), $(round(result/analyt_sol*100, digits = 4))%"

    # Visualization
    if visualize
        # Generate a graphical display of resultants
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            push!(scalars, ("m$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            push!(scalars, ("em$nc", fld.values))
        end
        vtkwrite("vis_scordelis_lo_t3ff-$(n)-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            push!(scalars, ("n$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            push!(scalars, ("en$nc", fld.values))
        end
        vtkwrite("vis_scordelis_lo_t3ff-$(n)-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("q$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("eq$nc", fld.values))
        end
        vtkwrite("vis_scordelis_lo_t3ff-o-$(n)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

        # vtkwrite("vis_scordelis_lo_t3ff-$(n)-uur.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])

        # scattersysvec!(dchi, (L/8)/maximum(abs.(U)).*U)
        # update_rotation_field!(Rfield0, dchi)
        # plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #     #plot_nodes(fens),
        #     plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
        #     plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
        #     dims = 1)
        # pl = render(plots)
    end

    result
end

function _execute_q4rnt(formul, n = 8, visualize = true)
    tolerance = R/n/10
    # fens, fes = T3blockrand(40/360*2*pi,L/2,n,n);
    fens, fes = Q4block(40/360*2*pi,L/2,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*(cos(a)-1))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    ocsys = CSys(3, 3, cylindrical!)
    
    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, GaussRule(2, 2), thickness),  mater)
    # femm.drilling_stiffness_scale = 1.0e-4
    # femm.mult_el_size = 5/12
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # rigid diaphragm
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,3,5]
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
    applyebc!(dchi)
    numberdofs!(dchi);

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    Kff = stiffness(femm, massem, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[sin(40/360*2*pi)*25 sin(40/360*2*pi)*25 L/2 L/2 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(Float64[0, 0, -90, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 2);
    
    # Solve
     Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)
    
    result =   dchi.values[nl, 3][1]
    @info "Solution: $(result), $(round(result/analyt_sol*100, digits = 4))%"

    # Visualization
    if visualize
        # Generate a graphical display of resultants
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            push!(scalars, ("m$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            push!(scalars, ("em$nc", fld.values))
        end
        vtkwrite("vis_scordelis_lo_q4rnt-$(n)-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            push!(scalars, ("n$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            push!(scalars, ("en$nc", fld.values))
        end
        vtkwrite("vis_scordelis_lo_q4rnt-$(n)-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("q$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("eq$nc", fld.values))
        end
        vtkwrite("vis_scordelis_lo_q4rnt-o-$(n)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

        # vtkwrite("vis_scordelis_lo_q4rnt-$(n)-uur.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])

        # scattersysvec!(dchi, (L/8)/maximum(abs.(U)).*U)
        # update_rotation_field!(Rfield0, dchi)
        # plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #     #plot_nodes(fens),
        #     plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
        #     plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
        #     dims = 1)
        # pl = render(plots)
    end

    result
end

function test_t3ff(ns = [16, 32], visualize = true)
    formul = FEMMShellT3FFModule
    @info "Scordelis-Lo shell, formulation=$(formul)"
    results = []
    for n in ns
        v = _execute_t3ff(formul, n, visualize)
        push!(results, v/(-0.3020)*100)
    end
    return ns, results
end

function test_q4rnt(ns = [16, 32], visualize = true)
    formul = FEMMShellQ4RNTModule
    @info "Scordelis-Lo shell, formulation=$(formul)"
    results = []
    for n in ns
        v = _execute_q4rnt(formul, n, visualize)
        push!(results, v/(-0.3020)*100)
    end
    return ns, results
end

function allrun()
    println("#####################################################")
    println("# test_t3ff ")
    test_t3ff()
    println("#####################################################")
    println("# test_q4rnt ")
    test_q4rnt()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
