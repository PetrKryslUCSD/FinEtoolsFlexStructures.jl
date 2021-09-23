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
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3DSGAModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

using Infiltrator

function _execute_dsg_model(formul, n = 8, visualize = true)
    E = 6.825e7;
    nu = 0.3;
    thickness  =  0.04;
    # analytical solution for the vertical deflection under the load
    analyt_sol = 0.093;
    R = 10.0;

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

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 R R 0 0], inflate = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, -1, 0, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    nl = selectnode(fens; box = Float64[R R 0 0 0 0], inflate = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[1, 0, 0, 0, 0, 0]);
    F += distribloads(lfemm, geom0, dchi, fi, 3);


    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    resultpercent =  dchi.values[nl, 1][1]*100
    @info "Solution: $(round(resultpercent/analyt_sol, digits = 4))%"

    # Generate a graphical display of resultants
    spherical!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) = begin
        r = vec(XYZ); 
        csmatout[:, 3] .= vec(r)/norm(vec(r))
        csmatout[:, 2] .= (0.0, 0.0, 1.0)
        cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
        cross3!(view(csmatout, :, 2), view(csmatout, :, 3), view(csmatout, :, 1))
        return csmatout
    end
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

function test_convergence(formul)
    @info "Hemisphere, formulation=$(formul)"
    for n in [2, 4, 8, 16, 32, 64, 128, 256]
        _execute_dsg_model(formul, n, false)
    end
    return true
end

end # module

using .hemisphere_open_examples
using FinEtoolsFlexStructures.FEMMShellT3DSGAModule
hemisphere_open_examples.test_convergence(FEMMShellT3DSGAModule)

