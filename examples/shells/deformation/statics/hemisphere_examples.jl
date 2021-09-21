"""
The pinched hemisphere benchmark
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

In our example the shell is closed at the top, and we avoid the use
of triangles on the axis of symmetry.

Macneal RH, Harder RL (1985) A proposed standard set of problems to test
finite element accuracy. Finite Elements in Analysis and Design 1: 3-20.
"""
module hemisphere_examples

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3DSGAModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

using Infiltrator

function _execute_dsg_model(formul, n = 2, visualize = true)
    E = 6.825e7;
    nu = 0.3;
    thickness  =  0.04;
    # analytical solution for the vertical deflection under the load
    analyt_sol = 0.0924;
    R = 10.0;

    tolerance = R/n/1000
    fens, fes = Q4spheren(R, n)
    fens, fes = Q4toT3(fens, fes)

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
    # top
    l1 = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 R R 0 0], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, -1, 0, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    nl = selectnode(fens; box = Float64[R R 0 0 0 0], tolerance = tolerance)
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
    for n in [2, 4, 8, 16, 32, 64]
        _execute_dsg_model(formul, n, false)
    end
    return true
end

end # module

using .hemisphere_examples
using FinEtoolsFlexStructures.FEMMShellT3DSGAModule
hemisphere_examples.test_convergence(FEMMShellT3DSGAModule)

