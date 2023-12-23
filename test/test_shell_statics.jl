
module scordelis_lo_dsg3_verification

using Test
using LinearAlgebra
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!

function _execute(n = 8, visualize = true)
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    analyt_sol = -0.3024
    # Parameters:
    E = 4.32e8
    nu = 0.0
    thickness = 0.25 # geometrical dimensions are in feet
    R = 25.0
    L = 50.0
    formul = FEMMShellT3FFModule

    tolerance = R / n / 1000
    fens, fes = T3block(40 / 360 * 2 * pi, L / 2, n, n)
    fens.xyz = xyz3(fens)
    for i = 1:count(fens)
        a = fens.xyz[i, 1]
        y = fens.xyz[i, 2]
        fens.xyz[i, :] .= (R * sin(a), y, R * (cos(a) - 1))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    femm.mult_el_size = 0.2
    femm.drilling_stiffness_scale = 1.0
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # rigid diaphragm
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1, 3, 5]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf L / 2 L / 2 -Inf Inf], inflate = tolerance)
    for i in [2, 4, 6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1, 5, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Midpoint of the free edge
    nl = selectnode(
        fens;
        box = Float64[sin(40 / 360 * 2 * pi) * 25 sin(40 / 360 * 2 * pi) * 25 L / 2 L / 2 -Inf Inf],
        inflate = tolerance,
    )
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt[0, 0, -90, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 3)

    # Solve
    solve_blocked!(dchi, K, F)
    resultpercent = dchi.values[nl, 3][1] / analyt_sol * 100

    return resultpercent
end

function test_convergence()
    formul = FEMMShellT3FFModule
    results = [
        66.54771615057949,
        85.54615143134853,
        89.85075281481419,
        92.50616661644985,
        95.40469210310079,
        97.65376880486126,
    ]
    for (n, res) in zip([4, 8, 10, 12, 16, 24], results)
        v = _execute(n, false)
        @test isapprox(res, v, rtol = 1.0e-4)
        # @show v
    end
    return true
end

test_convergence()

end # module

"""
The geometry consists of a “hook” in the form of a curved strip rigidly clamped
at one end and loaded with a unit in-plane shear along the width at the other
end. It has two circular segments that are connected at the tangent point. The
smaller segment has a mean radius of 0.3556 m (14 inches) and spans 60° from
the clamped end to the tangent point. The larger segment spans 150° from the
tangent point to the free end and has a mean radius of 1.1684 m (46 inches).
The hook is 0.0508 m (2 inches) thick and 0.508 m (20 inches) wide, modeled as
linear elastic with an elastic modulus of 22.77 MPa (3300 psi) and a Poisson's
ratio of 0.35. In most tests the shear force is applied through the use of a
distributing coupling constraint. The coupling constraint provides coupling
between a reference node on which the load is prescribed and the nodes located
on the free end. The distributed nodal loads on the free end are equivalent to
a uniformly distributed load of 8.7563 N/m (0.05 lb/in). In two of the tests an
equivalent shear force is applied as a distributed shear traction instead.


"""
module raasch_dsg3_verification

using Test
using LinearAlgebra
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!

function _execute(input = "raasch_s4_1x9.inp", visualize = true)
    E = 3300.0
    nu = 0.35
    thickness = 2.0
    tolerance = thickness / 2
    # analytical solution for the vertical deflection under the load
    analyt_sol = 5.02
    R = 46.0
    formul = FEMMShellT3FFModule

    output = import_ABAQUS(joinpath(input))
    fens = output["fens"]
    fes = output["fesets"][1]

    connected = findunconnnodes(fens, fes)
    fens, new_numbering = compactnodes(fens, connected)
    fes = renumberconn!(fes, new_numbering)

    fens, fes = Q4toT3(fens, fes)

    # plots = cat(plot_space_box([[0 0 -R/2]; [R/2 R/2 R/2]]),
    #     plot_nodes(fens),
    #     plot_midsurface(fens, fes);
    # dims = 1)
    # pl = render(plots)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Report
    # @info "Raasch hook, formulation=$(formul)"
    # @info "Mesh: $input"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    femm.mult_el_size = 0.2
    femm.drilling_stiffness_scale = 1.0
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # Clamped end
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, l1, true, i)
    end

    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Load
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, box = [97.9615 97.9615 -16 -16 0 20], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, 0, 0.05, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 3)

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    nl = selectnode(fens; box = Float64[97.9615 97.9615 -16 -16 0 0], inflate = tolerance)
    targetu = dchi.values[nl, 3][1]
    # @info "Target: $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

    return targetu / analyt_sol * 100
end

function test_convergence()
    formul = FEMMShellT3FFModule
    results = [91.7059961843231, 95.9355786538892, 97.19276899988246, 98.38896641657374]
    for (m, res) in zip(["1x9", "3x18", "5x36", "10x72"], results)
        v = _execute("raasch_s4_" * m * ".inp", false)
        @test isapprox(res, v, rtol = 1.0e-4)
        # @show v
    end
    return true
end

test_convergence()

end # module

"""
The initially twisted cantilever beam is one of the standard test
problems for verifying the finite-element accuracy [1]. The beam is
clamped at one end and loaded either with unit in-plane or 
unit out-of-plane force at the other. The centroidal axis of the beam is
straight at the undeformed  configuration, while its cross-sections are
twisted about the centroidal axis from 0 at the clamped end to pi/2 at
the free end. 

Reference:
Zupan D, Saje M (2004) On "A proposed standard set of problems to test
finite element accuracy": the twisted beam. Finite Elements in Analysis
and Design 40: 1445-1451.  

Thicker cross section (t = 0.32)
#     Loading in the Z direction
dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
#     Loading in the Y direction
dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;

Thinner cross section (t = 0.0032)
#     Loading in the Z direction
dir = 3; uex = 0.005256;
#     Loading in the Y direction
dir = 2; uex = 0.001294;

MacNeal,  R. H., and R. L. Harder, “A Proposed Standard Set of Problems to Test Finite Element Accuracy,” Finite Elements in Analysis Design, vol. 11, pp. 3–20, 1985.

Simo,  J. C., D. D. Fox, and M. S. Rifai, “On a Stress Resultant Geometrically Exact Shell Model. Part II: The Linear Theory; Computational Aspects,” Computational Methods in Applied Mechanical Engineering, vol. 73, pp. 53–92, 1989.


"""
module twisted_beam_dsg3_verification

using Test
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!

params_thicker_dir_3 = (t = 0.32, force = 1.0, dir = 3, uex = 0.005424534868469)
params_thicker_dir_2 = (t = 0.32, force = 1.0, dir = 2, uex = 0.001753248285256)

params_thinner_dir_3 = (t = 0.0032, force = 1.0e-6, dir = 3, uex = 0.005256)
params_thinner_dir_2 = (t = 0.0032, force = 1.0e-6, dir = 2, uex = 0.001294)


function _execute(
    t = 0.32,
    force = 1.0,
    dir = 3,
    uex = 0.005424534868469,
    nL = 8,
    nW = 2,
    visualize = true,
)
    E = 0.29e8
    nu = 0.22
    W = 1.1
    L = 12.0
    formul = FEMMShellT3FFModule

    tolerance = W / nW / 100
    fens, fes = T3block(L, W, nL, nW, :a)
    fens.xyz = xyz3(fens)
    for i = 1:count(fens)
        a = fens.xyz[i, 1] / L * (pi / 2)
        y = fens.xyz[i, 2] - (W / 2)
        z = fens.xyz[i, 3]
        fens.xyz[i, :] = [fens.xyz[i, 1], y * cos(a) - z * sin(a), y * sin(a) + z * cos(a)]
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), t), mater)
    # femm.mult_el_size = 0.2
    # femm.drilling_stiffness_scale = 1.0
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # Clamped end
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i = 1:6
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Load
    nl = selectnode(fens; box = Float64[L L 0 0 0 0], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    v = FFlt[0, 0, 0, 0, 0, 0]
    v[dir] = force
    fi = ForceIntensity(v)
    F = distribloads(lfemm, geom0, dchi, fi, 3)

    # Solve
    solve_blocked!(dchi, K, F)
    result = dchi.values[nl, dir][1] / uex * 100
    return result
end

function test_convergence()
    formul = FEMMShellT3FFModule
    results = [
        39.709921740907355,
        68.87306876326497,
        86.01944734315117,
        95.04101960524827,
        53.10262177376127,
        83.8593790803426,
        94.91359387874728,
        98.21549248655576,
        48.16757753755567,
        79.43420077873479,
        92.54464819755955,
        96.85008269135751,
        50.577029703967334,
        80.34160167730624,
        92.48675665271801,
        96.7096641005938,
    ]
    for n in [2, 4, 8, 16]
        v = _execute(params_thicker_dir_2..., 2 * n, n, false)
        # println(v, ",")
        @test isapprox(v, popat!(results, 1), rtol = 1.0e-3)
    end
    for n in [2, 4, 8, 16]
        v = _execute(params_thicker_dir_3..., 2 * n, n, false)
        # println(v, ",")
        @test isapprox(v, popat!(results, 1), rtol = 1.0e-3)
    end
    for n in [2, 4, 8, 16]
        v = _execute(params_thinner_dir_2..., 2 * n, n, false)
        # println(v, ",")
        @test isapprox(v, popat!(results, 1), rtol = 1.0e-3)
    end
    for n in [2, 4, 8, 16]
        v = _execute(params_thinner_dir_3..., 2 * n, n, false)
        # println(v, ",")
        @test isapprox(v, popat!(results, 1), rtol = 1.0e-3)
    end
    return true
end

test_convergence()

end # module

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
module LE5_Z_cantilever_dsg3_verification

using Test
using LinearAlgebra
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!


function _execute(input = "nle5xf3c.inp", nrefs = 0, visualize = true)
    E = 210e9
    nu = 0.3
    L = 10.0
    thickness = 0.1
    formul = FEMMShellT3FFModule

    tolerance = thickness / 1000
    output = import_ABAQUS(joinpath(dirname(@__FILE__()), input))
    fens = output["fens"]
    fes = output["fesets"][1]

    connected = findunconnnodes(fens, fes)
    fens, new_numbering = compactnodes(fens, connected)
    fes = renumberconn!(fes, new_numbering)

    for r = 1:nrefs
        fens, fes = T3refine(fens, fes)
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)


    # Report
    # @info "Mesh: $input, nrefs = $nrefs"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    femm.mult_el_size = 0.2
    femm.drilling_stiffness_scale = 1.0
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))


    # plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
    #     plot_nodes(fens),
    #     plot_midsurface(fens, fes);
    # dims = 1)
    # pl = render(plots)

    # Apply EBC's
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1, 2, 3]
        setebc!(dchi, l1, true, i)
    end

    applyebc!(dchi)
    numberdofs!(dchi)
    # @show dchi.nfreedofs

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Load
    nl = selectnode(fens; box = Float64[10.0 10.0 1.0 1.0 0 0], tolerance = tolerance)
    loadbdry1 = FESetP1(reshape(nl, 1, 1))
    lfemm1 = FEMMBase(IntegDomain(loadbdry1, PointRule()))
    fi1 = ForceIntensity(FFlt[0, 0, +0.6e6, 0, 0, 0])
    nl = selectnode(fens; box = Float64[10.0 10.0 -1.0 -1.0 0 0], tolerance = tolerance)
    loadbdry2 = FESetP1(reshape(nl, 1, 1))
    lfemm2 = FEMMBase(IntegDomain(loadbdry2, PointRule()))
    fi2 = ForceIntensity(FFlt[0, 0, -0.6e6, 0, 0, 0])
    F =
        distribloads(lfemm1, geom0, dchi, fi1, 3) +
        distribloads(lfemm2, geom0, dchi, fi2, 3)

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    return minimum(dchi.values[:, 3]), maximum(dchi.values[:, 3])
end

function test_convergence()
    formul = FEMMShellT3FFModule
    # @info "LE5 Z-cantilever, formulation=$(formul)"
    for n in [0]
        res = _execute("nle5xf3c.inp", n, false)
        @test res[1] ≈ (-0.01568028415401719, 0.015401663810490641)[1]
        @test res[2] ≈ (-0.01568028415401719, 0.015401663810490641)[2]
    end
    return true
end

# res = test_st3dsg("nle5xf3c.inp", 0, false)
# @test res[1] ≈ (-0.01568028415401719, 0.015401663810490641)[1]
# @test res[2] ≈ (-0.01568028415401719, 0.015401663810490641)[2]

test_convergence()

end # module


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
centerpoint of the roof moves upward under the self-weight (downwardly directed)
load. Perhaps this is one reason why the problem is not straightforward
numerically. 
"""
module barrel_vault_test

using Test
using LinearAlgebra
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

function _execute(input = "barrelvault_s3r_fineirreg.inp", visualize = true)
    E = 3.0e6
    nu = 0.0
    thickness = 3.0
    tolerance = thickness / 20
    # analytical solution for the vertical deflection under the load
    analyt_sol = -0.3024 * 12
    R = 25.0 * 12
    L = 50.0 * 12
    formul = FEMMShellT3FFModule

    output = import_ABAQUS(joinpath(dirname(@__FILE__()), input))
    fens = output["fens"]
    fes = output["fesets"][1]

    # boundingbox(fens.xyz)
    connected = findunconnnodes(fens, fes)
    fens, new_numbering = compactnodes(fens, connected)
    fes = renumberconn!(fes, new_numbering)

    fens, fes = mergenodes(fens, fes, thickness / 10)

    # plots = cat(plot_space_box([[0 0 -R/2]; [R/2 R/2 R/2]]),
    #     plot_nodes(fens),
    #     plot_midsurface(fens, fes);
    # dims = 1)
    # pl = render(plots)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    # @info "Mesh: $input"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    femm.mult_el_size = 0.2
    femm.drilling_stiffness_scale = 0.1
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # vtkwrite("barrel_vault_test-mesh.vtu", fens, fes)
    # bfes = meshboundary(fes)
    # vtkwrite("barrel_vault_test-boundary-mesh.vtu", fens, bfes)

    # Apply EBC's
    # Rigid diaphragm end
    l1 = selectnode(fens; box = Float64[-Inf Inf -Inf Inf 0 0], inflate = tolerance)
    for i in [1, 2]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Z
    l1 = selectnode(fens; box = Float64[-Inf Inf -Inf Inf L / 2 L / 2], inflate = tolerance)
    for i in [3, 4, 5]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y (apex)
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [2, 4, 6]
        setebc!(dchi, l1, true, i)
    end

    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Load
    # Midpoint of the free edge
    nl = selectnode(
        fens;
        box = Float64[-Inf Inf sin(40 / 360 * 2 * pi) * R sin(40 / 360 * 2 * pi) * R L / 2 L /
                                                                                       2],
        inflate = tolerance,
    )
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt[-0.625, 0, 0, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 3)

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)

    targetu = dchi.values[nl, 1][1]
    # @info "Solution: $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

    # Generate a graphical display of resultants
    cylindrical!(
        csmatout::FFltMat,
        XYZ::FFltMat,
        tangents::FFltMat,
        feid::FInt,
        qpid::FInt,
    ) = begin
        r = -vec(XYZ)
        r[3] = 0.0
        csmatout[:, 3] .= vec(r) / norm(vec(r))
        csmatout[:, 2] .= (0.0, 0.0, 1.0)
        cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
        return csmatout
    end
    ocsys = CSys(3, 3, cylindrical!)
    test_results = [
        (-1520.6449167366522, 14.067403309095397)
        (-73.8262145426215, 425.93651541819503)
        (-0.005341121492547284, 953.0929383629322)
    ]
    scalars = []
    for nc = 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
        # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("m$nc", fld.values))
        @test isapprox(minimum(fld.values), test_results[nc][1], rtol = 0.01)
        @test isapprox(maximum(fld.values), test_results[nc][2], rtol = 0.01)
    end
    # vtkwrite("$(input)-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    test_results = [
        (-306.83173146926623, 309.5860993742647),
        (-1011.4832977998705, 2167.0403478574167),
        (-687.3500043290137, 69.38703021678862),
    ]
    scalars = []
    for nc = 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
        # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("n$nc", fld.values))
        @test isapprox(minimum(fld.values), test_results[nc][1], rtol = 0.01)
        @test isapprox(maximum(fld.values), test_results[nc][2], rtol = 0.01)
    end
    # vtkwrite("$(input)-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    test_results = [
        (-13.784764125688811, 21.165312065421237),
        (-7.4216963152916255, 25.679801383392967),
    ]
    scalars = []
    for nc = 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
        # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("q$nc", fld.values))
        @test isapprox(minimum(fld.values), test_results[nc][1], rtol = 0.01)
        @test isapprox(maximum(fld.values), test_results[nc][2], rtol = 0.01)
    end
    # vtkwrite("$(input)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (R / 2) / maximum(abs.(U)) .* U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(
        plot_space_box([[0 0 -R]; [R R R]]),
        plot_nodes(fens),
        plot_midsurface(
            fens,
            fes;
            x = geom0.values,
            u = dchi.values[:, 1:3],
            R = Rfield0.values,
        );
        dims = 1,
    )
    pl = render(plots)
    return true
end

function test_convergence()

    # @info "Scordelis-Lo Abaqus model"

    _execute("barrelvault_stri3_irreg.inp", false)

    # _execute("barrelvault_s3r_fineirreg.inp", false)

    return true
end

test_convergence()

end # module

