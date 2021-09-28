# Pinched cylinder with diagphram supports and concentrated force
module scordelis_lo_dsg3_verification

using Test
using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3DSGAModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

using Infiltrator

function _execute(n = 8, visualize = true)
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    analyt_sol=-0.3024;
    # Parameters:
    E=4.32e8;
    nu=0.0;
    thickness = 0.25; # geometrical dimensions are in feet
    R = 25.0;
    L = 50.0;
    formul = FEMMShellT3DSGAModule
    
    tolerance = R/n/1000
    fens, fes = T3block(40/360*2*pi,L/2,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*(cos(a)-1))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    femm.drilling_stiffness_scale = 1.0
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

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[sin(40/360*2*pi)*25 sin(40/360*2*pi)*25 L/2 L/2 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt[0, 0, -90, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    resultpercent = dchi.values[nl, 3][1]/analyt_sol*100

    return resultpercent
end

function test_convergence()
    formul = FEMMShellT3DSGAModule
    results = [
    66.6606197093463,                                                
    85.5997636990392,                                                
    89.88947890745521,                                               
    92.53516614150621,                                               
    95.4225137166831, 
    97.66243080996651
    ]
    for (n, res) in  zip([4, 8, 10, 12, 16, 24], results)
        v = _execute(n, false)
        @test isapprox(res, v, rtol = 1.0e-4)
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
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3DSGAModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

using Infiltrator

function _execute(input = "raasch_s4_1x9.inp", visualize = true)
    E = 3300.0;
    nu = 0.35;
    thickness  =  2.0;
    tolerance = thickness/2
    # analytical solution for the vertical deflection under the load
    analyt_sol = 5.02;
    R = 46.0;
    formul = FEMMShellT3DSGAModule

    output = import_ABAQUS(joinpath(input))
    fens = output["fens"]
    fes = output["fesets"][1]

    connected = findunconnnodes(fens, fes);
    fens, new_numbering = compactnodes(fens, connected);
    fes = renumberconn!(fes, new_numbering);

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
    femm.drilling_stiffness_scale = 1.0
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped end
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, box = [97.9615 97.9615 -16 -16 0 20], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, 0, 0.05, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    nl = selectnode(fens; box = Float64[97.9615 97.9615 -16 -16 0 0], inflate = tolerance)
    targetu =  dchi.values[nl, 3][1]
    # @info "Target: $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

    return targetu/analyt_sol*100
end
 
function test_convergence()
    formul = FEMMShellT3DSGAModule
    results = [
    95.96774650732868, 
    96.68065217811291, 
    97.39437142020671, 
    98.4519858066279
    ]
    for (m, res) in zip(["1x9", "3x18", "5x36", "10x72"], results)
        v = _execute("raasch_s4_" * m * ".inp", false)
        @test isapprox(res, v, rtol = 1.0e-4)
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
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3DSGAModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json, plot_triads

params_thicker_dir_3 = (t =  0.32, force = 1.0, dir = 3, uex = 0.005424534868469); 
params_thicker_dir_2 = (t =  0.32, force = 1.0, dir = 2, uex = 0.001753248285256); 

params_thinner_dir_3 = (t =  0.0032, force = 1.0e-6, dir = 3, uex = 0.005256); 
params_thinner_dir_2 = (t =  0.0032, force = 1.0e-6, dir = 2, uex = 0.001294); 


function _execute(t = 0.32, force = 1.0, dir = 3, uex = 0.005424534868469, nL = 8, nW = 2, visualize = true)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.0;
    formul = FEMMShellT3DSGAModule
    
    tolerance = W/nW/100
    fens, fes = T3block(L,W,nL,nW,:a);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i,1]/L*(pi/2); y=fens.xyz[i,2]-(W/2); z=fens.xyz[i,3];
        fens.xyz[i,:]=[fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), t), mater)
    femm.drilling_stiffness_scale = 1.0
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped end
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in 1:6
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[L L 0 0 0 0], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    v = FFlt[0, 0, 0, 0, 0, 0]
    v[dir] = force
    fi = ForceIntensity(v);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    result =  dchi.values[nl, dir][1]/uex*100
    return result
end

function test_convergence()
    formul = FEMMShellT3DSGAModule
    results = [
    60.729833319894325,
    75.03868479181727,
    87.50963379377527,
    95.47124053670485,
    84.65337521150289,
    91.88041248379042,
    96.63969262242233,
    98.65588185555086,
    78.01013136700996,
    88.55708071163137,
    93.80585717074035,
    96.77822361996257,
    85.16655215203586,
    89.84038466760585,
    93.94781632610085,
    96.7274995250634,
    ]
    for n in [2, 4, 8, 16, ]
        v = _execute(params_thicker_dir_2..., 2*n, n, false)
        # @show v
        @test isapprox(v, popat!(results, 1), rtol = 1.0e-3)
    end
    for n in [2, 4, 8, 16, ]
        v = _execute(params_thicker_dir_3..., 2*n, n, false)
        # @show v
        @test isapprox(v, popat!(results, 1), rtol = 1.0e-3)
    end
    for n in [2, 4, 8, 16, ]
        v = _execute(params_thinner_dir_2..., 2*n, n, false)
        # @show v
        @test isapprox(v, popat!(results, 1), rtol = 1.0e-3)
    end
    for n in [2, 4, 8, 16, ]
        v = _execute(params_thinner_dir_3..., 2*n, n, false)
        # @show v
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
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3DSGAModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json



function _execute(input = "nle5xf3c.inp", nrefs = 0, visualize = true)
    E = 210e9;
    nu = 0.3;
    L = 10.0;
    thickness = 0.1
    formul = FEMMShellT3DSGAModule

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
    
    
    # Report
        # @info "Mesh: $input, nrefs = $nrefs"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    femm.drilling_stiffness_scale = 1.0
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
    # @show dchi.nfreedofs

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    
    # Load
    nl = selectnode(fens; box = Float64[10.0 10.0 1.0 1.0 0 0], tolerance = tolerance)
    loadbdry1 = FESetP1(reshape(nl, 1, 1))
    lfemm1 = FEMMBase(IntegDomain(loadbdry1, PointRule()))
    fi1 = ForceIntensity(FFlt[0, 0, +0.6e6, 0, 0, 0]);
    nl = selectnode(fens; box = Float64[10.0 10.0 -1.0 -1.0 0 0], tolerance = tolerance)
    loadbdry2 = FESetP1(reshape(nl, 1, 1))
    lfemm2 = FEMMBase(IntegDomain(loadbdry2, PointRule()))
    fi2 = ForceIntensity(FFlt[0, 0, -0.6e6, 0, 0, 0]);
    F = distribloads(lfemm1, geom0, dchi, fi1, 3) + distribloads(lfemm2, geom0, dchi, fi2, 3);

    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    return minimum(dchi.values[:, 3]), maximum(dchi.values[:, 3])
end

function test_convergence()
    formul = FEMMShellT3DSGAModule
    # @info "LE5 Z-cantilever, formulation=$(formul)"
    for n in [0,  ]
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

 