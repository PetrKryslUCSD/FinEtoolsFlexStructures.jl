# Pinched cylinder with diagphram supports and concentrated force
module scordelis_lo_dsg3_verification

using Test
using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3DSGOModule
using FinEtoolsFlexStructures.FEMMShellT3DSGModule
using FinEtoolsFlexStructures.FEMMShellT3DSGICModule
using FinEtoolsFlexStructures.FEMMShellCSDSG3Module
using FinEtoolsFlexStructures.FEMMShellIsoPModule
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

using Infiltrator

function test_st3dsg(args...)
    return _execute_dsg_model(FEMMShellT3DSGModule, args...)
end

function test_st3dsgic(args...)
  return _execute_dsg_model(FEMMShellT3DSGICModule, args...)
end

function test_dsg3(args...)
  return _execute_dsg_model(FEMMShellT3DSGOModule, args...)
end

function test_csdsg3(args...)
  return _execute_dsg_model(FEMMShellCSDSG3Module, args...)
end

function _execute_dsg_model(formul, n = 8, visualize = true)
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    analyt_sol=-0.3024;
    # Parameters:
    E=4.32e8;
    nu=0.0;
    thickness = 0.25; # geometrical dimensions are in feet
    R = 25.0;
    L = 50.0;
    
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

function test_convergence(t, results)
    for (n, res) in  zip([4, 8, 10, 12, 16, 24], results)
        v =  t(n, false)
        @test res ≈ v
    end
    return true
end
                                                            
 
test_convergence(test_dsg3, [64.3167759233874,
85.36891297631232,        
90.2486220974747,
93.32599689114275,                                                              
96.82533584021382,                                                              
99.91041942923614])
test_convergence(test_st3dsgic, [68.93400706852107,                             
89.31902528662164,                                                              
94.05006970072284,                                                              
97.03680208373054,                                                              
100.46625686052543,                                                             
103.64317681989752])
test_convergence(test_st3dsg, [66.64995973373549,                       
    85.59194662963257,                                                              
    89.88220350622672,                                                              
    92.52849163603474,                                                              
    95.41700393607626,                                                              
    97.65867120146895])
test_convergence(test_csdsg3, [67.6211603156003,                       
    87.47670927679906,                                                              
    92.07919359232324,                                                              
    94.99649878356443,                                                              
    98.38626631439227,                                                              
    101.66676977204703])

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
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3DSGOModule
using FinEtoolsFlexStructures.FEMMShellT3DSGICModule
using FinEtoolsFlexStructures.FEMMShellT3DSGModule
using FinEtoolsFlexStructures.FEMMShellCSDSG3Module
using FinEtoolsFlexStructures.FEMMShellIsoPModule
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

using Infiltrator

function test_st3dsg(args...)
    return _execute_dsg_model(FEMMShellT3DSGModule, args...)
end

function test_st3dsgic(args...)
  return _execute_dsg_model(FEMMShellT3DSGICModule, args...)
end

function test_dsg3(args...)
  return _execute_dsg_model(FEMMShellT3DSGOModule, args...)
end

function test_csdsg3(args...)
  return _execute_dsg_model(FEMMShellCSDSG3Module, args...)
end

function _execute_dsg_model(formul, input = "raasch_s4_1x9.inp", visualize = true)
    E = 3300.0;
    nu = 0.35;
    thickness  =  2.0;
    tolerance = thickness/2
    # analytical solution for the vertical deflection under the load
    analyt_sol = 5.02;
    R = 46.0;

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
 
function test_convergence(t, results)
    for (m, res) in zip(["1x9", "3x18", "5x36", "10x72"], results)
        v = t("raasch_s4_" * m * ".inp", false)
        @test res ≈ v
    end
    return true
end

 
test_convergence(test_dsg3, [197.07225021939522, 122.29375814892506, 131.7397746827515, 166.9855842073799])
test_convergence(test_st3dsgic, [98.93041135486558, 107.56727708008175, 121.0279541750594, 160.55183716024845])
test_convergence(test_st3dsg, [96.51445461618886, 96.48046436369553, 97.30091092942541, 98.41519147885221])
test_convergence(test_csdsg3, [90.09510193432982, 101.2652927601732, 114.4122860334214, 155.13736181762695])
       
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
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3DSGOModule
using FinEtoolsFlexStructures.FEMMShellT3DSGICModule
using FinEtoolsFlexStructures.FEMMShellT3DSGModule
using FinEtoolsFlexStructures.FEMMShellCSDSG3Module
using FinEtoolsFlexStructures.FEMMShellIsoPModule
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json, plot_triads

params_thicker_dir_3 = (t =  0.32, force = 1.0, dir = 3, uex = 0.005424534868469); 
params_thicker_dir_2 = (t =  0.32, force = 1.0, dir = 2, uex = 0.001753248285256); 

params_thinner_dir_3 = (t =  0.0032, force = 1.0e-6, dir = 3, uex = 0.005256); 
params_thinner_dir_2 = (t =  0.0032, force = 1.0e-6, dir = 2, uex = 0.001294); 


function test_st3dsg(args...)
    return _execute_dsg_model(FEMMShellT3DSGModule, args...)
end

function test_st3dsgic(args...)
  return _execute_dsg_model(FEMMShellT3DSGICModule, args...)
end

function test_dsg3(args...)
  return _execute_dsg_model(FEMMShellT3DSGOModule, args...)
end

function _execute_dsg_model(formul, t = 0.32, force = 1.0, dir = 3, uex = 0.005424534868469, nL = 8, nW = 2, visualize = true)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.0;
    
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

function test_convergence(t, results)
    for n in [2, 4, 8, 16, ]
        v = t(params_thicker_dir_2..., 2*n, n, false)
        @test isapprox(v, popat!(results, 1), rtol = 1.0e-4)
    end
    for n in [2, 4, 8, 16, ]
        v = t(params_thicker_dir_3..., 2*n, n, false)
        @test isapprox(v, popat!(results, 1), rtol = 1.0e-4)
    end
    for n in [2, 4, 8, 16, ]
        v = t(params_thinner_dir_2..., 2*n, n, false)
        @test isapprox(v, popat!(results, 1), rtol = 1.0e-4)
    end
    for n in [2, 4, 8, 16, ]
        v = t(params_thinner_dir_3..., 2*n, n, false)
        @test isapprox(v, popat!(results, 1), rtol = 1.0e-4)
    end
    return true
end

test_convergence(test_st3dsg, [54.64381168194862                                  
73.44131219368768                                  
87.32404900158092                                  
95.46324216331342                                  
79.62634140519276                                  
91.14293575791413                                  
96.58985595218654                                  
98.66116466318049                                  
68.78329703527952
85.06739918101329                                  
92.79651734583757                               
96.51004023953843
79.23340078402887
88.4444050673403
93.60212106778023                                 
96.63919884543652 ])

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
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3DSGOModule
using FinEtoolsFlexStructures.FEMMShellT3DSGICModule
using FinEtoolsFlexStructures.FEMMShellT3DSGModule
using FinEtoolsFlexStructures.FEMMShellCSDSG3Module
using FinEtoolsFlexStructures.FEMMShellIsoPModule
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json



function test_st3dsg(args...)
    return _execute_dsg_model(FEMMShellT3DSGModule, args...)
end

function test_st3dsgic(args...)
  return _execute_dsg_model(FEMMShellT3DSGICModule, args...)
end

function test_dsg3(args...)
  return _execute_dsg_model(FEMMShellT3DSGOModule, args...)
end

function test_csdsg3(args...)
  return _execute_dsg_model(FEMMShellCSDSG3Module, args...)
end

function _execute_dsg_model(formul, input = "nle5xf3c.inp", nrefs = 0, visualize = true)
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
    
    
    # Report
        # @info "Mesh: $input, nrefs = $nrefs"

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

function test_convergence(t)
    @info "LE5 Z-cantilever, formulation=$(t)"
    for n in [0, 1, 2, 3, 4, ]
        t("nle5xf3c.inp", n, false)
    end
    return true
end

res = test_st3dsg("nle5xf3c.inp", 0, false)
@test res[1] ≈ (-0.01568028415401719, 0.015401663810490641)[1]
@test res[2] ≈ (-0.01568028415401719, 0.015401663810490641)[2]

end # module

 