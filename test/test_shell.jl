# Pinched cylinder with diagphram supports and concentrated force
module scordelis_lo_test

using Test
using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellDSG3Module
using FinEtoolsFlexStructures.FEMMShellDSG3IFModule
using FinEtoolsFlexStructures.FEMMShellDSG3IModule
using FinEtoolsFlexStructures.FEMMShellCSDSG3Module
using FinEtoolsFlexStructures.FEMMShellIsoPModule
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

using Infiltrator

function test_dsg3if(args...)
    return _execute_dsg_model(FEMMShellDSG3IFModule, args...)
end

function test_dsg3i(args...)
  return _execute_dsg_model(FEMMShellDSG3IModule, args...)
end

function test_dsg3(args...)
  return _execute_dsg_model(FEMMShellDSG3Module, args...)
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
test_convergence(test_dsg3i, [68.93400706852107,                             
89.31902528662164,                                                              
94.05006970072284,                                                              
97.03680208373054,                                                              
100.46625686052543,                                                             
103.64317681989752])
test_convergence(test_dsg3if, [66.68238024904946,                       
85.61660951581757,                                                              
89.90356898669276,                                                              
92.54697050425355,                                                              
95.43101599843146,                                                              
97.66729527571941])
test_convergence(test_csdsg3, [67.6211603156003,                       
    87.47670927679906,                                                              
    92.07919359232324,                                                              
    94.99649878356443,                                                              
    98.38626631439227,                                                              
    101.66676977204703])
end # module
 