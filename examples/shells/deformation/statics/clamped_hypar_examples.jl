"""
From: CMES, vol.49, no.2, pp.81-110, 2009

The problem considered in this section is that of a hyperbolic paraboloid shell,
clamped along one side and free on three edges and loaded by self-weight
(Figure 16). This is a pure bending dominated problem and known to be a very
hard test for locking behaviour as suggested in References (Chapelle and Bathe,
1998; Bathe, Iosilevich, and Chapelle, 2000). 

Table 2 in Bathe, Iosilevich, and Chapelle, 2000 (computed with MITC16).
t/L   Strain energy  Displacement
1/100 1.6790e-3 9.3355e-5
1/1000 1.1013e-2 6.3941e-3
1/10,000 8.9867e-2 5.2988e-1

See also the table in  
An Improved Quadrilateral Flat Element with Drilling
Degrees of Freedom for Shell Structural Analysis
H. Nguyen-Van1 , N. Mai-Duy1 and T. Tran-Cong1
CMES, vol.49, no.2, pp.81-110, 2009.

The shell geometry is described by the
equation: z = x^2 −y^2 ; (x,y) ∈ − L/2 ; L/2
"""
module clamped_hypar_examples

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3ODSGModule
using FinEtoolsFlexStructures.FEMMShellT3DSGICModule
using FinEtoolsFlexStructures.FEMMShellT3DSGModule
# using FinEtoolsFlexStructures.FEMMShellT3Module
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function _execute_dsg_model(formul, tL_ratio = 1/100, g = 80*0.1^0, analyt_sol=-9.3355e-5, n = 32, visualize = false)
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    # Parameters:
    E = 2.0e11;
    nu = 0.3;
    L = 1.0
    thickness = tL_ratio * L
    # Bathe, Iosilevich, and Chapelle (2000) with a refined mesh of
    # high-order element MITC16

    @info "Mesh: $n elements per side"

    tolerance = L/n/1000
    fens, fes = T3block(L,L,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        x=fens.xyz[i, 1]-L/2; y=fens.xyz[i, 2]-L/2;
        fens.xyz[i, :] .= (x, y, x^2 - y^2)
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
    # Clamped edge
    l1 = selectnode(fens; box = Float64[-L/2 -L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[L/2 L/2 0 0 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    # computeforce!(forceout, XYZ, tangents, fe_label) = let
    #     n = cross(tangents[:, 1], tangents[:, 2])
    #     n = n / norm(n)
    #     forceout .= 0.0
    #     forceout[3] = -g * n[3]
    #     return forceout
    # end
    # fi = ForceIntensity(FFlt, 6, computeforce!)
    fi = ForceIntensity(FFlt[0, 0, -g, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    targetu = dchi.values[nl, 3][1]
    @info "Solution: $(round(targetu/analyt_sol, digits = 4)*100)%"

        # Visualization
    if visualize
        scattersysvec!(dchi, (L/4)/abs(targetu).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end
    return targetu
end

function test_st3dsg(args...)
    return _execute_dsg_model(FEMMShellT3DSGModule, args...)
end

function test_st3dsgic(args...)
  return _execute_dsg_model(FEMMShellT3DSGICModule, args...)
end

function test_dsg3(args...)
  return _execute_dsg_model(FEMMShellT3ODSGModule, args...)
end

function test_csdsg3(args...)
  return _execute_dsg_model(FEMMShellCSDSG3Module, args...)
end

function test_q4sri(args...)
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    # Parameters:
    E = 2.0e11;
    nu = 0.3;
    L = 1.0
    # Bathe, Iosilevich, and Chapelle (2000) with a refined mesh of
    # high-order element MITC16

    thickness = L/100; 
    analyt_sol=-9.3355e-5;
    g = 80*0.1^0

    thickness = L/1000; 
    analyt_sol=-6.3941e-5;
    g = 80*0.1^3
    
    formul = FEMMShellQ4SRIModule
    # Report
    @info "Clamped hypar, t/L=$(thickness/L), formulation=$(formul)"
    @info "Mesh: $n elements per side"

    tolerance = L/n/1000
    fens, fes = Q4block(L,L,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        x=fens.xyz[i, 1]-L/2; y=fens.xyz[i, 2]-L/2;
        fens.xyz[i, :] .= (x, y, x^2 - y^2)
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.FEMMShellQ4SRI(IntegDomain(fes, GaussRule(2, 2), thickness), mater)
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped edge
    l1 = selectnode(fens; box = Float64[-L/2 -L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[L/2 L/2 0 0 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    # computeforce!(forceout, XYZ, tangents, fe_label) = let
    #     n = cross(tangents[:, 1], tangents[:, 2])
    #     n = n / norm(n)
    #     forceout .= 0.0
    #     forceout[3] = -g * n[3]
    #     return forceout
    # end
    # fi = ForceIntensity(FFlt, 6, computeforce!)
    fi = ForceIntensity(FFlt[0, 0, -g, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    targetu = dchi.values[nl, 3][1]
    @info "Target: $targetu,  $(round(targetu/analyt_sol, digits = 4)*100)%"

    # Visualization
    scattersysvec!(dchi, (L/4)/abs(targetu).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_convergence(t)
    
    tL_ratio = 1/100; 
    g = 80*0.1^0
    analyt_sol=-9.3355e-5;

    # tL_ratio = 1/1000; 
    # g = 80*0.1^3#
    # analyt_sol=-6.3941e-5;

    @info "Clamped hypar, t/L=$(tL_ratio), formulation=$(t)"
       
    for n in [2, 4, 8, 16, 32, 64]
        t(tL_ratio, g, analyt_sol, n, false)
    end
    return true
end


function test_dependence_on_thickness(t)
    
    tL_ratios = [1/100, 1/1000, 1/10000]; 
    gs = [80*0.1^0, 80*0.1^1, 80*0.1^2]
    analyt_sols = [-9.3355e-5, -6.3941e-3, -5.2988e-1];

    @info "Clamped hypar, t/L=$(tL_ratio), formulation=$(t)"
    
    for (tL_ratio, g, analyt_sol) in zip(tL_ratios, gs, analyt_sols)
        results = []
        for n in [4, 8, 16, 32, 64, 128, 256, ]
            r = t(tL_ratio, g, analyt_sol, n, false)
            push!(results, r/analyt_sol)
        end   
        @show results
    end

    return true
end

# using Gnuplot

# ns = 1 ./ [4, 8, 16, 32, 64, 128, 256, 512]
# results = [0.8571252448599641, 0.8665242340654647, 0.9178272403908829, 0.9613361100770298, 0.9854577511292909, 0.9954878630275767, 0.9992902287331279, 1.0007228233464933]    
# @gp ns results "with lp"      :-           
# results = [1.0075109452933262, 0.918680499178048, 0.9336875438959243, 0.9558938174582949, 0.9762891817704417, 0.9902400506734798, 0.9968806281803028, 0.9993845760976003]  
# @gp :- ns results "with lp"    :-                             
# results = [1.2091137985344178, 1.0002021765631814, 0.9641514478788299, 0.971027447847773, 0.9811117321314542, 0.9889056459333478, 0.9953395935946254, 0.999163974846612]  
# @gp :- ns results "with lp"


function allrun()
    println("#####################################################")
    println("# test_dsg3 ")
    test_dsg3()
    println("#####################################################")
    println("# test_q4sri ")
    test_q4sri()
    println("#####################################################")
    println("# test_dsg3_convergence  ")
    test_dsg3_convergence()
    println("#####################################################")
    println("# test_q4sri_convergence  ")
    test_q4sri_convergence()
    return true
end # function allrun

end # module

using .clamped_hypar_examples
m = clamped_hypar_examples
# m.test_convergence(m.test_dsg3)
# m.test_convergence(m.test_st3dsgic)
# m.test_convergence(m.test_st3dsg)
m.test_dependence_on_thickness(m.test_st3dsgic)