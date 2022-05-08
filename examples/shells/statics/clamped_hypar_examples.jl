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
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

function _execute(tL_ratio = 1/100, g = 80*0.1^0, analyt_sol=-9.3355e-5, n = 32, visualize = false)
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    # Parameters:
    E = 2.0e11;
    nu = 0.3;
    L = 1.0
    thickness = tL_ratio * L
    # Bathe, Iosilevich, and Chapelle (2000) with a refined mesh of
    # high-order element MITC16
    formul = FEMMShellT3FFModule

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
    femm.drilling_stiffness_scale = 0.1
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
    @info "Displacement Solution: $(round(targetu/analyt_sol, digits = 4)*100)%"
    @info "Strain Energy Solution: $(U'*K*U/2)"

        # Visualization
    if visualize
        vtkwrite("clamped_hypar_examples-geometry.vtu", fens, fes)
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

function test_convergence()
    
    tL_ratios = [1/100, 1/1000, 1/10000]; 
    gs = [80*0.1^0, 80*0.1^1, 80*0.1^2]
    analyt_sols = [-9.3355e-5, -6.3941e-3, -5.2988e-1];
    
    for (tL_ratio, g, analyt_sol) in zip(tL_ratios, gs, analyt_sols)
        @info "Clamped hypar, t/L=$(tL_ratio)"
        results = Float64[]
        for n in [4, 8, 16, 32, 64, 128, 256, ]
        # for n in [4, 8, 16, 32, 64, ]
            r = _execute(tL_ratio, g, analyt_sol, n, false)
            push!(results, r/analyt_sol)
        end   
        @show results
    end

    return true
end

# using Gnuplot

# ns = 1 ./ [4, 8, 16, 32, 64, 128, 256, ]
# results = [0.8839348674712617, 0.8750949157612452, 0.9199805802913757, 0.9619508790573108, 0.9856803572479892, 0.9955727622499687, 0.9993169031485688]   
# @gp ns results "with lp"      :-           
# results = [1.0429797613488236, 0.9314984628085947, 0.9365905225801154, 0.9565506281799385, 0.9764476699285441, 0.9902805329751646, 0.9968920296205528]   
# @gp :- ns results "with lp"    :-                             
# results = [1.251888013432877, 1.0155533090845452, 0.9678060658415124, 0.9718188010061173, 0.9812934066246979, 0.9889499817887738, 0.9953521405300628] 
# @gp :- ns results "with lp"

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