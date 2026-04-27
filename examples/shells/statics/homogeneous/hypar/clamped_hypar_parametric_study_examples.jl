"""
From: CMES, vol.49, no.2, pp.81-110, 2009

The problem considered in this section is that of a hyperbolic paraboloid shell,
clamped along one side and free on three edges and loaded by self-weight
(Figure 16). This is a pure bending dominated problem and known to be a very
hard test for locking behaviour as suggested in References (Chapelle and Bathe,
1998; Bathe, Iosilevich, and Chapelle, 2000). 

Table 2 in Bathe, Iosilevich, and Chapelle, 2000 (computed with MITC16).
mult_el_size   Strain energy  Displacement
1/100 1.6790e-3 9.3355e-5
1/1000 1.1013e-2 6.3941e-3
1/10,000 8.9867e-2 5.2988e-1

See also the table in  
An Improved Quadrilateral Flat Element with Drilling
Degrees of Freedom for Shell Structural Analysis
H. Nguyen-Van , N. Mai-Duy and T. Tran-Cong
CMES, vol.49, no.2, pp.81-110, 2009.

The shell geometry is described by the
equation: z = x^2 −y^2 ; (x,y) ∈ − L/2 ; L/2
"""
module clamped_hypar_parametric_study

using PGFPlotsX
using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.FEMMShellQ4RSModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

# Parameters:
const E = 2.0e11;
const nu = 0.3
const L = 1.0

function _execute_t3ff(tL_ratio = 1/100, g = 80*0.1^0, mult_el_size = 1.0, n = 32, visualize = false)
    thickness = tL_ratio * L
    formul = FEMMShellT3FFModule

    # @info "Mesh: $n elements per side"

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
    femm.mult_el_size = mult_el_size
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

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    Kff = stiffness(femm, massem, geom0, u0, Rfield0, dchi);

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
    # fi = ForceIntensity(Float64, 6, computeforce!)
    fi = ForceIntensity(Float64[0, 0, -g, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 3);
    
    # Solve
    Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)
    
    targetu = dchi.values[nl, 3][1]
    targetenergy = Uf'*Ff/2
    # @info "Displacement Solution: $(round(targetu/analyt_sol, digits = 4)*100)%"
    # @info "Strain Energy Solution: $(Uf'*Ff/2)"

        # Visualization
    if visualize
        vtkwrite("clamped_hypar_parametric_study-geometry.vtu", fens, fes)
        scattersysvec!(dchi, (L/4)/abs(targetu).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end
    return targetu, targetenergy
end

function _execute_q4rs(tL_ratio = 1/100, g = 80*0.1^0, mult_el_size = 1.0, n = 32, visualize = false)
    thickness = tL_ratio * L
    formul = FEMMShellQ4RSModule

    # @info "Mesh: $n elements per side"

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
    femm = formul.make(IntegDomain(fes, 
            GaussRule(2, 2),
            thickness), mater)
    femm.mult_el_size = mult_el_size
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

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    Kff = stiffness(femm, massem, geom0, u0, Rfield0, dchi);

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
    # fi = ForceIntensity(Float64, 6, computeforce!)
    fi = ForceIntensity(Float64[0, 0, -g, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 3);
    
    # Solve
    Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)
    
    targetu = dchi.values[nl, 3][1]
    targetenergy = Uf'*Ff/2
    # @info "Displacement Solution: $(round(targetu/analyt_sol, digits = 4)*100)%"
    # @info "Strain Energy Solution: $(Uf'*Ff/2)"

        # Visualization
    if visualize
        vtkwrite("clamped_hypar_parametric_study-geometry.vtu", fens, fes)
        scattersysvec!(dchi, (L/4)/abs(targetu).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end
    return targetu, targetenergy
end

const tL_ratios = [1/100, 1/1000, 1/10000]; 
const gs = [80*0.1^0, 80*0.1^1, 80*0.1^2]
const u_sols = [-9.3355e-5, -6.3941e-3, -5.2988e-1];
const energy_sols = [1.6790e-3, 1.1013e-2, 8.9867e-2];
const case_data = zip(tL_ratios, gs, u_sols, energy_sols)
const ns = [4, 8, 16, 32, 64, 128]
const mult_el_size = [0.0, 0.1, 0.2, 0.5, 1.0]
const colors = ["black", "red", "blue", "green", "magenta"]

function test_convergence_t3ff(tL_ratio, g, u_sol, energy_sol)
    @info "Clamped hypar, T3FF elements, t/L=$(tL_ratio)"
    overall_results = []
    for mult_el_size in [0.0, 0.1, 0.2, 0.5, 1.0]
        results = []
        for n in ns
            r = _execute_t3ff(tL_ratio, g, mult_el_size, n, false)
            push!(results, r)
        end   
        push!(overall_results, (mult_el_size = mult_el_size, results = results))
    end
    return overall_results
end

function start_case()
    objects = []
    return objects
end

function plot_case!(objects, j, nu_errs)
    @pgf p = PGFPlotsX.Plot(
        {
            color = colors[j],
            line_width = 1.0,
        },
        Coordinates([v for v in zip(L ./ ns, abs.(nu_errs))])
    )
    push!(objects, p)
    push!(objects, LegendEntry("$(mult_el_size[j])"))
end

function display_case(tL_ratio, objects)
    @pgf ax = Axis(
        {
            xlabel = "Normalized element size [ND]",
            ylabel = "Normalized approximate error [ND]",
            # ymin = -1.1e-5,
            # ymax = 1.1e-5,
            # xmin = 0.0,
            # xmax = 1.0,
            xmode = "log",
            ymode = "log",
            yminorgrids = "true",
            grid = "both",
            legend_style = {
                at = Coordinate(0.5, 1.05),
                anchor = "south",
                legend_columns = -1
            },
        },
        objects...
    )

    display(ax)
    pgfsave("tL_ratio=$(tL_ratio).pdf", ax)
end

function test_convergence_all()
    # for c in case_data
    #     @info "--------------------------------------------------"
    #     @info "Convergence study for case t/L=$(c[1]) "
    #     @show c
    #     or = test_convergence_t3ff(c...)
    #     for (m, r) in or
    #         @info "  mult_el_size=$(m)"
    #         aprox_u_sols = [r[i][1] for i in eachindex(r)]
    #         approx_energy_sols = [r[i][2] for i in eachindex(r)]
    #         nu_errs = [abs(aprox_u_sols[i+1] - aprox_u_sols[i]) / aprox_u_sols[end] for i in 1:length(aprox_u_sols)-1]
    #         @info "  Displ. solutions: $(round.(aprox_u_sols, digits = 7))"
    #         # @info "  Displ. norm. solutions: $(round.(aprox_u_sols ./ c[3], digits = 4))"
    #         @info "  Displ. norm. approximate errors: $(round.(nu_errs, digits = 4))"
    #         # @info "  Energy convergence: $(nenergy_sols)"
    #     end
    # end
    for c in case_data
        @info "--------------------------------------------------"
        @info "Convergence study for case t/L=$(c[1]) "
        or = test_convergence_q4rs(c...)
        objects = start_case()
        for (j, (m, r)) in enumerate(or)
            @info "  mult_el_size=$(m)"
            aprox_u_sols = [r[i][1] for i in eachindex(r)]
            approx_energy_sols = [r[i][2] for i in eachindex(r)]
            nu_errs = [abs(aprox_u_sols[i+1] - aprox_u_sols[i]) / aprox_u_sols[end] for i in 1:length(aprox_u_sols)-1]
            @info "  Displ. solutions: $(round.(aprox_u_sols, digits = 7))"
            # @info "  Displ. norm. solutions: $(round.(aprox_u_sols ./ c[3], digits = 4))"
            @info "  Displ. norm. approximate errors: $(round.(nu_errs, digits = 4))"
            # @info "  Energy convergence: $(nenergy_sols)"
            plot_case!(objects, j, nu_errs)
        end
        display_case(c[1], objects)
    end
    return true
end

function test_convergence_q4rs(tL_ratio, g, u_sol, energy_sol)
    @info "Clamped hypar, Q4RS elements, t/L=$(tL_ratio)"
    overall_results = []
    for mult_el_size in [0.0, 0.1, 0.2, 0.5, 1.0]
        results = []
        for n in ns
            r = _execute_q4rs(tL_ratio, g, mult_el_size, n, false)
            push!(results, r)
        end   
        push!(overall_results, (mult_el_size = mult_el_size, results = results))
    end
    return overall_results
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
    println("# test_convergence_all ")
    test_convergence_all()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
