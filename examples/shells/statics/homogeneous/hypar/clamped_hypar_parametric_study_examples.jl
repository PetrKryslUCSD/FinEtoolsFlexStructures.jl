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

mutable struct _CSEval{T}
    tang::Matrix{T}
end

function (o::_CSEval)(csmatout, XYZ, tangents, feid, qpid)
    o.tang[:, 1] .= (1.0, 0.0, 2.0 .* XYZ[1])
    o.tang[:, 2] .= (0.0, 1.0, -2.0 .* XYZ[2])
    cross3!(view(csmatout, :, 3), view(o.tang, :, 1), view(o.tang, :, 2))
    csmatout[:, 3] .= vec(view(csmatout, :, 3))/norm(vec(view(csmatout, :, 3)))
    csmatout[:, 1] .= o.tang[:, 1]
    csmatout[:, 1] .= vec(view(csmatout, :, 1))/norm(vec(view(csmatout, :, 1)))
    cross3!(view(csmatout, :, 2), view(csmatout, :, 3), view(csmatout, :, 1))
    return csmatout
end


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

function _execute_q4rs(tL_ratio = 1/100, g = 80*0.1^0, mes_fun=:linfract, mult_el_size = 1.0, n = 32, visualize = false)
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
    
    cse = _CSEval(zeros(3, 2))
    ocsys = CSys(3, 3, (csmatout, XYZ, tangents, feid, qpid) -> cse(csmatout, XYZ, tangents, feid, qpid))
    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, 
            GaussRule(2, 2),
            thickness), ocsys, mater,
            if mes_fun == :linfract
                (t, h) -> 1 / (1 + mult_el_size * h/t)
            elseif mes_fun == :powfract
                (t, h) -> 1 / (1 + mult_el_size * (h/t)^(5/4))
            elseif mes_fun == :sqrtlinfract
                (t, h) -> 1 / (1 + mult_el_size * sqrt(h/t))
            elseif mes_fun == :sqrtquadfract
                (t, h) -> 1 / sqrt(1 + mult_el_size * (h/t)^2)
            elseif mes_fun == :quadfract
                (t, h) -> 1 / (1 + mult_el_size * (h/t)^2)
            end
            # (t, h) -> 1 / (1 + mult_el_size * h^2/t^2))
            # (t, h) -> t^2 / (t^2 + mult_el_size * h^2))
    )
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
    vtkwrite("debug-normals.vtu", fens, fes; vectors = [("normals", deepcopy(femm._normals[:, 1:3]))])
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

struct CaseData
    tL_ratio::Float64
    g::Float64
    u_sol::Float64
    energy_sol::Float64
end


const case_data = [
    CaseData(1/100, 80*0.1^0, -9.3355e-5, 1.6790e-3),
    CaseData(1/1000, 80*0.1^1, -6.3941e-3, 1.1013e-2),
    CaseData(1/10000, 80*0.1^2, -5.2988e-1, 8.9867e-2)
    ]

const colors = ["black", "red", "blue", "green", "magenta", "cyan"]
marks = ["x", "o", "square", "+", "diamond", "star"]

function start_case()
    objects = []
    return objects
end

function plot_case_approx_error!(objects, j, nu_errs, ns, mult_el_sizes)
    @pgf p = PGFPlotsX.Plot(
        {
            color = colors[j],
            mark = "$(marks[j])",
            line_width = 1.0,
        },
        Coordinates([v for v in zip(L ./ ns, abs.(nu_errs))])
    )
    push!(objects, p)
    push!(objects, LegendEntry("$(mult_el_sizes[j])"))
end

function display_case_approx_error(objects, tL_ratio, mes_fun)
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
    pgfsave("approx_error-$(mes_fun)-tL_ratio=$(tL_ratio).pdf", ax)
end

function test_convergence_all(cases, ns = [4, 8, 16, 32, 64, 128, ], mes_fun = :linfract, mult_el_sizes = [0.0, 0.05, 0.1, 0.2, 0.4])
    for c in cases
        @info "--------------------------------------------------"
        @info "Convergence study for case t/L=$(c.tL_ratio) "
        or = _test_convergence_q4rs(c, ns, mes_fun, mult_el_sizes)
        objects_errors = start_case()
        objects_displacements = start_case()
        for (j, (_, m, r)) in enumerate(or)
            @info "  mult_el_size=$(m)"
            aprox_u_sols = [r[i][1] for i in eachindex(r)]
            approx_energy_sols = [r[i][2] for i in eachindex(r)]
            nu_errs = [abs(aprox_u_sols[i+1] - aprox_u_sols[i]) / aprox_u_sols[end] for i in 1:length(aprox_u_sols)-1]
            nus = [aprox_u_sols[i] / c.u_sol for i in 1:length(aprox_u_sols)]
            @info "  Displ. solutions: $(round.(aprox_u_sols, digits = 7))"
            # @info "  Displ. norm. solutions: $(round.(aprox_u_sols ./ c[3], digits = 4))"
            @info "  Displ. norm. approximate errors: $(round.(nu_errs, digits = 4))"
            # @info "  Energy convergence: $(nenergy_sols)"
            plot_case_approx_error!(objects_errors, j, nu_errs, ns, mult_el_sizes)
            plot_case_result!(objects_displacements, j, nus, ns, mult_el_sizes)
        end
        display_case_approx_error(objects_errors, c.tL_ratio, mes_fun)
        display_case_result(objects_displacements, c.tL_ratio, mes_fun)
    end
    return true
end

function _test_convergence_q4rs(c, ns, mes_fun, mult_el_sizes)
    @info "Clamped hypar, Q4RS elements, t/L=$(c.tL_ratio)"
    all_results = []
    for mult_el_size in mult_el_sizes
        results = []
        for n in ns
            r = _execute_q4rs(c.tL_ratio, c.g, mes_fun, mult_el_size, n, false)
            push!(results, r)
        end   
        push!(all_results, (case = c, mult_el_size = mult_el_size, results = results))
    end
    return all_results
end

function plot_case_result!(objects, j, res, ns, mult_el_sizes)
    @pgf p = PGFPlotsX.Plot(
        {
            color = colors[j],
            mark = "$(marks[j])",
            line_width = 1.0,
        },
        Coordinates([v for v in zip(L ./ ns, res)])
    )
    push!(objects, p)
    push!(objects, LegendEntry("$(mult_el_sizes[j])"))
end

function display_case_result(objects, tL_ratio, mes_fun)
    @pgf ax = Axis(
        {
            xlabel = "Normalized element size [ND]",
            ylabel = "Result, normalized [ND]",
            # ymin = -1.1e-5,
            # ymax = 1.1e-5,
            # xmin = 0.0,
            # xmax = 1.0,
            xmode = "log",
            # ymode = "log",
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
    pgfsave("result-$(mes_fun)-tL_ratio=$(tL_ratio).pdf", ax)
end

using FinEtools.AlgoBaseModule: richextrapol

function try_rich_all(c, ns = [4, 8, 16, 32, 64, 128, ], mes_fun = :linfract, mult_el_sizes = [0.0, 0.05, 0.1, 0.2, 0.4])
    @info "--------------------------------------------------"
    @info "Convergence study for case t/L=$(c.tL_ratio), ns=$(ns), u_sol=$(c.u_sol), energy_sol=$(c.energy_sol)"
    or = _test_convergence_q4rs(c, ns, mes_fun, mult_el_sizes)
    _try_rich(or)
    return true
end

function _try_rich(all_results)
    start = 1
    for j in eachindex(all_results)
        c = all_results[j].case
        mult_el_size = all_results[j].mult_el_size
        results = all_results[j].results
        aprox_u_sols = [results[s][1] for s in start:start+2]
        try
            er = richextrapol(aprox_u_sols, [4.0, 2.0, 1.0])
            @info "mult_el_size=$(mult_el_size): $(er[1]), i.e. $(round(er[1]/c.u_sol, digits = 4)*100)%"
        catch e
            @warn "Richardson extrapolation failed for mult_el_size=$(mult_el_size)\n  $(aprox_u_sols)"
        end
    end
end

function allrun()
    ns = [4, 8, 16, 32, 64, 128, ]
    mes_fun = :linfract
    mes_fun = :quadfract
    mes_fun = :sqrtquadfract
    mes_fun = :sqrtlinfract
    mes_fun = :powfract
    mult_el_sizes = [0.0, 0.05, 0.1, 0.2, 0.4, 1.0]
    println("#####################################################")
    println("# test_convergence_all ")
    test_convergence_all(case_data, ns, mes_fun, mult_el_sizes)
    println("#####################################################")
    println("# try_rich_all ")
    try_rich_all(case_data[1], [4, 8, 16], mes_fun, mult_el_sizes)
    try_rich_all(case_data[1], [8, 16, 32], mes_fun, mult_el_sizes)
    try_rich_all(case_data[1], [16, 32, 64], mes_fun, mult_el_sizes)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
