"""
Spherical cap, transient vibration.

The structure is loaded with a suddenly applied pressure.
The downward deflection at the centre is monitored.
One quarter of the entire shell is modeled.
"""
module spherical_cap_expl_examples

using LinearAlgebra
using SparseArrays
using Arpack
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.AssemblyModule
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using SymRCM
using VisualStructures: default_layout_3d, plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using PlotlyJS
# using Gnuplot; @gp "clear"
using PGFPlotsX
using CSV
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection, vtkwrite
using ThreadedSparseCSR
using UnicodePlots
using InteractiveUtils
using BenchmarkTools

const E = 10.5e6*phun("psi")
const nu = 0.3;
const rho = 2.45e-4*phun("lbf*s^2/in^4")
const R = 22.27*phun("in")
const h = 0.41*phun("in")
const Rxy = (R)*sin(26.67/180*pi)
const q = 600*phun("psi")
const ksi = 0.0
const omegad = 2*pi*1000
const color = "red"
const tend = 0.001 #6e-4
const drilling_stiffness_scale = 1.0
const visualize = true

function parloop_csr!(M, K, ksi, U0, V0, tend, dt, force!, peek, nthr)
    U = deepcopy(U0)
    V = deepcopy(U0)
    A = deepcopy(U0)
    F = deepcopy(U0)
    E = deepcopy(U0)
    C = deepcopy(U0)
    C .= (ksi*2*omegad) .* vec(diag(M))
    invMC = deepcopy(U0)
    invMC .= 1.0 ./ (vec(diag(M)) .+ (dt/2) .* C)
    nsteps = Int64(round(tend/dt))
    if nsteps*dt < tend
        dt = tend / (nsteps+1)
    end

    nth = (nthr == 0 ? Base.Threads.nthreads() : nthr)
    @info "$nth threads used"
    
    t = 0.0
    # Initial Conditions
    @. U = U0; @. V = V0
    A .= invMC .* force!(F, t);
    peek(0, U, V, t)
    @time for step in 1:nsteps
        # Displacement update
        @. U += dt*V + ((dt^2)/2)*A; 
        # External loading
        force!(F, t);
        # Add elastic restoring forces
        F .-= ThreadedSparseCSR.bmul!(E, K, U)
        @inbounds @simd for i in eachindex(U)
            _Fi = F[i]; _Ai = A[i]; _Vi = V[i]
            # Add damping forces
            _Fi -= C[i] * (_Vi + (dt/2) * _Ai)
            # Compute the new acceleration.
            _A1i = invMC[i] * _Fi
            A[i] = _A1i
            # Update the velocity
            V[i] += (dt/2)* (_Ai + _A1i);
        end
        t = t + dt
        peek(step, U, V, t)
    end
end

function _execute_parallel_csr(n = 64, nthr = 0)
    color = "red"
    tolerance = min(R, h) / n  / 100
    
    tolerance = R/n/1000
    fens, fes = Q4circlen(Rxy, n)
    fens, fes = Q4toT3(fens, fes)
    
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        x=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (x, y, sqrt(R^2-x^2-y^2))
    end

    let
        fens1 = FENodeSet(fens.xyz ./ phun("in"))
        vtkwrite("spherical_cap_expl_$n.vtu", fens1, fes)
    end

    # bfes = meshboundary(fes)
    # candidates = connectednodes(bfes)
    # fens, fes = mergenodes(fens, fes, tolerance, candidates)
    bfes = meshboundary(fes)
    @info "Mesh $(count(fens)) nodes, $(count(fes)) elements"

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, TriRule(1)))
    C = connectionmatrix(femm, count(fens))
    perm = symrcm(C)
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), h), mater)
    # Set up
    femm.drilling_mass_scale = 1.0
    femm.drilling_stiffness_scale = drilling_stiffness_scale

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    lx = selectnode(fens; box = Float64[0 0 -Inf Rxy-10*tolerance -Inf Inf], inflate = tolerance)
    ly = selectnode(fens; box = Float64[-Inf Rxy-10*tolerance 0 0 -Inf Inf], inflate = tolerance)
    lo = setdiff(connectednodes(meshboundary(fes)), vcat(lx, ly))
    for i in [2, 4, 6]
        setebc!(dchi, ly, true, i)
    end
    for i in [1, 5, 6]
        setebc!(dchi, lx, true, i)
    end
    for i in 1:6
        setebc!(dchi, lo, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi, perm);
    # numberdofs!(dchi);

    # Assemble the system matrix
    FEMMShellT3FFModule.associategeometry!(femm, geom0)
    SM = FinEtoolsFlexStructures.AssemblyModule
    K = FEMMShellT3FFModule.stiffness(femm, SM.SysmatAssemblerSparseCSRSymm(0.0), geom0, u0, Rfield0, dchi);
    M = FEMMShellT3FFModule.mass(femm, SysmatAssemblerSparseDiag(), geom0, dchi);
    
    # Solve
    function pwr(K, M)
        invM = fill(0.0, size(M, 1))
        invM .= 1.0 ./ (vec(diag(M)))
        v = rand(size(M, 1))
        w = fill(0.0, size(M, 1))
        for i in 1:30
            ThreadedSparseCSR.bmul!(w, K, v)
            wn = norm(w)
            w .*= (1.0/wn)
            v .= invM .* w
            vn = norm(v)
            v .*= (1.0/vn)
        end
        sqrt((v' * (K * v)) / (v' * M * v))
    end
    @time omega_max = pwr(K, M)
    @show omega_max

    # @time evals, evecs, nconv = eigs(Symmetric(K), Symmetric(M); nev=1, which=:LM, explicittransform=:none)
    # @show sqrt(evals)
    # @show omega_max = sqrt(evals[1])

    @show dt = Float64(0.9* 2/omega_max) * (sqrt(1+ksi^2) - ksi)

    U0 = gathersysvec(dchi)
    V0 = deepcopy(U0)
    
    
    cpoint = selectnode(fens; nearestto=[0 0 R])[1]
    
    
    cpointdof = dchi.dofnums[cpoint, 3]
    
    
    function computetrac!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        r = vec(XYZ); 
        r .= vec(r)/norm(vec(r))
        forceout[1:3] .= -r*q
        forceout[4:6] .= 0.0
        return forceout
    end

    # Distributed loading on the surface of the shell
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt, 6, computetrac!);
    Fmag = distribloads(lfemm, geom0, dchi, fi, 2);

    function force!(F, t)
        F .= Fmag
        return F
    end

    nsteps = Int(round(tend/dt))
    cdeflections = fill(0.0, nsteps+1)
    kinetic_energy = fill(0.0, nsteps+1)
    displacements = []
    nbtw = Int(round(nsteps/100))


    peek(step, U, V, t) = begin
        cdeflections[step+1] = U[cpointdof]
        kinetic_energy[step+1] = 1/2 * V' * M * V
        if rem(step+1, nbtw) == 0
            push!(displacements, deepcopy(U))
        end
        nothing
    end

    @info "$nsteps steps"
    parloop_csr!(M, K, ksi, U0, V0, nsteps*dt, dt, force!, peek, nthr)

    savecsv("spherical_cap_expl_$n-KE-$(drilling_stiffness_scale)", t = collect(0.0:dt:(nsteps*dt))./phun("milli*s"), v = kinetic_energy/phun("lbf*in"))
    
    if visualize
        reference_Bathe = [   0.5310   -0.0008
         10.6195    0.0069
         16.9912    0.0167
         23.8938    0.0302
         30.2655    0.0445
         32.3894    0.0503
         35.5752    0.0513
         41.9469    0.0471
         46.7257    0.0434
         51.5044    0.0471
         55.2212    0.0548
         60.0000    0.0667
         63.1858    0.0754
         67.4336    0.0855
         71.1504    0.0902
         78.5841    0.0855
         82.8319    0.0780
         85.4867    0.0630
         89.2035    0.0460
         94.5133    0.0241
        101.4159   -0.0066
        107.2566   -0.0293
        109.9115   -0.0439
        112.5664   -0.0497
        112.5664   -0.0500
        115.2212   -0.0515
        120.0000   -0.0470
        1.20000e+02  -3.98675e-02
        1.24194e+02  -1.69095e-02
        1.28387e+02   4.63576e-03
        1.34194e+02   3.43046e-02
        1.38065e+02   5.44371e-02
        1.42903e+02   6.92715e-02
        1.46452e+02   7.98675e-02
        1.50323e+02   8.76380e-02
        1.60323e+02   9.47020e-02
        1.70323e+02   9.82340e-02
        1.75806e+02   9.36424e-02
        1.78387e+02   8.62252e-02
        1.79677e+02   7.49227e-02
        1.83548e+02   5.93819e-02
        1.88387e+02   4.94923e-02
        1.93226e+02   3.25386e-02
        1.99355e+02   8.52097e-03  ]
        reference_Bathe[:, 1] .*= 6e-4/120/phun("Milli*s")
        reference_Bathe[:, 2] .*= -1
        reference_Belytschko = [
        0   0.001290322580645
        0.000036507936508  -0.003225806451613
        0.000065079365079  -0.010967741935484
        0.000095238095238  -0.022580645161290
        0.000123809523810  -0.038387096774194
        0.000155555555556  -0.045806451612903
        0.000193650793651  -0.039354838709677
        0.000217460317460  -0.034838709677419
        0.000255555555556  -0.041935483870968
        0.000287301587302  -0.054193548387097
        0.000315873015873  -0.070000000000000
        0.000355555555556  -0.089677419354839
        0.000401587301587  -0.092903225806452
        0.000439682539683  -0.082258064516129
        0.000466666666667  -0.058064516129032
        0.000482539682540  -0.037096774193548
        0.000498412698413  -0.012903225806452
        0.000526984126984   0.017096774193548
        0.000549206349206   0.034193548387097
        0.000582539682540   0.041612903225806
        0.000598412698413   0.038387096774194]  
        reference_Belytschko[:, 1] .*= 1/phun("Milli*s") 
        reference_Tabiei = [
        0.000003183023873   0.000649350649351
        0.000049336870027  -0.005844155844156
        0.000084350132626  -0.017207792207792
        0.000084350132626  -0.017532467532468
        0.000114588859416  -0.033766233766234
        0.000144827586207  -0.046103896103896
        0.000186206896552  -0.043181818181818
        0.000206896551724  -0.041883116883117
        0.000233952254642  -0.047077922077922
        0.000267374005305  -0.062987012987013
        0.000297612732095  -0.079220779220779
        0.000332625994695  -0.089610389610390
        0.000367639257294  -0.078246753246753
        0.000401061007958  -0.055844155844156
        0.000434482758621  -0.028896103896104
        0.000458355437666  -0.005194805194805
        0.000479045092838   0.017207792207792
        0.000507692307692   0.037012987012987
        0.000539522546419   0.047402597402597
        0.000572944297082   0.035389610389610
        0.000585676392573   0.025649350649351
        0.000600000000000   0.015909090909091
        ]
        reference_Tabiei[:, 1] .*= 1/phun("Milli*s") 
        
        
        # reference[:, 2] .*= phun("in")
        # @gp  "reset" :-

        # @gp  :- "set terminal cairolatex standalone pdf size 16cm,10.5cm " :-
        # @gp  :- "set output 'fire_severity1.tex' "  :-
        # @gp "clear"
        # @gp "set key top left box  font \"Times-Roman\"" :-
        # @gp  :- reference_Bathe[:, 1] reference_Bathe[:, 2] " lw 2 lc rgb 'black' with points title 'Bathe' "  :-
        # @gp  :- reference_Belytschko[:, 1] reference_Belytschko[:, 2] " lw 2 lc rgb 'black' with points title 'Belytschko' "  :-
        # @gp  :- reference_Tabiei[:, 1] reference_Tabiei[:, 2] " lw 2 lc rgb 'black' with points title 'Tabiei' "  :-
        # @gp  :- reference_S4_impl[1]/phun("Milli*s")  reference_S4_impl[2] " lw 2 lc rgb 'black' with points title 'S4 impl' "  :-
        # @gp  :- reference_S4_expl[1]/phun("Milli*s")  reference_S4_expl[2] " lw 2 lc rgb 'black' with points title 'S4 expl' "  :-
        # @gp  :- collect(0.0:dt:(nsteps*dt))/phun("Milli*s")  cdeflections/phun("in") " lw 2 lc rgb '$color' with lines title 'Present' "  :-

        # @gp  :- "set tics font \"Times-Roman\"" :-
        # @gp  :- "set xlabel 'Time [ms]' font \"Times-Roman\"" :-
        # @gp  :- "set ylabel 'Deflection [in]' font \"Times-Roman\"" :-
        # @gp  :- "set title 'Spherical cap, deflection at the center' font \"Times-Roman\""
        # # @gp  :- "set output "  
        # save(term="cairolatex pdf input color dashed size 5in,3.3in", output="spherical_cap_expl_$n.tex")


        objects = []

        xdata, ydata = collect(0.0:dt:(nsteps*dt))/phun("Milli*s"),  cdeflections/phun("in")
        @pgf p = PGFPlotsX.Plot(
        {
        color = "black",
        line_width  = 1.0,
        },
        Coordinates([v for v in  zip(xdata, ydata)])
        )
        push!(objects, p)
        push!(objects, LegendEntry("Present"))

        xdata, ydata = reference_Bathe[:, 1], reference_Bathe[:, 2]
        @pgf p = PGFPlotsX.Plot(
        {
        only_marks,
        color = "black",
        line_width  = 1.0,
        mark = "diamond",
        },
        Coordinates([v for v in  zip(xdata, ydata)])
        )
        push!(objects, p)
        push!(objects, LegendEntry("Bathe"))

        ft = CSV.File("$(@__DIR__)/S4_impl-t.csv")  
        fw = CSV.File("$(@__DIR__)/S4_impl-w.csv")  
        xdata, ydata = ft["t"][1:47:end]/phun("Milli*s"), fw["w"][1:47:end]
        @pgf p = PGFPlotsX.Plot(
        {
        only_marks,
        color = "black",
        line_width  = 1.0,
        mark = "square",
        },
        Coordinates([v for v in  zip(xdata, ydata)])
        )
        push!(objects, p)
        push!(objects, LegendEntry("A/impl"))

        ft = CSV.File("$(@__DIR__)/S4_expl-t.csv")  
        fw = CSV.File("$(@__DIR__)/S4_expl-w.csv")  
        xdata, ydata = ft["t"][1:43:end]/phun("Milli*s"), fw["w"][1:43:end]
        @pgf p = PGFPlotsX.Plot(
        {
        only_marks,
        color = "black",
        line_width  = 1.0,
        mark = "star",
        },
        Coordinates([v for v in  zip(xdata, ydata)])
        )
        push!(objects, p)
        push!(objects, LegendEntry("A/expl"))

       
        @pgf ax = Axis(
        {
        xlabel = "Time [ms]",
        ylabel = "Apex deflection [in]",
        # ymin = -1.1e-5,
        # ymax = 1.1e-5,
        xmin = 0.0,
        xmax = 1.0,
        xmode = "linear", 
        ymode = "linear",
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
        pgfsave("spherical_cap_expl_$n.pdf", ax)


        # Visualization
        @info "Dumping visualization"
        times = Float64[]
        vectors = []
        for i in 1:length(displacements)
            scattersysvec!(dchi, displacements[i])
            push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
            push!(times, i*dt*nbtw)
        end
        vtkwritecollection("spherical_cap_expl_$n", fens, fes, times; vectors = vectors)
    end
end

function test_parallel_csr(ns = [4*64], nthr = 0)
    @info "Clamped cylinder: parallel CSR"
    for n in ns
        _execute_parallel_csr(n, nthr)
    end
    return true
end

function allrun(ns = [22])
    println("#####################################################")
    println("# test_parallel_csr ")
    test_parallel_csr(ns)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing

