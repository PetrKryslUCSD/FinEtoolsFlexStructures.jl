"""


The cylinder is banged at the initial time (i.e. it is given an initial velocity).
"""
module plate_with_crack_examples

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
using Gnuplot; # @gp "clear"
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection
using ThreadedSparseCSR
using UnicodePlots
using InteractiveUtils
using BenchmarkTools
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection, vtkwrite

const E = 72.7*phun("GPa");
const nu = 0.33;
const rho = 2700.0*phun("kg/m^3")
const d1 = 500*phun("mm")
const d2 = 50*phun("mm")
const d3 = 750*phun("mm")
const d4 = 250*phun("mm")
const thickness = 5*phun("mm")
const ksi = 0.0
const omegad = 1000*phun("Hz")
const carrier_frequency = 75*phun("kilo*Hz")
const modulation_frequency = carrier_frequency/4
const forcepatchradius = 20*phun("mm")
const forcemagnitude = 1*phun("N")
const color = "red"
const tend = 0.3*phun("milli*s")
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
    peek(0, U, t)
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
        peek(step, U, t)
    end
end

function _execute_parallel_csr(nref = 2, nthr = 0, color = "red")
    tolerance = d2/nref/100

    d5 = (d1 - d2) / 2
    d6 = (d3 - d4) 
    nd5 = 5
    xs = sort(unique(vcat(
        linearspace(0.0, d5, nd5), 
        linearspace(d5, d5+d2, 2), 
        linearspace(d5+d2, d1, nd5)
        )))
    ys = collect(linearspace(0.0, d4, nd5))
    
    fens1, fes1 = T3blockx(xs, ys); # Mesh
    for r in 1:nref
        fens1, fes1 = T3refine(fens1, fes1)
    end
    @show count(fens1), count(fes1)

    ys = collect(linearspace(d4, d3, nd5*2))
    
    fens2, fes2 = T3blockx(xs, ys); # Mesh
    for r in 1:nref
        fens2, fes2 = T3refine(fens2, fes2)
    end
    @show count(fens2), count(fes2)

    fens, fes1, fes2 = mergemeshes(fens1, fes1,  fens2, fes2, 0.0)
    fes = cat(fes1, fes2)
    @show count(fens), count(fes)

    offset = max(d5 / 4 / nref / 10, 2*tolerance)
    l1 = selectnode(fens; box = [0 d5-offset d4 d4], inflate = tolerance)
    l2 = selectnode(fens; box = [d5+d2+offset d1 d4 d4], inflate = tolerance)


    
    candidates = vcat(l1, l2)
    fens, fes = mergenodes(fens, fes, tolerance, candidates)
    bfes = meshboundary(fes)
    @info "Mesh $(count(fens)) nodes, $(count(fes)) elements"

    @show count(fens)
    vtkwrite("plate_with_crack.vtu", fens, fes)
    vtkwrite("plate_with_crack-boundary.vtu", fens, bfes)

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, TriRule(1)))
    C = connectionmatrix(femm, count(fens))
    perm = symrcm(C)
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), thickness), mater)
    # Set up
    femm.drilling_mass_scale = 1.0
    femm.drilling_stiffness_scale = 1.0

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # No EBC's
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
    @show omega_max = max(omega_max, 20*2*pi*carrier_frequency)
    @show dt = Float64(0.9* 2/omega_max) * (sqrt(1+ksi^2) - ksi)

    U0 = gathersysvec(dchi)
    V0 = deepcopy(U0)
    
    mpoint = selectnode(fens; nearestto=[barlength/2 barside 0.0])[1]
    cpoint = selectnode(fens; nearestto=[barlength/2 0 0])[1]
    
    mpointdof = dchi.dofnums[mpoint, 3]
    cpointdof = dchi.dofnums[cpoint, 3]
    
    # Four cycles of the carrier frequency

    function computetrac!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        dx = XYZ[1] - fens.xyz[mpoint, 1]
        dy = XYZ[2] - fens.xyz[mpoint, 2]
        dz = XYZ[3] - fens.xyz[mpoint, 3]
        forceout[1:2] .= 0.0
        forceout[3] = forcemagnitude*exp(-20*sqrt(dx^2+dy^2+dz^2)/forcepatchradius)
        forceout[4:6] .= 0.0
        return forceout
    end

    # Sinusoidal loading on the surface of the shell
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt, 6, computetrac!);
    Fmag = distribloads(lfemm, geom0, dchi, fi, 2);

    function force!(F, t)
        mul = 0.0
        if t <= 4/carrier_frequency
            mul = 0.5 * (1 - cos(2*pi*modulation_frequency*t)) * sin(2*pi*carrier_frequency*t)
        end
        F .= mul .* Fmag
        return F
    end
    
    nsteps = Int(round(tend/dt))
    cdeflections = fill(0.0, nsteps+1)
    displacements = []
    nbtw = Int(round(nsteps/100))


    peek(step, U, t) = begin
        cdeflections[step+1] = U[cpointdof]
        if rem(step+1, nbtw) == 0
            push!(displacements, deepcopy(U))
        end
        nothing
    end

    @info "$nsteps steps"
    parloop_csr!(M, K, ksi, U0, V0, nsteps*dt, dt, force!, peek, nthr)
    
    if visualize
        # @gp  "set terminal windows 0 "  :-
        # @gp "clear"
        @gp  :- collect(0.0:dt:(nsteps*dt)) cdeflections " lw 2 lc rgb '$color' with lines title 'Deflection at the center' "  :-

        @gp  :- "set xlabel 'Time'" :-
        @gp  :- "set ylabel 'Deflection'" :-
        @gp  :- "set title 'Free-floating plate'"


        # Visualization
        @info "Dumping visualization"
        times = Float64[]
        vectors = []
        for i in 1:length(displacements)
            scattersysvec!(dchi, displacements[i])
            push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
            push!(times, i*dt*nbtw)
        end
        vtkwritecollection("plate_with_crack_$nref", fens, fes, times; vectors = vectors)
    end
end

function test_parallel_csr(nrefs = [4], nthr = 0, color = "red")
    @info "Cracked Plate, nrefs = $nrefs: parallel CSR"
    for nref in nrefs
        _execute_parallel_csr(nref, nthr, color)
    end
    return true
end

function allrun(nrefs = [2], nthr = 0)
    println("#####################################################")
    println("# test_parallel_csr ")
    test_parallel_csr(nrefs, nthr, "blue")
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
