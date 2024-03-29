"""
Angle bar, transient vibration.

The bar is excited by a Hann-windowed force pulse.

Reference:
Ostachowicz W, Kudela P, Krawczuk M, Zak A. Guided Waves in Structures for SHM: The Time - domain Spectral Element Method. A John Wiley & Sons, Ltd., 2011.
"""
module angle_bar_examples

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
using Gnuplot; #@gp "clear"
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection
using SparseMatricesCSR
using ThreadedSparseCSR
using VibrationGEPHelpers: gep_largest

# using InteractiveUtils
# using BenchmarkTools
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection, vtkwrite

const E = 72.7*phun("GPa");
const nu = 0.33;
const rho = 2700.0*phun("kg/m^3")
const barside = 250*phun("mm")
const barlength = 1000*phun("mm")
const barthickness = 10*phun("mm")
const ksi = 0.0
const omegad = 1000*phun("Hz")
const carrier_frequency = 75*phun("kilo*Hz")
const modulation_frequency = carrier_frequency/4
const totalforce = 1*phun("N")
const forcepatchradius = 2*phun("mm")
const forcedensity = 2*totalforce/(pi*forcepatchradius^2)
const tend = 0.3*phun("milli*s")
const visualize = true
const color = "magenta"
const distributedforce = true

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

function _execute_parallel_csr(nref = 2, nthr = 0)
        
    tolerance = barthickness/nref/10
    X = [0.0 0.0 0.0;
    0.0 barside 0.0;
    0.0 0.0 barside;
    ]
    l2fens = FENodeSet(X);
    conn = [2 1; 1 3; ]
    l2fes = FESetL2(conn);
    nLayers = 4
    ex = function (x, k)
        x[1] += k*barlength/nLayers
        x
    end
    fens,fes = Q4extrudeL2(l2fens, l2fes, nLayers, ex); # Mesh
    fens, fes = Q4toT3(fens, fes)
    for r in 1:nref
        fens, fes = T3refine(fens, fes)
    end
    @show count(fens)
    vtkwrite("angle_bar.vtu", fens, fes)

    bfes = meshboundary(fes)
    candidates = connectednodes(bfes)
    fens, fes = mergenodes(fens, fes, tolerance, candidates)
    bfes = meshboundary(fes)
    @info "Mesh $(count(fens)) nodes, $(count(fes)) elements"

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, TriRule(1)))
    C = connectionmatrix(femm, count(fens))
    perm = symrcm(C)
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), barthickness), mater)
    # Set up
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
    K = FEMMShellT3FFModule.stiffness(femm, geom0, u0, Rfield0, dchi);
    K = SparseMatricesCSR.sparsecsr(findnz(K)..., size(K)...)
    M = FEMMShellT3FFModule.mass(femm, SysmatAssemblerSparseDiag(), geom0, dchi);
    
    # Solve
    @time omega_max = pwr_largest(K, M)
    @show omega_max = max(omega_max, 20*2*pi*carrier_frequency)
    @show dt = Float64(0.9* 2/omega_max) * (sqrt(1+ksi^2) - ksi)

    U0 = gathersysvec(dchi)
    V0 = deepcopy(U0)
    
    mpoint = selectnode(fens; nearestto=[barlength/2 barside 0.0])[1]
    cpoint = selectnode(fens; nearestto=[barlength/2 0 0])[1]
    
    
    points = Dict(
    "Load" => selectnode(fens; nearestto=[barlength/2 barside 0.0])[1],
    "Crease" => selectnode(fens; nearestto=[barlength/2 0 0])[1],
    )
    
    pointdofs = Dict()
    for k in keys(points)
        pointdofs[k] = dchi.dofnums[points[k], 3]
    end
        
    function computetrac!(forceout, XYZ, tangents, feid, qpid)
        dx = XYZ[1] - fens.xyz[points["Load"], 1]
        dy = XYZ[2] - fens.xyz[points["Load"], 2]
        dz = XYZ[3] - fens.xyz[points["Load"], 3]
        d = sqrt(dx^2+dy^2+dz^2)
        forceout .= 0.0
        if d < forcepatchradius
            forceout[3] = -forcedensity
        end
        return forceout
    end

    if distributedforce
        lfemm = FEMMBase(IntegDomain(fes, TriRule(6)))
        fi = ForceIntensity(Float64, 6, computetrac!);
        Fmag = distribloads(lfemm, geom0, dchi, fi, 2);
    else
        Fmag = fill(0.0, nfreedofs(dchi))
        Fmag[pointdofs["Load"]] = -totalforce
    end

    @show sum(Fmag), -totalforce

    function force!(F, t)
        mul = 0.0
        if t <= 1/modulation_frequency
            mul = 0.5 * (1 - cos(2*pi*modulation_frequency*t)) * sin(2*pi*carrier_frequency*t)
        end
        F .= mul .* Fmag
        return F
    end
    
    nsteps = Int(round(tend/dt))
    pointdeflections = Dict()
    for k in keys(points)
        pointdeflections[k] = fill(0.0, nsteps+1)
    end
    displacements = []
    nbtw = Int(round(nsteps/100))


    peek(step, U, t) = begin
        for k in keys(points)
            b = pointdeflections[k]
            b[step+1] = U[pointdofs[k]]
        end
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
        # @gp  :- collect(0.0:dt:(nsteps*dt)) cdeflections " lw 2 lc rgb '$color' with lines title 'Deflection at the center' "  :-
        for k in sort(string.(keys(points)))
            @gp  :- collect(0.0:dt:(nsteps*dt)) pointdeflections[k] " lw 2 lc rgb '$color' with lines title 'Deflection at $(string(k))' "  :-
        end

        @gp  :- "set xlabel 'Time'" :-
        @gp  :- "set ylabel 'Deflection'" :-
        @gp  :- "set title 'Free-floating plate'"


        # Visualization
        @info "Dumping visualization"
        times = Float64[]
        vectors = []
        for i in eachindex(displacements)
            scattersysvec!(dchi, displacements[i])
            push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
            push!(times, i*dt*nbtw)
        end
        vtkwritecollection("angle_bar_$nref", fens, fes, times; vectors = vectors)
    end
end

function test_parallel_csr(nrefs = [4], nthr = 0)
    @info "Angle bar, nrefs = $nrefs: parallel CSR"
    for nref in nrefs
        _execute_parallel_csr(nref, nthr)
    end
    return true
end

function allrun(nrefs = [7], nthr = 0)
    println("#####################################################")
    println("# test_parallel_csr ")
    test_parallel_csr(nrefs, nthr)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
