"""
Aluminum plate with force at the center. One quarter is modeled.

The plate is excited with the transducer at the initial time (i.e. distributed
loading mimicking a concentrated force is applied).
"""
module plate_expl_examples

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
using Gnuplot;  @gp "clear"
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection
using ThreadedSparseCSR
using UnicodePlots
using InteractiveUtils
using BenchmarkTools
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection, vtkwrite

# E = 68*phun("GPa");
# nu = 0.33;
# rho = 2660.0*phun("kg/m^3")
# thickness = 10*phun("mm")
# D = E * thickness^3 / 12 / (1 - nu^2)
# cp = (omega) -> sqrt(omega) * (D / rho / thickness)^(1/4)
# cg = (omega) -> 2*cp(omega)
# cg(75000)

# Aluminium Alloy
const E = 68*phun("GPa");
const nu = 0.33;
const rho = 2660.0*phun("kg/m^3")
const d1 = 1000*phun("mm")
const dD = 106*phun("mm")
const dF = dD + 197*phun("mm")
const thickness = 1.0*phun("mm") # the book says 1.0
const D = E * thickness^3 / 12 / (1 - nu^2)
const cp = (omega) -> sqrt(omega) * (D / rho / thickness)^(1/4)
const cg = (omega) -> 2*cp(omega)
const ksi = 0.01
const omegad = 10000*phun("Hz")
const carrier_frequency = 35*phun("kilo*Hz")
const modulation_frequency = 7*phun("kilo*Hz")
const totalforce = 1*phun("N")
const forcepatchradius = 20*phun("mm")
const forcemagnitude = totalforce/(pi*forcepatchradius^2)
const tend = 1.0*phun("milli*s")
const visualize = true
const color = "magenta"
const distributedforce = !true

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

function _execute(nref = 2, nthr = 0)
    tolerance = d1/nref/100

    xs = sort(unique(vcat(
        linearspace(0.0, dD, 2), 
        linearspace(dD, dF, 3), 
        linearspace(dF, d1/2, 3), 
        )))
    ys = collect(xs)
    
    fens, fes = T3blockx(xs, ys); # Mesh
    for r in 1:nref
        fens, fes = T3refine(fens, fes)
    end
    
    fens.xyz = xyz3(fens)

    @info "Mesh $(count(fens)) nodes, $(count(fes)) elements"

    # vtkwrite("plate_expl.vtu", fens, fes)
    # vtkwrite("plate_expl-boundary.vtu", fens, bfes)

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

    # Symmetry EBC's
    lx = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    ly = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [2, 4, 6]
        setebc!(dchi, ly, true, i)
    end
    for i in [1, 5, 6]
        setebc!(dchi, lx, true, i)
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
    @show omega_max = max(omega_max, 10*2*pi*carrier_frequency)
    @show dt = Float64(0.9* 2/omega_max) * (sqrt(1+ksi^2) - ksi)

    U0 = gathersysvec(dchi)
    V0 = deepcopy(U0)
    
    points = Dict(
    "A" => selectnode(fens; nearestto=[0.0 0.0 0.0])[1],
    "B" => selectnode(fens; nearestto=[dD 0.0 0.0])[1],
    "C" => selectnode(fens; nearestto=[dF 0.0 0.0])[1],
    "D" => selectnode(fens; nearestto=[dD dD 0.0])[1],
    "E" => selectnode(fens; nearestto=[dD dF 0.0])[1],
    "F" => selectnode(fens; nearestto=[dF dF 0.0])[1],
    )
    
    pointdofs = Dict()
    for k in keys(points)
        pointdofs[k] = dchi.dofnums[points[k], 3]
    end
    
    # Four cycles of the carrier frequency

    function computetrac!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        dx = XYZ[1] - fens.xyz[points["A"], 1]
        dy = XYZ[2] - fens.xyz[points["A"], 2]
        dz = XYZ[3] - fens.xyz[points["A"], 3]
        d = sqrt(dx^2+dy^2+dz^2)
        forceout .= 0.0
        if d < forcepatchradius
            forceout[3] = forcemagnitude
        end
        return forceout
    end

    if distributedforce
        lfemm = FEMMBase(IntegDomain(fes, TriRule(6)))
        fi = ForceIntensity(FFlt, 6, computetrac!);
        Fmag = distribloads(lfemm, geom0, dchi, fi, 2);
    else
        Fmag = fill(0.0, dchi.nfreedofs)
        Fmag[pointdofs["A"]] = -totalforce/4 # only a quarter of the plate is modeled
    end

    @show sum(Fmag), -totalforce/4

    function force!(F, t)
        mul = 0.0
        if t <= 4/carrier_frequency
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
        
        # for k in sort(string.(keys(points)))
        #     @gp  :- collect(0.0:dt:(nsteps*dt)) pointdeflections[k] " lw 2 lc rgb '$color' with lines title 'Deflection at $(string(k))' "  :-
        # end
        
        for (i, k) in enumerate(["C"])
            @gp  :- collect(0.0:dt:(nsteps*dt))./phun("milli*s") pointdeflections[k] " lw 2 lc rgb '$color' pt $i with lp title 'Deflection at $(string(k))' "  :-
        end
        
        @gp  :- "set xlabel 'Time [ms]'" :-
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
        vtkwritecollection("plate_expl_$nref", fens, fes, times; vectors = vectors)
    end
end

function test(nrefs = [5], nthr = 0)
    @info "Cracked Plate, nrefs = $nrefs: parallel CSR"
    for nref in nrefs
        _execute(nref, nthr)
    end
    return true
end

function allrun(nrefs = [6], nthr = 0)
    println("#####################################################")
    println("# test ")
    test(nrefs, nthr)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
