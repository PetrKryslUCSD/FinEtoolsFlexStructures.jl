"""
Half pipe of aluminum, with an open through crack.

The pipe corner is excited by a concentrated force in the form of a
Hann-windowed sinusoidal pulse.

The crack configuration is not exactly the same one as in the book.
It is not described unambiguously in the book.

Reference:
Ostachowicz W, Kudela P, Krawczuk M, Zak A. Guided Waves in Structures for SHM: The Time - domain Spectral Element Method. A John Wiley & Sons, Ltd., 2011.
"""
module cracked_half_pipe_examples

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
using Gnuplot;  # @gp "clear"
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection
using ThreadedSparseCSR
using GEPHelpers: pwr_largest
# using InteractiveUtils
# using BenchmarkTools
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection, vtkwrite


const E = 72.7*phun("GPa");
const nu = 0.33;
const rho = 2700.0*phun("kg/m^3")
const d1 = pi*phun("rad")
const d2 = 6.0/180*pi*phun("rad")
const d5 = (180-40.0)/180*pi*phun("rad") - d2/2
const d6 = d1 - d5 - d2
const nd5 = 19
const nd2 = 2
const nd6 = 2
const nd3 = 25
const nd4 = 25
const halfcirc = 500*phun("mm")
const radius = halfcirc / pi
const d3 = 1000*phun("mm")
const d4 = d3/2
const thickness = 10*phun("mm")
const ksi = 0.0
const omegad = 1000*phun("Hz")
const carrier_frequency = 75*phun("kilo*Hz")
const modulation_frequency = carrier_frequency/4
const totalforce = 1*phun("N")
const forcepatchradius = 3*phun("mm")
const forcedensity = totalforce/(pi*forcepatchradius^2)
const distributedforce = !true
const tend = 0.5*phun("milli*s")

const visualize = true
const visualizeclear = true
const visualizevtk = !true
const color = "black"

# function _pwr(K, M, maxit = 30, rtol = 1/10000)
#     invM = fill(0.0, size(M, 1))
#     invM .= 1.0 ./ (vec(diag(M)))
#     v = rand(size(M, 1))
#     w = fill(0.0, size(M, 1))
#     everyn = Int(round(maxit / 50)) + 1
#     lambda = lambdap = 0.0
#     for i in 1:maxit
#         ThreadedSparseCSR.bmul!(w, K, v)
#         wn = norm(w)
#         w .*= (1.0/wn)
#         v .= invM .* w
#         vn = norm(v)
#         v .*= (1.0/vn)
#         if i % everyn  == 0
#             lambda = sqrt((v' * (K * v)) / (v' * M * v))
#             # @show i, abs(lambda - lambdap) / lambda
#             if abs(lambda - lambdap) / lambda  < rtol
#                 break
#             end
#             lambdap = lambda
#         end
#     end
#     return lambda
# end

function _cd_loop!(M, K, ksi, U0, V0, tend, dt, force!, peek)
    # Central difference integration loop, for mass-proportional Rayleigh damping.
    U = similar(U0)
    V = similar(U0)
    A = similar(U0)
    F = similar(U0)
    E = similar(U0)
    C = similar(U0)
    C .= (ksi*2*omegad) .* vec(diag(M))
    invMC = similar(U0)
    invMC .= 1.0 ./ (vec(diag(M)) .+ (dt/2) .* C)
    nsteps = Int64(round(tend/dt))
    if nsteps*dt < tend
        dt = tend / (nsteps+1)
    end
    dt2_2 = ((dt^2)/2)
    dt_2 = (dt/2)

    t = 0.0
    @. U = U0; @. V = V0; # Initial Conditions
    A .= invMC .* force!(F, t); # Compute initial acceleration
    peek(0, U, V, t)
    for step in 1:nsteps
        t = t + dt
        @. U += dt*V + dt2_2*A;   # Displacement update
        force!(F, t);  # External forces are computed
        ThreadedSparseCSR.bmul!(E, K, U)  # Evaluate elastic restoring forces
        @. F -= E + C * (V + dt_2 * A)  # Calculate total force
        @. V += dt_2 * A; # Update the velocity, part one: add old acceleration
        @. A = invMC * F # Compute the acceleration
        @. V += dt_2 * A; # Update the velocity, part two: add new acceleration
        peek(step, U, V, t)
    end
end

function _execute(nref = 2, nthr = 0, color = "red")
    tolerance = min(d2/nd2, d3/nd3, d4/nd4, d5/nd5)/nref/10

    xs = sort(unique(vcat(
        linearspace(0.0, d5, nd5), 
        linearspace(d5, d5+d2, nd2), 
        linearspace(d5+d2, d1, nd6)
        )))
    ys = collect(linearspace(0.0, d4, nd4))
    
    fens1, fes1 = T3blockx(xs, ys); # Mesh
    for r in 1:nref
        fens1, fes1 = T3refine(fens1, fes1)
    end
    
    ys = collect(linearspace(d4, d3, nd3))
    
    fens2, fes2 = T3blockx(xs, ys); # Mesh
    for r in 1:nref
        fens2, fes2 = T3refine(fens2, fes2)
    end
    
    fens, fes1, fes2 = mergemeshes(fens1, fes1,  fens2, fes2, 0.0)
    fes = cat(fes1, fes2)
    
    offset = min(d2 / nd2 / nref / 10, 4*tolerance)
    l1 = selectnode(fens; box = [0 d5+offset d4 d4], inflate = tolerance)
    l2 = selectnode(fens; box = [d5+d2-offset d1 d4 d4], inflate = tolerance)
    la = selectnode(fens; box = [-Inf Inf d4 d4], inflate = tolerance)

    
    candidates = vcat(l1, l2)
    fens, fes = mergenodes(fens, fes, tolerance, candidates)
    bfes = meshboundary(fes)
    @info "Mesh $(count(fens)) nodes, $(count(fes)) elements"

    fens.xyz = xyz3(fens)

    for i in 1:count(fens)
        a, y = fens.xyz[i, 1:2]
        fens.xyz[i, :] .= cos(pi-a)*radius, y, sin(pi-a)*radius
    end

    vtkwrite("cracked_half_pipe.vtu", fens, fes)
    vtkwrite("cracked_half_pipe-boundary.vtu", fens, bfes)

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, TriRule(1)))
    C = connectionmatrix(femm, count(fens))
    perm = symrcm(C)
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), thickness), mater)
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
    K = FEMMShellT3FFModule.stiffness(femm, SM.SysmatAssemblerSparseCSRSymm(0.0), geom0, u0, Rfield0, dchi);
    M = FEMMShellT3FFModule.mass(femm, SysmatAssemblerSparseDiag(), geom0, dchi);
    
    omega_max = pwr_largest(K, M)
    omega_max = max(omega_max, 20*2*pi*carrier_frequency)
    @show dt = Float64(0.99 * 2/omega_max)

    U0 = gathersysvec(dchi)
    V0 = deepcopy(U0)
    
    mpoint = selectnode(fens; nearestto=[radius 0 0.0])[1]
    cpoint = selectnode(fens; nearestto=[radius d3 0])[1]
    
    mpointdof1 = dchi.dofnums[mpoint, 1]
    mpointdof3 = dchi.dofnums[mpoint, 3]
    cpointdof = dchi.dofnums[cpoint, 3]
    
    # Four cycles of the carrier frequency

    function computetrac!(forceout, XYZ, tangents, feid, qpid)
        dx = XYZ[1] - fens.xyz[mpoint, 1]
        dy = XYZ[2] - fens.xyz[mpoint, 2]
        dz = XYZ[3] - fens.xyz[mpoint, 3]
        d = (dx^2+dy^2+dz^2)/forcepatchradius
        forceout[1:3] .= 0.0
        if d < 1.0
            forceout[1] = forcedensity*exp(-20*d^2/forcepatchradius^2)
        end
        forceout[4:6] .= 0.0
        return forceout
    end

    if distributedforce
        lfemm = FEMMBase(IntegDomain(fes, TriRule(6)))
        fi = ForceIntensity(Float64, 6, computetrac!);
        Fmag = distribloads(lfemm, geom0, dchi, fi, 2);
    else
        Fmag = fill(0.0, nfreedofs(dchi))
        Fmag[mpointdof1] = -totalforce/4 # only a quarter of the plate is modeled
    end

    function force!(F, t)
        mul = 0.0
        if t <= 1/modulation_frequency
            mul = 0.5 * (1 - cos(2*pi*modulation_frequency*t)) * sin(2*pi*carrier_frequency*t)
        end
        F .= mul .* Fmag
        return F
    end

    nsteps = Int(round(tend/dt))
    cdeflections = fill(0.0, nsteps+1)
    mdeflections = fill(0.0, nsteps+1)
    displacements = []
    nbtw = Int(round(nsteps/100))


    peek(step, U, V, t) = begin
        cdeflections[step+1] = U[cpointdof]
        mdeflections[step+1] = U[mpointdof1]
        if rem(step+1, nbtw) == 0
            push!(displacements, deepcopy(U))
        end
        nothing
    end

    @info "$nsteps steps"
    @time _cd_loop!(M, K, ksi, U0, V0, nsteps*dt, dt, force!, peek)
        
    if visualize
        # @gp  "set terminal windows 0 "  :-
        # @gp "clear"
        @gp  :- collect(0.0:dt:(nsteps*dt)) cdeflections " lw 2 lc rgb '$color' with lines title 'Deflection at the center' "  :-
        @gp  :- collect(0.0:dt:(nsteps*dt)) mdeflections " lw 4 lc rgb '$color' with lines title 'Deflection at the source' "  :-
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
        vtkwritecollection("cracked_half_pipe_$nref", fens, fes, times; vectors = vectors)
    end
end

function test(nrefs = [4], nthr = 0, color = "red")
    @info "Cracked half pipe, nrefs = $nrefs: parallel CSR"
    for nref in nrefs
        _execute(nref, nthr, color)
    end
    return true
end

function allrun(nrefs = [4], nthr = 0, color = "red")
    println("#####################################################")
    println("# test ")
    test(nrefs, nthr, color)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
