"""
Clamped cylinder, transient vibration.

The cylinder is banged at the initial time (i.e. it is given an initial velocity).
"""
module clamp_cyl_transient_examples

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
using ThreadedSparseCSR
using UnicodePlots
using InteractiveUtils
using BenchmarkTools

const E = 200e9;
const nu = 0.3;
const rho = 7800.0
const R = 0.1
const L = 0.8
const ksi = 0.0
const omegad = 2*pi*1000
const color = "red"
const omegaf = 2*pi*10^5
const tend = 25 * (2*pi/omegaf)
const visualize = true
const visualizeclear = true
const visualizevtk = !true
const distributedforce = !true

# Solve
function pwr(K, M, maxit = 30)
    invM = fill(0.0, size(M, 1))
    invM .= 1.0 ./ (vec(diag(M)))
    v = rand(size(M, 1))
    w = fill(0.0, size(M, 1))
    for i in 1:maxit
        ThreadedSparseCSR.bmul!(w, K, v)
        wn = norm(w)
        w .*= (1.0/wn)
        v .= invM .* w
        vn = norm(v)
        v .*= (1.0/vn)
    end
    sqrt((v' * (K * v)) / (v' * M * v))
end

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

cylindrical!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) = begin
    r = vec(XYZ); r[2] = 0.0;
    csmatout[:, 3] .= vec(r)/norm(vec(r))
    csmatout[:, 2] .= (0.0, 1.0, 0.0) #  this is along the axis
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

function _set_up_model(n = 64, thickness = 0.01, drilling_stiffness_scale = 1.0)
    tolerance = min(R, L) / n  / 100
    
    tolerance = R/n/1000
    fens, fes = T3block(360.0, L, n, 4*n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a/180*pi), y-L/2, R*cos(a/180*pi))
    end

    bfes = meshboundary(fes)
    candidates = connectednodes(bfes)
    fens, fes = mergenodes(fens, fes, tolerance, candidates)
    bfes = meshboundary(fes)


    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, TriRule(1)))
    C = connectionmatrix(femm, count(fens))
    perm = symrcm(C)
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
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
    l1 = connectednodes(meshboundary(fes))
    for i in 1:6
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi, perm);
    # numberdofs!(dchi);

    # Assemble the system matrix
    FEMMShellT3FFModule.associategeometry!(femm, geom0)
    SM = FinEtoolsFlexStructures.AssemblyModule
    K = FEMMShellT3FFModule.stiffness(femm, SM.SysmatAssemblerSparseCSRSymm(0.0), geom0, u0, Rfield0, dchi);
    M = FEMMShellT3FFModule.mass(femm, SysmatAssemblerSparseDiag(), geom0, dchi);
    
    return fens, fes, geom0, dchi, K, M
end

function _execute_parallel_csr(n = 64, thickness = 0.01, nthr = 0)
    fens, fes, geom0, dchi, K, M = _set_up_model(n, thickness)
    
    @info "Mesh $(count(fens)) nodes, $(count(fes)) elements"

    # Compute the maximum eigenvalue to determine the stable time step
    omega_max = pwr(K, M)
    

    # @time evals, evecs, nconv = eigs(Symmetric(K), Symmetric(M); nev=1, which=:LM, explicittransform=:none)
    # @show sqrt(evals)
    # @show omega_max = sqrt(evals[1])

    dt = Float64(0.9* 2/omega_max) * (sqrt(1+ksi^2) - ksi)

    U0 = gathersysvec(dchi)
    V0 = deepcopy(U0)
    
    qpoint = selectnode(fens; nearestto=[0 -L/4 R])[1]
    cpoint = selectnode(fens; nearestto=[0 0 R])[1]
    
    qpointdof = dchi.dofnums[qpoint, 3]
    cpointdof = dchi.dofnums[cpoint, 3]
    cpointdof6 = dchi.dofnums[cpoint, 6]
    
    function computetrac!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        dx = XYZ[1] - fens.xyz[qpoint, 1]
        dy = XYZ[2] - fens.xyz[qpoint, 2]
        dz = XYZ[3] - fens.xyz[qpoint, 3]
        mag = 10^5*exp(-20/R*sqrt(dx^2+dy^2+dz^2))
        r = vec(XYZ); r[2] = 0.0
        r .= vec(r)/norm(vec(r))
        theta = atan(r[3], r[1])
        forceout[1:3] .= r*mag
        forceout[4:6] .= 0.0
        return forceout
    end

    # Sinusoidal loading on the surface of the shell
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt, 6, computetrac!);
    Fmag = distribloads(lfemm, geom0, dchi, fi, 2);

    function force!(F, t)
        if t < 2 * (2*pi/omegaf)
            F .= sin(omegaf*t) .* Fmag
        else
            F .= 0.0
        end
        # @show norm(F)
        return F
    end

    # Suddenly applied constant force at a node
    # Fmag = fill(0.0, dchi.nfreedofs)
    # Fmag[qpointdof] = -1.0
    # function force!(F, t)
    #     F .= Fmag
    #     return F
    # end
    
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
        if visualizevtk
            @info "Dumping visualization"
            times = Float64[]
            vectors = []
            for i in 1:length(displacements)
                scattersysvec!(dchi, displacements[i])
                push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
                push!(times, i*dt*nbtw)
            end
            vtkwritecollection("clamp_cyl_forc_damp_expl_$n", fens, fes, times; vectors = vectors)
        end
    end
end

function _test_stable_time_step(n = 64, thickness = 0.01, drilling_stiffness_scale = 1.0)
    fens, fes, geom0, dchi, K, M = _set_up_model(n, thickness, drilling_stiffness_scale)
    
    omega_max = pwr(K, M, 300)
    

    # evals, evecs, nconv = eigs(Symmetric(K), Symmetric(M); nev=1, which=:LM, explicittransform=:none)
    
    # omega_max = sqrt(evals[1])

    dt = Float64(0.9* 2/omega_max) * (sqrt(1+ksi^2) - ksi)
    return dt
end

function test_parallel_csr(ns = [4*64], nthr = 0)
    @info "Clamped cylinder: transient vibration"
    for n in ns
        _execute_parallel_csr(n, 0.01, nthr)
    end
    return true
end


function test_stable_time_step()
    @info "Clamped cylinder: test stable time step"
    results = []
    for t in [0.01, 0.002, 0.001, 0.0001, 0.00001]
        dts = []
        for d in [0.01, 0.1, 1.0, 1.5, 2.0, 2.5, 5.0, 10.0]
            dt = _test_stable_time_step(2*64, t, d)
            push!(dts, dt)
        end
        push!(results, (t, dts))
    end
    return results
end


function allrun(ns = [2*64])
     println("#####################################################")
    println("# test_parallel_csr ")
    test_parallel_csr(ns)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
