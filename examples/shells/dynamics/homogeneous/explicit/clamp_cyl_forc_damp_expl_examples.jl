"""
Clamped cylinder, transient vibration.

The cylinder is banged at the initial time (i.e. it is given an initial velocity).
"""
module clamp_cyl_forc_damp_expl_examples

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
using SparseMatricesCSR
using GEPHelpers: pwr_largest

# using InteractiveUtils
# using BenchmarkTools

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


function loop!(M, K, ksi, U0, V0, tend, dt, force!, peek)
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
        F .-=  mul!(E, K, U)
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

struct ThreadBuffer{KVT}
    kcolumns::KVT
    uv::SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true} 
    result::Vector{Float64}
end

struct ThreadBufferVOM{RT}
    rang::RT
    result::Vector{Float64}
end

# function parmul!(R, tbs)
#     tasks = [];
#     for th in 1:length(tbs)
# tb = tbs[th]
# push!(tasks, Threads.@spawn begin 
#     mul!(tb.result, tb.kcolumns, tb.uv)
# end);
#     end
#     Threads.wait(tasks[1]);
#     R .= tbs[1].result
#     for th in 2:length(tbs)
# Threads.wait(tasks[th]);
# R .+= tbs[th].result
#     end
#     R
# end

function parmul!(R, tbs)
    Threads.@threads for i in eachindex(tbs)
        if i == 1
            mul!(R, tbs[i].kcolumns, tbs[i].uv)
        else
            mul!(tbs[i].result, tbs[i].kcolumns, tbs[i].uv)
        end
    end
    for i in 2:length(tbs)
        R .+= tbs[i].result
    end
    R
end

function parloop!(M, K, ksi, U0, V0, tend, dt, force!, peek, nthr)
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
    chunk = Int(floor(length(U) / nth))
    threadbuffs = ThreadBuffer[];
    for th in 1:nth
        colrange = th < nth ? (chunk*(th-1)+1:chunk*(th)+1-1) : (chunk*(th-1)+1:length(U))
        push!(threadbuffs, ThreadBuffer(K[:, colrange], view(U, colrange), deepcopy(U)));
    end

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
        F .-=  parmul!(E, threadbuffs)
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

function parvommul!(R, vom, P, threadbuffs, U)
    Threads.@threads for i in eachindex(threadbuffs)
        if i == 1
            vommul!(R, vom, P, threadbuffs[i].rang, U)
        else
            vommul!(threadbuffs[i].result, vom, P, threadbuffs[i].rang, U)
        end
    end
    for i in 2:length(threadbuffs)
        R .+= threadbuffs[i].result
    end
    R
end

function vommul!(R, vom, P, rang, U)
    z = zero(eltype(R))
    R .= z
    # v = Kasm.buffer
    # P = Kasm.elem_mat_dim
    # r = fill(z, P)
    # u = fill(z, P)
    for i in rang
        d = vom[i][1]
        m = vom[i][2]
        @inbounds for i in 1:P
            acc = z
            @inbounds for j in 1:P
                acc += m[i, j] * U[d[j]]
            end
            R[d[i]] += acc
        end
    end
    R
end

function parloop_vom!(M, Kasm, threadbuffs, ksi, U0, V0, tend, dt, force!, peek, nthr)
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

    vom = Kasm.buffer
    P = Kasm.elem_mat_dim
    
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
        F .-= parvommul!(E, vom, P, threadbuffs, U)
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

cylindrical!(csmatout, XYZ, tangents, feid, qpid) = begin
    r = vec(XYZ); r[2] = 0.0;
    csmatout[:, 3] .= vec(r)/norm(vec(r))
    csmatout[:, 2] .= (0.0, 1.0, 0.0) #  this is along the axis
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end



function _execute_parallel(n = 64, thickness = 0.01, nthr = 0)
    color = "red"
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
    @info "Mesh $(count(fens)) nodes, $(count(fes)) elements"

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, TriRule(1)))
    C = connectionmatrix(femm, count(fens))
    perm = symrcm(C)
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    
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
    K = FEMMShellT3FFModule.stiffness(femm, SysmatAssemblerFFBlock(nfreedofs(dchi)), geom0, u0, Rfield0, dchi);
    M = FEMMShellT3FFModule.mass(femm, SysmatAssemblerFFBlock(SysmatAssemblerSparseDiag(), nfreedofs(dchi), nfreedofs(dchi)), geom0, dchi);
    I, J, V = findnz(K)
    # Use the CSR package.
    K = SparseMatricesCSR.sparsecsr(I, J, V, size(K)...)
    @show typeof(K), typeof(M)
    
    # Solve
    @time omega_max = pwr_largest(K, M)
    @show omega_max

    # @time evals, evecs, nconv = eigs(Symmetric(K), Symmetric(M); nev=1, which=:LM, explicittransform=:none)
    # @show sqrt(evals)
    # @show omega_max = sqrt(evals[1])

    @show dt = Float64(0.9* 2/omega_max) * (sqrt(1+ksi^2) - ksi)

    U0 = gathersysvec(dchi)
    V0 = deepcopy(U0)
    
    qpoint = selectnode(fens; nearestto=[0 -L/4 R])[1]
    cpoint = selectnode(fens; nearestto=[0 0 R])[1]
    
    qpointdof = dchi.dofnums[qpoint, 3]
    cpointdof = dchi.dofnums[cpoint, 3]
    cpointdof6 = dchi.dofnums[cpoint, 6]
    
    function computetrac!(forceout, XYZ, tangents, feid, qpid)
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
    fi = ForceIntensity(Float64, 6, computetrac!);
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
    # Fmag = fill(0.0, nfreedofs(dchi))
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
    parloop!(M, K, ksi, U0, V0, nsteps*dt, dt, force!, peek, nthr)
    
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
        for i in eachindex(displacements)
            scattersysvec!(dchi, displacements[i])
            push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
            push!(times, i*dt*nbtw)
        end
        vtkwritecollection("clamp_cyl_forc_damp_expl_$n", fens, fes, times; vectors = vectors)
    end
end

function _execute_parallel_csr(n = 64, thickness = 0.01, nthr = 0)
    color = "green"
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
    @info "Mesh $(count(fens)) nodes, $(count(fes)) elements"

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, TriRule(1)))
    C = connectionmatrix(femm, count(fens))
    perm = symrcm(C)
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    

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
    # K = FEMMShellT3FFModule.stiffness(femm, SM.SysmatAssemblerSparseCSRSymm(0.0), geom0, u0, Rfield0, dchi);
    # M = FEMMShellT3FFModule.mass(femm, SysmatAssemblerSparseDiag(), geom0, dchi);
    K = FEMMShellT3FFModule.stiffness(femm, SysmatAssemblerFFBlock(nfreedofs(dchi)), geom0, u0, Rfield0, dchi);
    M = FEMMShellT3FFModule.mass(femm, SysmatAssemblerFFBlock(SysmatAssemblerSparseDiag(), nfreedofs(dchi), nfreedofs(dchi)), geom0, dchi);
    I, J, V = findnz(K)
    # Use the CSR package.
    K = SparseMatricesCSR.sparsecsr(I, J, V, size(K)...)
    @show typeof(K), typeof(M)
    
    @time omega_max = pwr_largest(K, M)
    @show omega_max

    # @time evals, evecs, nconv = eigs(Symmetric(K), Symmetric(M); nev=1, which=:LM, explicittransform=:none)
    # @show sqrt(evals)
    # @show omega_max = sqrt(evals[1])

    @show dt = Float64(0.9* 2/omega_max) * (sqrt(1+ksi^2) - ksi)

    U0 = gathersysvec(dchi)
    V0 = deepcopy(U0)
    
    qpoint = selectnode(fens; nearestto=[0 -L/4 R])[1]
    cpoint = selectnode(fens; nearestto=[0 0 R])[1]
    
    qpointdof = dchi.dofnums[qpoint, 3]
    cpointdof = dchi.dofnums[cpoint, 3]
    cpointdof6 = dchi.dofnums[cpoint, 6]
    
    function computetrac!(forceout, XYZ, tangents, feid, qpid)
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

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Sinusoidal loading on the surface of the shell
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64, 6, computetrac!);
    Fmag = distribloads(lfemm, vassem, geom0, dchi, fi, 2);

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
    # Fmag = fill(0.0, nfreedofs(dchi))
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
        @info "Dumping visualization"
        times = Float64[]
        vectors = []
        for i in eachindex(displacements)
            scattersysvec!(dchi, displacements[i])
            push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
            push!(times, i*dt*nbtw)
        end
        vtkwritecollection("clamp_cyl_forc_damp_expl_$n", fens, fes, times; vectors = vectors)
    end
end

  function _execute_serial(n = 64, thickness = 0.01)
    color = "blue"  

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
    @info "Mesh $(count(fens)) nodes, $(count(fes)) elements"

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, TriRule(1)))
    C = connectionmatrix(femm, count(fens))
    perm = symrcm(C)
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    

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
    K = FEMMShellT3FFModule.stiffness(femm, SysmatAssemblerFFBlock(nfreedofs(dchi)), geom0, u0, Rfield0, dchi);
    M = FEMMShellT3FFModule.mass(femm, SysmatAssemblerFFBlock(SysmatAssemblerSparseDiag(), nfreedofs(dchi), nfreedofs(dchi)), geom0, dchi);
    
    I, J, V = findnz(K)
    # Use the CSR package.
    K = SparseMatricesCSR.sparsecsr(I, J, V, size(K)...)

    # Solve
    # function pwr(K, M)
    #     invM = fill(0.0, size(M, 1))
    #     invM .= 1.0 ./ (vec(diag(M)))
    #     v = rand(size(M, 1))
    #     w = fill(0.0, size(M, 1))
    #     for i in 1:40
    #         mul!(w, K, v)
    #         wn = norm(w)
    #         w .*= (1.0/wn)
    #         v .= invM .* w
    #         vn = norm(v)
    #         v .*= (1.0/vn)
    #     end
    #     # ev = sqrt((v' * K * v) / (v' * M * v)) # EXTREMELY SLOW
    #     mul!(w, K, v) # FAST
    #     ev = sqrt((w' * v) / (v' * (M * v))) # FAST
    #     return ev
    # end
    @info "Power iteration"
    @time omega_max = pwr_largest(K, M)
    @show omega_max
    @show dt = Float64(0.9* 2/omega_max) * (sqrt(1+ksi^2) - ksi)

    # @time evals, evecs, nconv = eigs(Symmetric(K), Symmetric(M); nev=1, which=:LM, explicittransform=:none)
    # # @show sqrt(evals)
    # @show omega_max = sqrt(evals[1])
    # @show dt = Float64(0.9* 2/omega_max) * (sqrt(1+ksi^2) - ksi)

    U0 = gathersysvec(dchi)
    V0 = deepcopy(U0)
    
    qpoint = selectnode(fens; nearestto=[0 -L/4 R])[1]
    cpoint = selectnode(fens; nearestto=[0 0 R])[1]
    
    qpointdof = dchi.dofnums[qpoint, 3]
    cpointdof = dchi.dofnums[cpoint, 3]
    cpointdof6 = dchi.dofnums[cpoint, 6]
    
    function computetrac!(forceout, XYZ, tangents, feid, qpid)
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

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Sinusoidal loading on the surface of the shell
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64, 6, computetrac!);
    Fmag = distribloads(lfemm, vassem, geom0, dchi, fi, 2);

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
    # Fmag = fill(0.0, nfreedofs(dchi))
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
    loop!(M, K, ksi, U0, V0, nsteps*dt, dt, force!, peek)
    
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
        for i in eachindex(displacements)
            scattersysvec!(dchi, displacements[i])
            push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
            push!(times, i*dt*nbtw)
        end
        vtkwritecollection("clamp_cyl_forc_damp_expl_$n", fens, fes, times; vectors = vectors)
    end
end

function _execute_parallel_vom(n = 64, thickness = 0.01, nthr = 0)
    color = "magenta"
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
    @info "Mesh $(count(fens)) nodes, $(count(fes)) elements"

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, TriRule(1)))
    C = connectionmatrix(femm, count(fens))
    perm = symrcm(C)
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    

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
    Kasm = SM.SysmatAssemblerVecOfMatSymm(0.0)
    FEMMShellT3FFModule.stiffness(femm, Kasm, geom0, u0, Rfield0, dchi);
    M = FEMMShellT3FFModule.mass(femm, SysmatAssemblerSparseDiag(), geom0, dchi);
    
    nth = (nthr == 0 ? Base.Threads.nthreads() : nthr)
    @info "$nth threads used"
    N = Kasm.buffer_length
    chunk = Int(floor(N / nth))
    threadbuffs = ThreadBufferVOM[];
    for th in 1:nth
        rang = th < nth ? (chunk*(th-1)+1:chunk*(th)+1-1) : (chunk*(th-1)+1:N)
        push!(threadbuffs, ThreadBufferVOM(rang, fill(0.0, nfreedofs(dchi))));
    end

    # Solve
    function pwr(vom, P, threadbuffs, M)
        invM = fill(0.0, size(M, 1))
        invM .= 1.0 ./ (vec(diag(M)))
        v = rand(size(M, 1))
        w = fill(0.0, size(M, 1))
        for i in 1:30
            parvommul!(w, vom, P, threadbuffs, v)
            wn = norm(w)
            w .*= (1.0/wn)
            v .= invM .* w
            vn = norm(v)
            v .*= (1.0/vn)
        end
        parvommul!(w, vom, P, threadbuffs, v)
        sqrt((v' * w) / (v' * M * v))
    end
    vom = Kasm.buffer
    P = Kasm.elem_mat_dim
    @time omega_max = pwr(vom, P, threadbuffs, M)
    @show omega_max

    # @time evals, evecs, nconv = eigs(Symmetric(K), Symmetric(M); nev=1, which=:LM, explicittransform=:none)
    # @show sqrt(evals)
    # @show omega_max = sqrt(evals[1])

    @show dt = Float64(0.9* 2/omega_max) * (sqrt(1+ksi^2) - ksi)

    U0 = gathersysvec(dchi)
    V0 = deepcopy(U0)
    
    qpoint = selectnode(fens; nearestto=[0 -L/4 R])[1]
    cpoint = selectnode(fens; nearestto=[0 0 R])[1]
    
    qpointdof = dchi.dofnums[qpoint, 3]
    cpointdof = dchi.dofnums[cpoint, 3]
    cpointdof6 = dchi.dofnums[cpoint, 6]
    
    function computetrac!(forceout, XYZ, tangents, feid, qpid)
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
    fi = ForceIntensity(Float64, 6, computetrac!);
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
    # Fmag = fill(0.0, nfreedofs(dchi))
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
    parloop_vom!(M, Kasm, threadbuffs, ksi, U0, V0, nsteps*dt, dt, force!, peek, nthr)
    
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
        for i in eachindex(displacements)
            scattersysvec!(dchi, displacements[i])
            push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
            push!(times, i*dt*nbtw)
        end
        vtkwritecollection("clamp_cyl_forc_damp_expl_$n", fens, fes, times; vectors = vectors)
    end
end

function test_serial(ns = [4*64])
    @info "Clamped cylinder: serial"
    for n in ns
        _execute_serial(n, 0.01)
    end
    return true
end

function test_parallel(ns = [4*64], nthr = 0)
    @info "Clamped cylinder: parallel"
    for n in ns
        _execute_parallel(n, 0.01, nthr)
    end
    return true
end

function test_parallel_csr(ns = [4*64], nthr = 0)
    @info "Clamped cylinder: parallel CSR"
    for n in ns
        _execute_parallel_csr(n, 0.01, nthr)
    end
    return true
end

function test_parallel_vom(ns = [4*64], nthr = 0)
    @info "Clamped cylinder: parallel VOM"
    for n in ns
        _execute_parallel_vom(n, 0.01, nthr)
    end
    return true
end

function _explore_csr(n = 64, thickness = 0.01, nthr = 0)
    color = "green"
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
    @info "Mesh $(count(fens)) nodes, $(count(fes)) elements"

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, TriRule(1)))
    C = connectionmatrix(femm, count(fens))
    perm = symrcm(C)
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    

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
    # @show dchi.dofnums

    # Assemble the system matrix
    FEMMShellT3FFModule.associategeometry!(femm, geom0)
    SM = FinEtoolsFlexStructures.AssemblyModule
    K = FEMMShellT3FFModule.stiffness(femm, SM.SysmatAssemblerSparseCSRSymm(0.0), geom0, u0, Rfield0, dchi);
    M = FEMMShellT3FFModule.mass(femm, SysmatAssemblerSparseDiag(), geom0, dchi);

    I, J, V = findnz(K)
    starts = fill(Inf, size(K, 1))
    finishes = fill(0, size(K, 1))
    nzs = fill(0, size(K, 1))

    for i in eachindex(I)
        nzs[I[i]] += 1
        starts[I[i]] = min(starts[I[i]], J[i])
        finishes[I[i]] = max(finishes[I[i]], J[i])
    end
    bw = finishes .- starts .+ 1
    display(lineplot(1:length(bw), bw))
        @show median(bw)
        @show minimum(bw)
        @show maximum(bw)
        @show median(nzs)
        @show minimum(nzs)
        @show maximum(nzs)
    display(spy(K))

    nothing
end

function allrun(ns = [64])
    println("#####################################################")
    println("# test_serial ")
    test_serial(ns)
    # println("#####################################################")
    # println("# test_parallel ")
    # test_parallel(ns)
    println("#####################################################")
    println("# test_parallel_csr ")
    test_parallel_csr(ns)
    # println("#####################################################")
    # println("# test_parallel_vom ")
    # test_parallel_vom(ns)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing

# #####################################################   
# # test_serial
# [ Info: Clamped cylinder: serial      
# [ Info: Mesh 262400 nodes, 524288 elements  
#   2.341943 seconds (40 allocations: 97.752 MiB)       
# omega_max = 1.41323e+07
# dt = Float64((9.00000e-01 * 2) / omega_max) * (sqrt(1 + ksi ^ 2) - ksi) = 1.27368e-07 
# [ Info: 1963 steps
# 146.272046 seconds (4.32 k allocations: 1.147 GiB, 0.05% gc time)     
# [ Info: Dumping visualization 
# ##################################################### 
# # test_parallel
# [ Info: Clamped cylinder: parallel    
# [ Info: Mesh 262400 nodes, 524288 elements    
#   2.376553 seconds (40 allocations: 97.752 MiB)       
# omega_max = 1.41333e+07       
# dt = Float64((9.00000e-01 * 2) / omega_max) * (sqrt(1 + ksi ^ 2) - ksi) = 1.27359e-07 
# [ Info: 1963 steps
# [ Info: 2 threads used  
#  84.598355 seconds (176.37 k allocations: 1.158 GiB, 0.16% gc time, 0.10% compilation time)   
# [ Info: Dumping visualization 
# ##################################################### 
# # test_parallel_csr
# [ Info: Clamped cylinder: parallel CSR
# [ Info: Mesh 262400 nodes, 524288 elements    
#   2.114689 seconds (118.45 k allocations: 103.899 MiB, 29.08% compilation time)       
# omega_max = 1.41333e+07       
# dt = Float64((9.00000e-01 * 2) / omega_max) * (sqrt(1 + ksi ^ 2) - ksi) = 1.27359e-07 
# [ Info: 1963 steps
# [ Info: 2 threads used
#  82.601808 seconds (20.03 k allocations: 1.148 GiB, 0.25% gc time)    
# [ Info: Dumping visualization 
# ##################################################### 
# # test_parallel_vom   
# [ Info: Clamped cylinder: parallel VOM
# [ Info: Mesh 262400 nodes, 524288 elements    
# [ Info: 2 threads used
#   4.372176 seconds (172.21 k allocations: 95.321 MiB, 1.67% compilation time) 
# omega_max = 1.41354e+07       
# dt = Float64((9.00000e-01 * 2) / omega_max) * (sqrt(1 + ksi ^ 2) - ksi) = 1.27340e-07 
# [ Info: 1963 steps
# [ Info: 2 threads used
# 268.287995 seconds (30.07 k allocations: 1.150 GiB, 0.04% gc time)    
# [ Info: Dumping visualization 
# true      


# ##################################################### 
# # test_serial 
# [ Info: Clamped cylinder: serial      
# [ Info: Mesh 262400 nodes, 524288 elements    
#   2.360393 seconds (40 allocations: 97.752 MiB)       
# omega_max = 1.41322e+07       
# dt = Float64((9.00000e-01 * 2) / omega_max) * (sqrt(1 + ksi ^ 2) - ksi) = 1.27369e-07 
# [ Info: 1963 steps    
# 147.455232 seconds (4.32 k allocations: 1.147 GiB, 0.06% gc time)     
# [ Info: Dumping visualization 
# ##################################################### 
# # test_parallel       
# [ Info: Clamped cylinder: parallel    
# [ Info: Mesh 262400 nodes, 524288 elements    
#   2.412289 seconds (40 allocations: 97.752 MiB)       
# omega_max = 1.41318e+07       
# dt = Float64((9.00000e-01 * 2) / omega_max) * (sqrt(1 + ksi ^ 2) - ksi) = 1.27373e-07 
# [ Info: 1963 steps    
# [ Info: 4 threads used
#  80.991000 seconds (205.23 k allocations: 1.160 GiB, 0.18% gc time, 0.13% compilation time)   
# [ Info: Dumping visualization 
# ##################################################### 
# # test_parallel_csr   
# [ Info: Clamped cylinder: parallel CSR
# [ Info: Mesh 262400 nodes, 524288 elements    
#   1.988305 seconds (158.39 k allocations: 106.171 MiB, 34.71% compilation time)       
# omega_max = 1.41331e+07       
# dt = Float64((9.00000e-01 * 2) / omega_max) * (sqrt(1 + ksi ^ 2) - ksi) = 1.27361e-07 
# [ Info: 1963 steps    
# [ Info: 4 threads used
#  72.222234 seconds (20.03 k allocations: 1.148 GiB, 0.35% gc time)    
# [ Info: Dumping visualization 
# ##################################################### 
# # test_parallel_vom   
# [ Info: Clamped cylinder: parallel VOM
# [ Info: Mesh 262400 nodes, 524288 elements    
# [ Info: 4 threads used
#   3.393447 seconds (172.64 k allocations: 95.362 MiB, 2.69% compilation time) 
# omega_max = 1.41302e+07       
# dt = Float64((9.00000e-01 * 2) / omega_max) * (sqrt(1 + ksi ^ 2) - ksi) = 1.27387e-07 
# [ Info: 1963 steps    
# [ Info: 4 threads used
# 220.470659 seconds (54.74 k allocations: 1.152 GiB, 0.05% gc time)    
# [ Info: Dumping visualization 
# true
