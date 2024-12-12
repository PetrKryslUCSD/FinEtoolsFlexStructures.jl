"""
Clamped cylinder, transient vibration.

The cylinder is banged at the initial time (i.e. it is given an initial velocity).
"""
module test_loops

using LinearAlgebra
using SparseArrays
using Arpack
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: default_layout_3d, plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using PlotlyJS
using Gnuplot
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection

using BenchmarkTools
using ProfileView

E = 200e9;
nu = 0.3;
rho = 7800.0
R = 0.1
L = 0.8
ksi = 0.0
omegad = 2*pi*1000
color = "blue"
omegaf = 2*pi*10^5
tend = 15 * (2*pi/omegaf)


function loop!(M, K, ksi, U0, V0, tend, dt, force!, peek)
    U1 = deepcopy(U0)
    V1 = deepcopy(U0)
    A0 = deepcopy(U0)
    A1 = deepcopy(U0)
    F = deepcopy(U0)
    E = deepcopy(U0)
    Vp = deepcopy(U0)
    C = deepcopy(U0)
    C .= (ksi*2*omegad) .* vec(diag(M))
    invMC = deepcopy(U0)
    invMC .= 1.0 ./ (vec(diag(M)) .+ (dt/2) .* C)
    nsteps = Int64(round(tend/dt))
    if nsteps*dt < tend
        dt = tend / (nsteps+1)
    end
    t = 0.0
    
    # @btime @. $U1 = $U0 + $dt*$V0 + (($dt^2)/2)*$A0; #1
    # @btime $force!($F, $t); #2
    # @btime $F .-=  mul!($E, $K, $U1) #3
    # @btime @. $Vp = $V0 + ($dt/2) * $A0 #4
    # @btime @. $F -= $C * $Vp #5
    # @btime $A1 .= $invMC .* $F #6
    # @btime @. $A1 = $invMC * $F #6
    # @btime @inbounds for i in eachindex($A1)
    #     $A1[i] = $invMC[i] * $F[i] #6
    # end
    # @btime @. $V1 = $V0 + ($dt/2)* ($A0+$A1); #7
    # @btime $peek(1, $U0, $t) #8

    A0 .= invMC .* force!(F, t);
    peek(0, U0, t)
    @time for step in 1:nsteps
        # Displacement update
        @. U1 = U0 + dt*V0 + ((dt^2)/2)*A0; 
        # External loading
        force!(F, t);
        # Add elastic restoring forces
        F .-=  mul!(E, K, U1)
        # Add damping forces
        @. Vp = V0 + (dt/2) * A0
        @. F -= C * Vp
        # Compute the new acceleration.
        A1 .= invMC .* F
        # # Update the velocity
        @. V1 = V0 + (dt/2)* (A0+A1);
        t = t + dt
        U0, U1 = U1, U0
        V0, V1 = V1, V0
        A0, A1 = A1, A0
        peek(step, U0, t)
    end
end

cylindrical!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) = begin
    r = vec(XYZ); r[2] = 0.0;
    csmatout[:, 3] .= vec(r)/norm(vec(r))
    csmatout[:, 2] .= (0.0, 1.0, 0.0) #  this is along the axis
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

struct ThreadBuffer
    # colrange::UnitRange{Int64}
    kcolumns::SparseMatrixCSC{Float64, Int64}
    uv::SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true} 
    result::Vector{Float64}
end

struct ThreadBuff
    colrange::UnitRange{Int64}
    result::Vector{Float64}
end

function _execute(n = 8, thickness = 0.01, visualize = true)

    # @info "Mesh: $n elements per side"
    # Mesh
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
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    # Set up
    femm.drilling_mass_scale = 1.0
    femm.drilling_stiffness_scale = 1.0

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
    numberdofs!(dchi);

# Assemble the system matrix
    FEMMShellT3FFModule.associategeometry!(femm, geom0)
    K = FEMMShellT3FFModule.stiffness(femm, geom0, u0, Rfield0, dchi);
    M = FEMMShellT3FFModule.mass(femm, SysmatAssemblerSparseDiag(), geom0, dchi);
    
    U0 = gathersysvec(dchi)
    R0 = deepcopy(U0)

    
    
    @show nth = Base.Threads.nthreads()
    chunk = Int(floor(dchi.nfreedofs / nth))

    threadbuffers = ThreadBuffer[];
    for th in 1:nth
        # This thread will work on a subset of the columns.
        colrange = th < nth ? (chunk*(th-1)+1:chunk*(th)+1-1) : (chunk*(th-1)+1:dchi.nfreedofs)
        push!(threadbuffers, ThreadBuffer(K[:, colrange], view(U0, colrange), deepcopy(R0)));
    end

    function parmul!(R, tbs)
        tasks = [];
        for th in eachindex(tbs)
            tb = tbs[th]
            push!(tasks, Threads.@spawn begin 
                mul!(tb.result, tb.kcolumns, tb.uv)
            end);
        end
        Threads.wait(tasks[1]);
        R .= tbs[1].result
        for th in 2:length(tbs)
            Threads.wait(tasks[th]);
            R .+= tbs[th].result
        end
        R
    end

    threadbuffs = ThreadBuff[];
    for th in 1:nth
        colrange = th < nth ? (chunk*(th-1)+1:chunk*(th)+1-1) : (chunk*(th-1)+1:dchi.nfreedofs)
        push!(threadbuffs, ThreadBuff(colrange, deepcopy(R0)));
    end

    function parcoomul!(R, IB, JB, VB, tbs, U)
        function innercoomul!(A, IB, JB, VB, C)
            nz = length(IB)
            @assert length(JB) == nz
            @assert length(VB) == nz
            @inbounds @simd for i in 1:nz
                ir = IB[i]
                jc = JB[i]
                A[ir] += VB[i] * C[jc]
            end
            nothing
        end
        tasks = [];
        for th in eachindex(tbs)
            tb = tbs[th]
            push!(tasks, Threads.@spawn begin 
                tb.result .= zero(eltype(tb.result))
                innercoomul!(tb.result, view(IB, tb.colrange), view(JB, tb.colrange), view(VB, tb.colrange), U)
            end);
        end
        Threads.wait(tasks[1]);
        R .= tbs[1].result
        for th in 2:length(tbs)
            Threads.wait(tasks[th]);
            R .+= tbs[th].result
        end
        R
    end

    function coomul!(A, IB, JB, VB, C)
        nz = length(IB)
        @assert length(JB) == nz
        @assert length(VB) == nz
        A .= zero(eltype(A))
        @inbounds @simd for i in 1:nz
            ir = IB[i]
            jc = JB[i]
            A[ir] += VB[i] * C[jc]
        end
    end

    IK, JK, VK = findnz(K)

    oU0 = deepcopy(U0)
    oU0 .= rand()

    U0 .= oU0
    mul!(R0, K, U0)
    R01 = deepcopy(R0)

    U0 .= oU0
    parmul!(R0, threadbuffers)
    @show norm(R0 - R01)/norm(R01)

    U0 .= oU0
    parcoomul!(R0, IK, JK, VK, threadbuffs, U0)
    @show norm(R0 - R01)/norm(R01)

    U0 .= oU0
    coomul!(R0, IK, JK, VK, U0)
    @show norm(R0 - R01)/norm(R01)

    @btime mul!($R0, $K, $U0)
    @btime $coomul!($R0, $IK, $JK, $VK, $U0)
    @btime $parmul!($R0, $threadbuffers)
    @btime $parcoomul!($R0, $IK, $JK, $VK, $threadbuffs, $U0)

    # oR0 = deepcopy(R0)
    # oR0 .= rand()
    # oU0 = deepcopy(U0)
    # oU0 .= rand()

    # R0 .= oR0; U0 .= oU0
    # @. R0 = U0 + 0.1*R0 
    # R01 = deepcopy(R0)

    # R0 .= oR0; U0 .= oU0
    # @inbounds @simd for i in eachindex(U0)
    #     R0[i] = U0[i] + 0.1*R0[i]
    # end
    # @show norm(R0 - R01)

    # R0 .= oR0; U0 .= oU0
    # Threads.@threads for i in eachindex(U0)
    #     @inbounds R0[i] = U0[i] + 0.1*R0[i]
    # end
    # @show norm(R0 - R01)

    # ranges = []
    # for th in 1:nth
    #     # This thread will work on a subset of the columns.
    #     colrange = th < nth ? (chunk*(th-1)+1:chunk*(th)+1-1) : (chunk*(th-1)+1:dchi.nfreedofs)
    #     push!(ranges, colrange)
    # end
    # R0 .= oR0; U0 .= oU0
    # Threads.@threads for r in ranges
    #     @inbounds for j in r
    #         R0[j] = U0[j] + 0.1*R0[j]
    #     end
    # end
    # @show norm(R0 - R01)

    # @btime @. $R0 = $U0 + 0.1*$R0 
    # @btime for i in eachindex($U0)
    #     $R0[i] = $U0[i] + 0.1*$R0[i]
    # end
    # @btime @inbounds @simd for i in eachindex($U0)
    #     $R0[i] = $U0[i] + 0.1*$R0[i]
    # end
    # @btime Threads.@threads for i in eachindex($U0)
    #     @inbounds $R0[i] = $U0[i] + 0.1*$R0[i]
    # end
    # ranges = []
    # for th in 1:2
    #     # This thread will work on a subset of the columns.
    #     colrange = th < nth ? (chunk*(th-1)+1:chunk*(th)+1-1) : (chunk*(th-1)+1:dchi.nfreedofs)
    #     push!(ranges, colrange)
    # end
    # @btime Threads.@threads for r in $ranges
    #     @inbounds for j in r
    #         $R0[j] = $U0[j] + 0.1*$R0[j]
    #     end
    # end

    nothing
end

function test(ns = [64])
    @info "Clamped cylinder vibration"
    for n in ns
        _execute(n, 0.01, !false)
    end
    return true
end

function allrun()
    println("#####################################################")
    println("# test_convergence ")
    test()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing

 using .Main.test_loops; Main.test_loops.test(64*[2, 3, 4, 5])                                                    





