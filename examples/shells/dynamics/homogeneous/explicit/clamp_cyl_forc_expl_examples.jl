"""
Clamped cylinder, transient vibration.

The cylinder is banged at the initial time (i.e. it is given an initial velocity).
"""
module clamp_cyl_forc_expl_examples

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
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

function loop!(M, K, U0, V0, tend, dt, force!, peek)
    U1 = deepcopy(U0)
    V1 = deepcopy(U0)
    A0 = deepcopy(U0)
    A1 = deepcopy(U0)
    F = deepcopy(U0)
    E = deepcopy(U0)
    invM = deepcopy(U0)
    invM .= 1.0 ./ vec(diag(M))
    nsteps = Int64(round(tend/dt))
    if nsteps*dt < tend
        dt = tend / (nsteps+1)
    end
    t = 0.0
    A0 .= M\force!(F, t);
    peek(0, U0, t)
    @time for step in 1:nsteps
        # External loading
        F = force!(F, t);
        # Displacement update
        @. U1 = U0 + dt*V0 + ((dt^2)/2)*A0; 
        #
        F .-=  mul!(E, K, U1)
        # Compute the new acceleration.
        A1 .= invM .* F
        # # Update the velocity
        @. V1 = V0 + (dt/2)* (A0+A1);
        t = t + dt
        U0, U1 = U1, U0
        V0, V1 = V1, V0
        A0, A1 = A1, A0
        peek(step, U0, t)
    end
end

E = 200e9;
nu = 0.3;
rho = 7800.0
R = 0.1
L = 0.8
omegaf = 2*pi*10^4
tend = 15 * 2*pi/omegaf

cylindrical!(csmatout, XYZ, tangents, feid, qpid) = begin
    r = vec(XYZ); r[2] = 0.0;
    csmatout[:, 3] .= vec(r)/norm(vec(r))
    csmatout[:, 2] .= (0.0, 1.0, 0.0) #  this is along the axis
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

function _execute(n = 8, thickness = 0.001, visualize = true)

    @info "Mesh: $n elements per side"
    # Mesh
    tolerance = min(R, L) / n  / 100
    
    tolerance = R/n/1000
    fens, fes = T3block(360.0, L, n, 4*n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a/180*pi), y-L/2, R*cos(a/180*pi))
    end
    fens, fes = mergenodes(fens, fes, tolerance)
    bfes = meshboundary(fes)
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    # Set up
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
    M = FEMMShellT3FFModule.mass(femm, geom0, dchi);
    # Check that the mass matrix is diagonal
    # Make sure the matrix is truly diagonal: delete all tiny off-diagonals
    I, J, V = findnz(M)
    for i in eachindex(I)
        if I[i] != J[i]
            V[i] = 0.0
        end
    end
    M = sparse(I, J, V, dchi.nfreedofs, dchi.nfreedofs)
    
# Solve
    function pwr(A, B)
        fa = cholesky(B)
        v = rand(size(B, 1))
        w = fill(0.0, size(B, 1))
        for i in 1:30
            w = (A * v)
            w = w / norm(w)
            v .= fa \ w
            v = v / norm(v)
        end
        sqrt((v' * A * v) / (v' * B * v))
    end
    @time omega_max = pwr(K, M)
    @show omega_max

    # @time evals, evecs, nconv = eigs(Symmetric(K), Symmetric(M); nev=1, which=:LM, explicittransform=:none)
    # @show sqrt(evals)
    # @show omega_max = sqrt(evals[1])

    @show dt = Float64(0.9* 2/omega_max)

    U0 = gathersysvec(dchi)
    V0 = deepcopy(U0)
    
    qpoint = selectnode(fens; nearestto=[0 -L/4 R])[1]
    cpoint = selectnode(fens; nearestto=[0 0 -R])[1]
    applyebc!(dchi)
    numberdofs!(dchi);
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

    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64, 6, computetrac!);
    Fmag = distribloads(lfemm, geom0, dchi, fi, 2);

    function force!(F, t)
        F .= sin(omegaf*t) .* Fmag
        # @show norm(F)
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
    loop!(M, K, U0, V0, nsteps*dt, dt, force!, peek)
    
    # @gp  "set terminal windows 0 "  :-

    color = "red"
    @gp  :- collect(0.0:dt:(nsteps*dt)) cdeflections " lw 2 lc rgb '$color' with lines title 'Deflection at the center' "  :-

    @gp  :- "set xlabel 'Time'" :-
    @gp  :- "set ylabel 'Deflection'" :-
    @gp  :- "set title 'Free-floating plate'"


# Visualization
    @show length(displacements)
    box = boundingbox(fens.xyz)
    box[1] -= 0.05*L
    box[2] += 0.05*L
    box[3] -= 0.05*L
    box[4] += 0.05*L
    box[5] -= 0.5*L
    box[6] += 0.5*L
    
    tbox = plot_space_box(reshape(box, 2, 3))
    tenv0 = plot_midsurface(fens, fes; x=geom0.values, u=0.0 .* dchi.values[:, 1:3], R=Rfield0.values, facecolor="rgb(155, 155, 155)", opacity=0.3);

    plots = cat(tbox, tenv0; dims=1)
    layout = default_layout_3d(;width=600, height=600)
    layout[:scene][:aspectmode] = "data"
    pl = render(plots; layout=layout, title = "Step ")
    sleep(2.5)
    
    scale = 10^6;
    for i in eachindex(displacements)
        scattersysvec!(dchi, scale .* displacements[i])
        Rfield1 = deepcopy(Rfield0)
        update_rotation_field!(Rfield1, dchi)
        tenv1 = plot_midsurface(fens, fes; x=geom0.values, u=dchi.values[:, 1:3], R=Rfield1.values, facecolor="rgb(50, 125, 125)",  lighting = attr(diffuse=0.3, fresnel=4, specular=0.5), lightposition = attr(x = 400, y = -100, z = -5));
        plots = cat(tbox, tenv0, tenv1; dims=1)
        pl.plot.layout[:title] = "Step $(i)"
        react!(pl, plots, pl.plot.layout)
        sleep(0.115)
    end

end

function test_convergence()
    @info "Clamped cylinder vibration"
    for n in [64] 
        _execute(n, 0.01, !false)
    end
    return true
end

function allrun()
    println("#####################################################")
    println("# test_convergence ")
    test_convergence()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing




