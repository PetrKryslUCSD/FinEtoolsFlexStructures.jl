# Simply supported square plate with uniform distributed load
# Note: In order to get fast convergence, the hard SS needs 
# to be adopted.
module razzaque_plate_udl_examples

using Arpack
using LinearAlgebra: norm
using Random: rand, Xoshiro
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.FEMMShellQ4RSModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule: VTKWrite
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using PGFPlotsX

const E = 1085.0
const nu = 0.31
const L = 10.0
const thickness = 0.01
const fz = 1.0
const Db = E / 12 / (1 - nu^2) * thickness^3
const wc = -7.945 * (fz * L^4 / 1000 / Db)

function _cartesian!(csmatout, XYZ, tangents, feid, qpid)
    csmatout[:, 1] .= (1.0, 0.0, 0.0)
    csmatout[:, 2] .= (0.0, 1.0, 0.0)
    csmatout[:, 3] .= (0.0, 0.0, 1.0)
    return csmatout
end

function _closest(loc, x, tol)
    mind = Inf
    c = 0 # no closest point (i.e. no point within tolerance) 
    for i in axes(x, 1)
        d = norm(vec(x[i, :]) - vec(loc))
        if d < tol && d < mind
            mind = d
            c = i
        end
    end
    return c
end

struct Idat
    nc::Int
    v::Vector{Float64}
    tol::Float64
end

function _inspector(idat, j, conn, x, out, loc)
    c = _closest(loc, x, idat.tol)
    if c > 0
        idat.v[c] = out[idat.nc]
    end
    return idat
end

function elementwise_arrays(femm, geom0, dchi, nc, tol)
    fes = finite_elements(femm)
    connectivity = zeros(Int, count(fes), 4)
    points = zeros(4 * count(fes), 3)
    values = zeros(4 * count(fes))
    idat = Idat(nc, zeros(4), tol)
    ecoords = zeros(4, 3)
    for i in eachindex(fes)
        gathervalues_asmat!(geom0, ecoords, fes.conn[i])
        inspectintegpoints(femm, geom0, dchi, [i], _inspector, idat, :shear)
        r = (i - 1) * 4 .+ (1:4)
        connectivity[i, :] .= r
        points[r, :] = ecoords
        values[r] = idat.v
        idat.v .= 0
    end
    return connectivity, points, values
end

function _execute_t3ff_full_model(
    n=2,
    visualize=true
)
    formul = FEMMShellT3FFModule
    
    tolerance = L / n / 1000
    fens, fes = T3block(L, L, n, n)
    fens.xyz = xyz3(fens)
    ncenter = selectnode(fens; box=Float64[L/2 L/2 L/2 L/2 -Inf Inf], tolerance=tolerance)
    y0l = selectnode(fens; box=Float64[0 L 0 0 -Inf Inf], tolerance=tolerance)
    y1l = selectnode(fens; box=Float64[0 L L L -Inf Inf], tolerance=tolerance)

    for j in eachindex(fens)
        s = sin(pi / 180 * 30) * L
        fens.xyz[j, 1] += s * fens.xyz[j, 2] / L
        fens.xyz[j, 2] *= cos(pi / 180 * 30)
    end


    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    ocsys = CSys(3, 3, _cartesian!)

    @info "Mesh: $n elements per side"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness),
        mater,
    )
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # edges parallel to Y, soft support: rotations in plane are free; hard simple support: rotation orthogonal to edge = 0
    dof = [1, 2, 3, 5, 6]
    for i in dof
        setebc!(dchi, y0l, true, i)
        setebc!(dchi, y1l, true, i)
    end
    
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Load
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64[0, 0, -fz, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 2)

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    targetu = dchi.values[ncenter, 3][1]
    @show targetu, wc
    @info "w(center): $(round(targetu, digits=8)),  $(round(targetu/wc, digits = 4)*100)%"
    epsrel = abs(targetu - wc) / abs(wc)
    @info "Digits of accuracy: $(-log10(epsrel))"
    # Visualization
    if visualize
        vtkwrite("razz-n=$(n)-uur.vtu", fens, fes; vectors=[("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        # ocsys = CSys(3)
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys=ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
            push!(scalars, ("m$nc", fld.values))
        end
        vtkwrite("razz-m.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys=ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
            push!(scalars, ("q$nc", fld.values))
        end
        vtkwrite("razz-n=$(n)-q.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    end
    return true
end

function _execute_q4rs_full_model(
    n=2,
    tL_ratio=0.01,
    simple_support=:hard,
    stab_alpha=0.1,
    mesh_distortion=:none,
    visualize=true
)
    formul = FEMMShellQ4RSModule
    thickness = L * tL_ratio
    D = E / 12 / (1 - nu^2) * thickness^3
    p = loading(tL_ratio)
    analyt_sol = w_analyt_sol(p, L, D)

    tolerance = L / n / 100
    fens, fes = Q4block(L, L, n, n)
    fens.xyz = xyz3(fens)
    fens.xyz[:, 1] .-= L / 2
    fens.xyz[:, 2] .-= L / 2
    ncenter = selectnode(fens; box=Float64[0 0 0 0 -Inf Inf], tolerance=tolerance)
    nm = selectnode(fens; box=Float64[-L / 2 + L / n -L / 2 + L / n -L / 2 + L / n -L / 2 + L / n -Inf Inf], inflate=tolerance)
    if mesh_distortion == :tweak
        shift = L / n / 4
        fens.xyz[nm[1], 1] -= shift
        fens.xyz[nm[1], 2] += shift
    elseif mesh_distortion == :randomize
        rng = Xoshiro(0)
        bfes = meshboundary(fes)
        cn = connectednodes(bfes)
        interior = setdiff(1:count(fens), cn)
        interior = setdiff(interior, ncenter)
        shift = L / n / 5
        for j in interior
            fens.xyz[j, 1] += 2 * (rand(rng) - 0.5) * shift
            fens.xyz[j, 2] += 2 * (rand(rng) - 0.5) * shift
        end
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    ocsys = CSys(3, 3, _cartesian!)

    # @info "Mesh: $n elements per side"

    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    ir = Simpson13Rule(2) # GaussRule(2, 2)
    femm = formul.make(IntegDomain(fes, ir, thickness),
        mater,
        (t, h) -> t^2 / (t^2 + stab_alpha * h^2),
    )
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # edges parallel to Y, soft support: rotations in plane are free; hard simple support: rotation orthogonal to edge = 0
    dof = [1, 2, 3, 6]
    simple_support == :hard && (dof = [1, 2, 3, 4, 6])
    l1 = selectnode(fens; box=Float64[-L / 2 -L / 2 -Inf Inf -Inf Inf], inflate=tolerance)
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box=Float64[+L / 2 +L / 2 -Inf Inf -Inf Inf], inflate=tolerance)
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    # edges parallel to X, soft support: rotations in plane are free; hard simple support: rotation orthogonal to edge = 0
    dof = [1, 2, 3, 6]
    simple_support == :hard && (dof = [1, 2, 3, 5, 6])
    l1 = selectnode(fens; box=Float64[-Inf Inf -L / 2 -L / 2 -Inf Inf], inflate=tolerance)
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box=Float64[-Inf Inf +L / 2 +L / 2 -Inf Inf], inflate=tolerance)
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    setebc!(dchi, 1:count(fens), true, 6)
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)
    # @info "stab_alpha = $stab_alpha, Condition number: $(_cond(matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:ff]))"

    # Load
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(Float64[0, 0, -p, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 2)

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    targetu = dchi.values[ncenter, 3][1]
    @info "w(center): $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"
    epsrel = abs(targetu - analyt_sol) / abs(analyt_sol)
    @info "Digits of accuracy: $(-log10(epsrel))"
    # Visualization
    if visualize
        vtkwrite("ss_sqpl_udl-q4rs-full-$(simple_support)-$(mesh_distortion)-r=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-uur.vtu", fens, fes; vectors=[("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        # ocsys = CSys(3)
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys=ocsys)
            push!(scalars, ("m$nc", fld.values))
        end
        vtkwrite("ss_sqpl_udl-q4rs-full-$(simple_support)-$(mesh_distortion)-r=$(tL_ratio)-s=$(stab_alpha)-m.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        scalars = []
        connectivity, points, values = elementwise_arrays(femm, geom0, dchi, 1, tolerance)
        push!(scalars, ("q1", values))
        connectivity, points, values = elementwise_arrays(femm, geom0, dchi, 2, tolerance)
        push!(scalars, ("q2", values))
        vtkwrite("ss_sqpl_udl-q4rs-full-$(simple_support)-$(mesh_distortion)-r=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-q.vtu", connectivity, points, VTKWrite.Q4; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        # scalars = []
        # for nc in 1:2
        # fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys=ocsys)
        # push!(scalars, ("q$nc", fld.values))
        # end
        # vtkwrite("ss_sqpl_udl-q4rs-full-$(simple_support)-r=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-q.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])

    end
    return targetu, analyt_sol
end

function test_convergence()
    formul = FEMMShellT3FFModule
    @info "Razzaque skew plate with uniform load, T3FF"
    visualize = true
    for n in [2, 4, 8, 16, ]
        @info "---------------------------------------------------------------------"
        _execute_t3ff_full_model(n, visualize)
        formul = FEMMShellT3FFModule
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
