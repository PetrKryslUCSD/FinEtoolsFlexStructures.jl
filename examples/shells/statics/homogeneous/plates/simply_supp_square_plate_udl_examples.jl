# Simply supported square plate with uniform distributed load
# Note: In order to get fast convergence, the hard SS needs 
# to be adopted.
module simply_supp_square_plate_udl_examples

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

const E = 30e6
const nu = 0.3
const L = 10.0

loading(tL_ratio) = 1.0e10 * (tL_ratio)^3

# From A Triangular Plate Bending Element Based on an Energy-Orthogonal Free, 
# by Felippa and Bergan, 1986, Table 3. 
w_analyt_sol(p, L, D) = -4.06235e-3 * p * L^4 / D

function _cond(A)
    emax = eigs(A, nev=1, which=:LM, tol=1e-3)[1]
    emin = eigs(A, nev=1, which=:SM, tol=1e-3)[1]
    return (emax[1]) / (emin[1])
end

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

# function _execute_t3ff_quarter_model(n = 2, tL_ratio = 0.01, visualize = true)
#     formul = FEMMShellT3FFModule
#     thickness = L * tL_ratio;
#     D = E*thickness^3/12/(1-nu^2)
#     p = loading(tL_ratio)
#     analyt_sol = w_analyt_sol(p, L, D)

#     tolerance = L/n/1000
#     fens, fes = T3block(L/2,L/2,n,n);
#     fens.xyz = xyz3(fens)

#     mater = MatDeforElastIso(DeforModelRed3D, E, nu)

#     ocsys = CSys(3, 3, _cartesian!)

#     @info "Mesh: $n elements per side"

#     sfes = FESetShellT3()
#     accepttodelegate(fes, sfes)
#     femm = formul.make(IntegDomain(fes, TriRule(1), thickness), 
#         mater,
#         (t, h) -> t^2 / (t^2 + 0.0 * h^2),
#     )
#     stiffness = formul.stiffness
#     associategeometry! = formul.associategeometry!

#     # Construct the requisite fields, geometry and displacement
#     # Initialize configuration variables
#     geom0 = NodalField(fens.xyz)
#     u0 = NodalField(zeros(size(fens.xyz,1), 3))
#     Rfield0 = initial_Rfield(fens)
#     dchi = NodalField(zeros(size(fens.xyz,1), 6))

#     # Apply EBC's
#     # plane of symmetry perpendicular to X
#     l1 = selectnode(fens; box = Float64[0 0 0 L/2 -Inf Inf], inflate = tolerance)
#     for i in [1,5,6]
#         setebc!(dchi, l1, true, i)
#     end
#     # plane of symmetry perpendicular to Y
#     l1 = selectnode(fens; box = Float64[0 L/2 0 0 -Inf Inf], inflate = tolerance)
#     for i in [2,4,6]
#         setebc!(dchi, l1, true, i)
#     end
#     # simple support
#     l1 = selectnode(fens; box = Float64[L/2 L/2 0 L/2 -Inf Inf], inflate = tolerance)
#     for i in [3,4,6]
#         setebc!(dchi, l1, true, i)
#     end
#     l1 = selectnode(fens; box = Float64[0 L/2 L/2 L/2 -Inf Inf], inflate = tolerance)
#     for i in [3,5,6]
#         setebc!(dchi, l1, true, i)
#     end
#     # in-plane, rotations
#     l1 = selectnode(fens; box = Float64[0 L/2 0 L/2 -Inf Inf], inflate = tolerance)
#     for i in [1, 2, 6]
#         setebc!(dchi, l1, true, i)
#     end
#     applyebc!(dchi)
#     numberdofs!(dchi);

#     # Assemble the system matrix
#     associategeometry!(femm, geom0)
#     K = stiffness(femm, geom0, u0, Rfield0, dchi);

#     # Load
#     nl = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], tolerance = tolerance)
#     lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
#     fi = ForceIntensity(Float64[0, 0, -p, 0, 0, 0]);
#     F = distribloads(lfemm, geom0, dchi, fi, 2);

#     # @infiltrate
#     # Solve
#     solve_blocked!(dchi, K, F)
#     targetu =  dchi.values[nl, 3][1]
#     @info "w(center): $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

#     # Visualization
#     if visualize
#         scalars = []
#         for nc in 1:2
#             fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
#             push!(scalars, ("q$nc", fld.values))
#             fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
#             push!(scalars, ("eq$nc", fld.values))
#         end
#         vtkwrite("simply_supp_square_plate_udl-$(n)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

#         scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
#         update_rotation_field!(Rfield0, dchi)
#         plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
#         #plot_nodes(fens),
#             plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
#             dims = 1)
#         pl = render(plots)
#     end
#     return true
# end

function _execute_q4rs_quarter_model(
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

    tolerance = L / n / 1000
    fens, fes = Q4block(L / 2, L / 2, n, n)
    fens.xyz[:, 1] .-= L / 2
    fens.xyz[:, 2] .-= L / 2
    fens.xyz = xyz3(fens)
    ncenter = selectnode(fens; box=Float64[0 0 0 0 -Inf Inf], tolerance=tolerance)
    nm = selectnode(fens; box=Float64[(-L / 2 + L / n / 2) (-L / 2 + L / n / 2) (-L / 2 + L / n / 2) (-L / 2 + L / n / 2) -Inf Inf], inflate=tolerance)
    if mesh_distortion == :tweak
        shift = L / 2 / n / 4
        fens.xyz[nm[1], 1] -= shift
        fens.xyz[nm[1], 2] += shift
    elseif mesh_distortion == :rand
        rng = Xoshiro(0)
        bfes = meshboundary(fes)
        cn = connectednodes(bfes)
        interior = setdiff(1:count(fens), cn)
        interior = setdiff(interior, ncenter)
        shift = L / 2 / n / 10
        for j in interior
            fens.xyz[j, 1] += 2 * (rand(rng) - 0.5) * shift
            fens.xyz[j, 2] += 2 * (rand(rng) - 0.5) * shift
        end
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    ocsys = CSys(3, 3, _cartesian!)

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
    # edges parallel to Y, soft support: rotations in plane are free; 
    # hard simple support: rotation orthogonal to edge = 0
    dof = [1, 2, 3, 6]
    simple_support == :hard && (dof = [1, 2, 3, 4, 6])
    l1 = selectnode(fens; box=Float64[(-L / 2) (-L / 2) (-L / 2) 0 -Inf Inf], inflate=tolerance)
    edges_parallel_to_y = (l1, dof)
    # edges parallel to X, soft support: rotations in plane are free; 
    # hard simple support: rotation orthogonal to edge = 0
    dof = [1, 2, 3, 6]
    simple_support == :hard && (dof = [1, 2, 3, 5, 6])
    l1 = selectnode(fens; box=Float64[(-L / 2) 0 (-L / 2) (-L / 2) -Inf Inf], inflate=tolerance)
    edges_parallel_to_x = (l1, dof)
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box=Float64[-Inf Inf 0 0 -Inf Inf], inflate=tolerance)
    symmetry_orthogonal_to_y = (l1, [1, 2, 4, 6])
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box=Float64[0 0 -Inf Inf -Inf Inf], inflate=tolerance)
    symmetry_orthogonal_to_x = (l1, [1, 2, 5, 6])
    for e in [edges_parallel_to_y, edges_parallel_to_x, symmetry_orthogonal_to_y, symmetry_orthogonal_to_x]
        l1, dof = e
        for i in dof
            setebc!(dchi, l1, true, i)
        end
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    AE = AbaqusExporter("ss-qua-$(simple_support)-$(mesh_distortion)-tL=$(tL_ratio)-n=$(n).inp")
    HEADING(AE, "Simply supported square plate with uniform distributed load")
    PART(AE, "part1")
    NODE(AE, fens.xyz)
    @assert nodesperelem(fes) == 4
    ELEMENT(
        AE,
        "S4",
        "AllElements",
        connasarray(fes),
    )
    NSET_NSET(AE, "edges_parallel_to_y", edges_parallel_to_y[1])
    NSET_NSET(AE, "edges_parallel_to_x", edges_parallel_to_x[1])
    NSET_NSET(AE, "symmetry_orthogonal_to_y", symmetry_orthogonal_to_y[1])
    NSET_NSET(AE, "symmetry_orthogonal_to_x", symmetry_orthogonal_to_x[1])
    SHELL_SECTION(AE, "elasticity", "AllElements", thickness)
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    ORIENTATION(AE, "GlobalOrientation", vec([1.0 0 0]), vec([0 1.0 0]))
    INSTANCE(AE, "INSTNC1", "PART1")
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    STEP_PERTURBATION_STATIC(AE)
    DLOAD(AE, "ASSEM1.INSTNC1.AllElements", [0, 0, -p])   
    for d in edges_parallel_to_y[2]
        BOUNDARY(AE, "ASSEM1.INSTNC1.edges_parallel_to_y", d)
    end
    for d in edges_parallel_to_x[2]
        BOUNDARY(AE, "ASSEM1.INSTNC1.edges_parallel_to_x", d)
    end
    for d in symmetry_orthogonal_to_y[2]
        BOUNDARY(AE, "ASSEM1.INSTNC1.symmetry_orthogonal_to_y", d)
    end
    for d in symmetry_orthogonal_to_x[2]
        BOUNDARY(AE, "ASSEM1.INSTNC1.symmetry_orthogonal_to_x", d)
    end
    FIELD_OUTPUT(AE, "PRESELECT", "LE, S, SF")
    END_STEP(AE)
    close(AE)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Load
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(Float64[0, 0, -p, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 2)

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    targetu = dchi.values[ncenter, 3][1]
    @info "w(center): $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

    # Visualization
    if visualize
        vtkwrite("sqssudlq-q4rs-$(simple_support)-$(mesh_distortion)-tL=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-uur.vtu", fens, fes; vectors=[("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        # ocsys = CSys(3)
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys=ocsys)
            push!(scalars, ("m$nc", fld.values))
        end
        vtkwrite("sqssudlq-q4rs-$(simple_support)-$(mesh_distortion)-tL=$(tL_ratio)-s=$(stab_alpha)-m.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        scalars = []
        connectivity, points, values = elementwise_arrays(femm, geom0, dchi, 1, tolerance)
        push!(scalars, ("q1", values))
        connectivity, points, values = elementwise_arrays(femm, geom0, dchi, 2, tolerance)
        push!(scalars, ("q2", values))
        vtkwrite("sqssudlq-q4rs-$(simple_support)-$(mesh_distortion)-tL=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-q.vtu", connectivity, points, VTKWrite.Q4; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        # scalars = []
        # for nc in 1:2
        # fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys=ocsys)
        # push!(scalars, ("q$nc", fld.values))
        # end
        # vtkwrite("sqssudlq-q4rs-$(simple_support)-tL=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-q.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])

    end
    return targetu, analyt_sol
end

# function test_convergence_quarter()
#     formul = FEMMShellT3FFModule
#     tL_ratio = 0.1
#     @info "Simply supported square plate with uniform load, T3FF"
#     @info "thickness/length = $tL_ratio"
#     # for n in [32,  ]
#     for n in [2, 4, 8, 16, 32, 64]
#         _execute_t3ff_quarter_model(n, tL_ratio, true)
#     end
#     formul = FEMMShellQ4RSModule
#     tL_ratio = 0.1
#     @info "Simply supported square plate with uniform load, Q4RS"
#     @info "thickness/length = $tL_ratio"
#     # for n in [32,  ]
#     for n in [2, 4, 8, 16, 32, 64]
#         _execute_Q4RS_quarter_model(n, tL_ratio, true)
#     end
#     return true
# end

function _execute_t3ff_full_model(
    n=2,
    tL_ratio=0.01,
    simple_support=:hard,
    stab_alpha=0.1,
    mesh_distortion=:none,
    visualize=true
)
    formul = FEMMShellT3FFModule
    thickness = L * tL_ratio
    D = E / 12 / (1 - nu^2) * thickness^3
    p = loading(tL_ratio)
    analyt_sol = w_analyt_sol(p, L, D)

    tolerance = L / n / 1000
    fens, fes = T3block(L, L, n, n)
    fens.xyz = xyz3(fens)
    fens.xyz[:, 1] .-= L / 2
    fens.xyz[:, 2] .-= L / 2
    ncenter = selectnode(fens; box=Float64[0 0 0 0 -Inf Inf], tolerance=tolerance)
    nm = selectnode(fens; box=Float64[-L / 2 + L / n -L / 2 + L / n -L / 2 + L / n -L / 2 + L / n -Inf Inf], inflate=tolerance)
    if mesh_distortion == :tweak
        shift = L / n / 4
        fens.xyz[nm[1], 1] -= shift
        fens.xyz[nm[1], 2] += shift
    elseif mesh_distortion == :rand
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

    @info "Mesh: $n elements per side"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness),
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
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Load
    nl = selectnode(fens; box=Float64[0 0 0 0 -Inf Inf], tolerance=tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64[0, 0, -p, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 2)

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    targetu = dchi.values[nl, 3][1]
    @info "w(center): $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"
    epsrel = abs(targetu - analyt_sol) / abs(analyt_sol)
    @info "Digits of accuracy: $(-log10(epsrel))"
    # Visualization
    if visualize
        vtkwrite("sqssudlq-t3ff-full-$(simple_support)-$(mesh_distortion)-tL=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-uur.vtu", fens, fes; vectors=[("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        # ocsys = CSys(3)
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys=ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
            push!(scalars, ("m$nc", fld.values))
        end
        vtkwrite("sqssudlq-t3ff-full-$(simple_support)-$(mesh_distortion)-tL=$(tL_ratio)-s=$(stab_alpha)-m.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys=ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
            push!(scalars, ("q$nc", fld.values))
        end
        vtkwrite("sqssudlq-t3ff-full-$(simple_support)-$(mesh_distortion)-tL=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-q.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
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
    elseif mesh_distortion == :rand
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
        vtkwrite("sqssudlq-q4rs-full-$(simple_support)-$(mesh_distortion)-tL=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-uur.vtu", fens, fes; vectors=[("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        # ocsys = CSys(3)
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys=ocsys)
            push!(scalars, ("m$nc", fld.values))
        end
        vtkwrite("sqssudlq-q4rs-full-$(simple_support)-$(mesh_distortion)-tL=$(tL_ratio)-s=$(stab_alpha)-m.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        scalars = []
        connectivity, points, values = elementwise_arrays(femm, geom0, dchi, 1, tolerance)
        push!(scalars, ("q1", values))
        connectivity, points, values = elementwise_arrays(femm, geom0, dchi, 2, tolerance)
        push!(scalars, ("q2", values))
        vtkwrite("sqssudlq-q4rs-full-$(simple_support)-$(mesh_distortion)-tL=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-q.vtu", connectivity, points, VTKWrite.Q4; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        # scalars = []
        # for nc in 1:2
        # fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys=ocsys)
        # push!(scalars, ("q$nc", fld.values))
        # end
        # vtkwrite("sqssudlq-q4rs-full-$(simple_support)-tL=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-q.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])

    end
    return targetu, analyt_sol
end

function test_convergence_full()
    formul = FEMMShellT3FFModule
    stab_alpha = 0.1
    visualize = false
    support = :soft
    tweak_mesh = false
    randomize_mesh = false
    tL_ratio = 0.00001
    tL_ratio = 0.001
    @info "Simply supported square plate with uniform load, T3FF"
    @info "thickness/length = $tL_ratio formulation=$(formul)"
    # for n in [32,  ]
    for n in [2, 4, 8, 16, 32, 64, 128, 256, 512]
        @info "---------------------------------------------------------------------"
        _execute_t3ff_full_model(n, tL_ratio, support, stab_alpha, visualize)
        formul = FEMMShellT3FFModule
    end
    formul = FEMMShellQ4RSModule
    tL_ratio = 0.00001
    tL_ratio = 0.001
    @info "Simply supported square plate with uniform load, Q4RS"
    @info "thickness/length = $tL_ratio formulation=$(formul)"
    # for n in [32,  ]
    for n in [2, 4, 8, 16, 32, 64, 128, 256, 512]
        @info "---------------------------------------------------------------------"
        _execute_q4rs_full_model(n, tL_ratio, support, stab_alpha, mesh_distortion, visualize)
    end
    return true
end

const colors = ["black", "red", "blue", "green", "magenta", "cyan", "orange", "purple", "brown", "pink"]
const marks = ["x", "o", "square", "+", "diamond", "star", "pentagon", "dtriangle", "oplus", "otimes"]

function start_case()
    objects = []
    return objects
end

function plot_case_approx_error!(objects, j, nu_errs, ns, stab_alphas)
    @pgf p = PGFPlotsX.Plot(
        {
            color = colors[j],
            mark = "$(marks[j])",
            line_width = 1.0,
        },
        Coordinates([v for v in zip(L ./ ns, abs.(nu_errs))])
    )
    push!(objects, p)
    push!(objects, LegendEntry("$(stab_alphas[j])"))
end

function display_case_approx_error(c, objects)
    @pgf ax = Axis(
        {
            xlabel = "Normalized element size [ND]",
            ylabel = "Normalized approximate error [ND]",
            # ymin = -1.1e-5,
            # ymax = 1.1e-5,
            # xmin = 0.0,
            # xmax = 1.0,
            xmode = "log",
            ymode = "log",
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
    pgfsave("$(c.support)-$(c.mesh_distortion)-approx_errors-tL_ratio=$(c.tL_ratio).pdf", ax)
end

function plot_case_result!(objects, j, res, ns, stab_alphas)
    @pgf p = PGFPlotsX.Plot(
        {
            color = colors[j],
            mark = "$(marks[j])",
            line_width = 1.0,
        },
        Coordinates([v for v in zip(L ./ ns, res)])
    )
    push!(objects, p)
    push!(objects, LegendEntry("$(stab_alphas[j])"))
end

function display_case_result(c, objects)
    @pgf ax = Axis(
        {
            xlabel = "Normalized element size [ND]",
            ylabel = "Result, normalized [ND]",
            # ymin = -1.1e-5,
            # ymax = 1.1e-5,
            # xmin = 0.0,
            # xmax = 1.0,
            xmode = "log",
            # ymode = "log",
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
    pgfsave("$(c.support)-$(c.mesh_distortion)-results-tL_ratio=$(c.tL_ratio).pdf", ax)
end

function test_convergence_full_thickness(tL_ratios=[1.0e-1, 1.0e-2, 1.0e-4, 1.0e-8])
    stab_alpha = 0.001
    visualize = false
    support = :hard
    mesh_distortion = :rand
    ns = [2, 4, 8, 16, 32, 64]
    formul = FEMMShellT3FFModule
    @info "Simply supported square plate with uniform load, T3FF --------------------"
    for tL_ratio in tL_ratios
        @info "thickness/length = $tL_ratio"
        for n in ns
            _execute_t3ff_full_model(n, tL_ratio, support, stab_alpha, visualize)
            formul = FEMMShellT3FFModule
        end
    end
    formul = FEMMShellQ4RSModule
    @info "Simply supported square plate with uniform load, Q4RS --------------------"
    for tL_ratio in tL_ratios
        @info "thickness/length = $tL_ratio"
        for n in ns
            _execute_q4rs_full_model(n, tL_ratio, support, stab_alpha, mesh_distortion, visualize)
        end
    end
    return true
end

function test_convergence_quarter_thickness(tL_ratios=[1.0e-1, 1.0e-4])
    visualize = true
    support = :hard
    # mesh_distortion = :rand
    for mesh_distortion in [:tweak, :rand, :none]
        @info "Mesh distortion: $mesh_distortion"
        ns = [16]
        for stab_alpha = [0, 0.1]
            @info "Simply supported square plate with uniform load, Q4RS, stab_alpha=$stab_alpha  ----------------"
            for tL_ratio in tL_ratios
                @info "thickness/length = $tL_ratio"
                for n in ns
                    _execute_q4rs_quarter_model(n, tL_ratio, support, stab_alpha, mesh_distortion, visualize)
                end
            end
        end
    end
    return true
end

function test_stabilization_q4rs(tL_ratios=[1.0e-2, 1.0e-4, 1.0e-8])
    visualize = true
    support = :hard
    mesh_distortion = :rand
    n = 16
    formul = FEMMShellQ4RSModule
    @info "Simply supported square plate with uniform load, Q4RS --------------------"
    for tL_ratio in tL_ratios
        @info ">>>>>>>>>>>>> thickness/length = $tL_ratio"
        for stab_alpha in [0.0, 0.000001, 0.1]
            @info "  stab_alpha = $stab_alpha"
            _execute_q4rs_full_model(n, tL_ratio, support, stab_alpha, mesh_distortion, visualize)
        end
    end
    return true
end

function test_stabilization_t3ff(tL_ratios=[1.0e-2, 1.0e-4, 1.0e-8])
    stab_alpha = 0.001
    visualize = true
    support = :hard
    mesh_distortion = :tweak
    n = 16
    formul = FEMMShellT3FFModule
    @info "Simply supported square plate with uniform load, T3FF --------------------"
    for tL_ratio in tL_ratios
        @info ">>>>>>>>>>>>> thickness/length = $tL_ratio"
        for stab_alpha in [0.0, 0.000001, 0.001, 0.1]
            @info "  stab_alpha = $stab_alpha"
            _execute_t3ff_full_model(n, tL_ratio, support, stab_alpha, mesh_distortion, visualize)
        end
    end
    return true
end


struct CaseData
    tL_ratio::Float64
    support::Symbol
    mesh_distortion::Symbol
    u_sol::Float64
end

const case_data = [
    CaseData(1 / 100, :hard, :rand, 0.0),
    CaseData(1 / 10000, :hard, :rand, 0.0),
    CaseData(1 / 1000000, :hard, :rand, 0.0),
    CaseData(1 / 100, :hard, :none, 0.0),
    CaseData(1 / 10000, :hard, :none, 0.0),
    CaseData(1 / 1000000, :hard, :none, 0.0)
]

function graph_convergence_q4rs(
    ns=[2, 4, 8, 16, 32, 64,],
    stab_alphas=[0.0, 0.000001, 0.001, 0.1]
)
    @info "--------------------------------------------------"
    @info "Convergence study for simply supported plate"
    for c in case_data
        @info "--------------------------------------------------"
        @info "Convergence study for case t/L=$(c.tL_ratio) "
        or = _convergence_data_q4rs(c, ns, stab_alphas)
        objects_errors = start_case()
        objects_displacements = start_case()
        for (j, (_, m, r)) in enumerate(or)
            @info "  stab_fun=$(m)"
            aprox_u_sols = [r[i][1] for i in eachindex(r)]
            exact_u_sols = [r[i][2] for i in eachindex(r)]
            nu_errs = [abs(aprox_u_sols[i+1] - aprox_u_sols[i]) / aprox_u_sols[end] for i in 1:length(aprox_u_sols)-1]
            nus = [aprox_u_sols[i] / exact_u_sols[i] for i in 1:length(aprox_u_sols)]
            @info "  Displ. solutions: $(round.(aprox_u_sols, digits = 7))"
            @info "  Displ. norm. approximate errors: $(round.(nu_errs, digits = 4))"
            plot_case_approx_error!(objects_errors, j, nu_errs, ns, stab_alphas)
            plot_case_result!(objects_displacements, j, nus, ns, stab_alphas)
        end
        display_case_approx_error(c, objects_errors)
        display_case_result(c, objects_displacements)
    end
    return true
end

function _convergence_data_q4rs(c, ns, stab_alphas)
    @info "Q4RS elements, t/L=$(c.tL_ratio), support=$(c.support), mesh_distortion=$(c.mesh_distortion)"
    visualize = false
    all_results = []
    for stab_alpha in stab_alphas
        @info "  stab_alpha = $stab_alpha"
        results = []
        for n in ns
            r = _execute_q4rs_full_model(n, c.tL_ratio, c.support, stab_alpha, c.mesh_distortion, visualize)
            push!(results, r)
        end
        push!(all_results, (case=c, stab_alpha=stab_alpha, results=results))
    end
    return all_results
end


function allrun()
    # println("#####################################################")
    # println("# test_stabilization_t3ff ")
    # test_stabilization_t3ff()
    # println("#####################################################")
    # println("# test_stabilization_q4rs ")
    # test_stabilization_q4rs()
     println("#####################################################")
    println("# test_convergence_quarter_thickness ")
    test_convergence_quarter_thickness()
    # println("#####################################################")
    # println("# test_convergence_full ")
    # test_convergence_full()
    # println("#####################################################")
    # println("# test_convergence_full_thickness ")
    # test_convergence_full_thickness()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
