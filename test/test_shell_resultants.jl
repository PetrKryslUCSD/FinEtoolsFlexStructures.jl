
module mcompostat0

# The barrel vault (Scordelis-Lo) roof is one of the benchmarks for linear elastic
# analysis of shells. 

# Here we test the computation of the resultant forces and moments.

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

# analytical solution for the vertical deflection and the midpoint of the
# free edge 
analyt_sol = -0.3024
# Parameters:
E = 4.32e8
nu = 0.0
thickness = 0.25 # geometrical dimensions are in feet
R = 25.0
L = 50.0

cylindrical!(csmatout, XYZ, tangents, feid, qpid) = begin
    r = vec(XYZ)
    r[2] = 0.0
    r[3] += R
    csmatout[:, 3] .= vec(r) / norm(vec(r))
    csmatout[:, 2] .= (0.0, 1.0, 0.0) #  this is along the axis
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

function test_homogeneous(visualize = true)
    formul = FEMMShellT3FFModule
    n = 8
    tolerance = R / n / 10
    fens, fes = T3block(40 / 360 * 2 * pi, L / 2, n, n, :b)
    fens.xyz = xyz3(fens)
    for i = 1:count(fens)
        a = fens.xyz[i, 1]
        y = fens.xyz[i, 2]
        fens.xyz[i, :] .= (R * sin(a), y, R * (cos(a) - 1))
    end

    # conn = fill(0, count(fes), 3)
    # for i in 1:count(fes)
    #     conn[i, :] .= fes.conn[i][2], fes.conn[i][3], fes.conn[i][1]
    # end
    # fromarray!(fes, conn)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    femm.drilling_stiffness_scale = 1.0e-4
    femm.mult_el_size = 5 / 12
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # rigid diaphragm
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1, 3, 5]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf L / 2 L / 2 -Inf Inf], inflate = tolerance)
    for i in [2, 4, 6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1, 5, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Midpoint of the free edge
    nl = selectnode(
        fens;
        box = Float64[sin(40 / 360 * 2 * pi) * 25 sin(40 / 360 * 2 * pi) * 25 L / 2 L / 2 -Inf Inf],
        inflate = tolerance,
    )
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity([0, 0, -90, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 3)

    # Solve
    U = K \ F
    scattersysvec!(dchi, U[:])
    result = dchi.values[nl, 3][1]
    # @info "Solution: $(result), $(round(result/analyt_sol*100, digits = 4))%"

    max_n = Float64[]
    min_n = Float64[]
    max_en = Float64[]
    min_en = Float64[]

    for nc = 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
        push!(max_n, maximum(fld.values))
        push!(min_n, minimum(fld.values))
        fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
        push!(max_en, maximum(fld.values))
        push!(min_en, minimum(fld.values))
    end

    max_m = Float64[]
    min_m = Float64[]
    max_em = Float64[]
    min_em = Float64[]

    for nc = 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
        push!(max_m, maximum(fld.values))
        push!(min_m, minimum(fld.values))
        fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
        push!(max_em, maximum(fld.values))
        push!(min_em, minimum(fld.values))
    end

    max_q = Float64[]
    min_q = Float64[]
    max_eq = Float64[]
    min_eq = Float64[]

    for nc = 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
        push!(max_q, maximum(fld.values))
        push!(min_q, minimum(fld.values))
        fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
        push!(max_eq, maximum(fld.values))
        push!(min_eq, minimum(fld.values))
    end

    if visualize
        scalars = []
        for nc = 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            push!(scalars, ("m$nc", fld.values))
            fld =
                elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            push!(scalars, ("em$nc", fld.values))
        end
        vtkwrite(
            "homo-$(n)-m.vtu",
            fens,
            fes;
            scalars = scalars,
            vectors = [("u", dchi.values[:, 1:3])],
        )
        scalars = []
        for nc = 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            push!(scalars, ("n$nc", fld.values))
            fld = elemfieldfromintegpoints(
                femm,
                geom0,
                dchi,
                :membrane,
                nc,
                outputcsys = ocsys,
            )
            push!(scalars, ("en$nc", fld.values))
        end
        vtkwrite(
            "homo-$(n)-n.vtu",
            fens,
            fes;
            scalars = scalars,
            vectors = [("u", dchi.values[:, 1:3])],
        )
        scalars = []
        for nc = 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("q$nc", fld.values))
            fld =
                elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("eq$nc", fld.values))
        end
        vtkwrite(
            "homo-o-$(n)-q.vtu",
            fens,
            fes;
            scalars = scalars,
            vectors = [("u", dchi.values[:, 1:3])],
        )
    end
    return max_n,
    min_n,
    max_en,
    min_en,
    max_m,
    min_m,
    max_em,
    min_em,
    max_q,
    min_q,
    max_eq,
    min_eq
end

function test_composite(visualize = true)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    n = 8
    tolerance = R / n / 10
    fens, fes = T3block(40 / 360 * 2 * pi, L / 2, n, n, :b)
    fens.xyz = xyz3(fens)
    for i = 1:count(fens)
        a = fens.xyz[i, 1]
        y = fens.xyz[i, 2]
        fens.xyz[i, :] .= (R * sin(a), y, R * (cos(a) - 1))
    end

    ocsys = CSys(3, 3, cylindrical!)
    mater = CM.lamina_material(E, nu)
    plies = [CM.Ply("ply_1", mater, thickness, 0)]

    layup = CM.CompositeLayup("sco-lo", plies, ocsys)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), layup)
    femm.drilling_stiffness_scale = 1.0e-4
    femm.mult_el_size = 5 / 12
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # rigid diaphragm
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1, 3, 5]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf L / 2 L / 2 -Inf Inf], inflate = tolerance)
    for i in [2, 4, 6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1, 5, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Midpoint of the free edge
    nl = selectnode(
        fens;
        box = Float64[sin(40 / 360 * 2 * pi) * 25 sin(40 / 360 * 2 * pi) * 25 L / 2 L / 2 -Inf Inf],
        inflate = tolerance,
    )
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity([0, 0, -90, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 3)

    # Solve
    U = K \ F
    scattersysvec!(dchi, U[:])
    result = dchi.values[nl, 3][1]
    # @info "Solution: $(result), $(round(result/analyt_sol*100, digits = 4))%"

    max_n = Float64[]
    min_n = Float64[]
    max_en = Float64[]
    min_en = Float64[]

    for nc = 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
        push!(max_n, maximum(fld.values))
        push!(min_n, minimum(fld.values))
        fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
        push!(max_en, maximum(fld.values))
        push!(min_en, minimum(fld.values))
    end

    max_m = Float64[]
    min_m = Float64[]
    max_em = Float64[]
    min_em = Float64[]

    for nc = 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
        push!(max_m, maximum(fld.values))
        push!(min_m, minimum(fld.values))
        fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
        push!(max_em, maximum(fld.values))
        push!(min_em, minimum(fld.values))
    end

    max_q = Float64[]
    min_q = Float64[]
    max_eq = Float64[]
    min_eq = Float64[]

    for nc = 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
        push!(max_q, maximum(fld.values))
        push!(min_q, minimum(fld.values))
        fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
        push!(max_eq, maximum(fld.values))
        push!(min_eq, minimum(fld.values))
    end

    if visualize
        scalars = []
        for nc = 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            push!(scalars, ("m$nc", fld.values))
            fld =
                elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            push!(scalars, ("em$nc", fld.values))
        end
        vtkwrite(
            "comp-$(n)-m.vtu",
            fens,
            fes;
            scalars = scalars,
            vectors = [("u", dchi.values[:, 1:3])],
        )
        scalars = []
        for nc = 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            push!(scalars, ("n$nc", fld.values))
            fld = elemfieldfromintegpoints(
                femm,
                geom0,
                dchi,
                :membrane,
                nc,
                outputcsys = ocsys,
            )
            push!(scalars, ("en$nc", fld.values))
        end
        vtkwrite(
            "comp-$(n)-n.vtu",
            fens,
            fes;
            scalars = scalars,
            vectors = [("u", dchi.values[:, 1:3])],
        )
        scalars = []
        for nc = 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("q$nc", fld.values))
            fld =
                elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("eq$nc", fld.values))
        end
        vtkwrite(
            "comp-o-$(n)-q.vtu",
            fens,
            fes;
            scalars = scalars,
            vectors = [("u", dchi.values[:, 1:3])],
        )
    end
    return max_n,
    min_n,
    max_en,
    min_en,
    max_m,
    min_m,
    max_em,
    min_em,
    max_q,
    min_q,
    max_eq,
    min_eq
end

end # module
using LinearAlgebra
using Test
using .mcompostat0
rh = mcompostat0.test_homogeneous()
rc = mcompostat0.test_composite()

for i = 1:length(rc)
    @test norm(rc[i] - rh[i]) / norm(rh[i]) < 1.0e56
end
