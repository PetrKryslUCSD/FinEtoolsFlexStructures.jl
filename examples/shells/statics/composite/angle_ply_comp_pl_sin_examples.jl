"""
A two variable reﬁ ned plate theory for laminated composite plates
Seung-Eock Kim a, Huu-Tai Thai a, Jaehong Lee b,* Composite Structures 89 (2009) 197–205

Example 2: [45/-45]s angle-ply composite plate under sinusoidal pressure loading.
The plate is simply supported on all edges (type SS-2) 
and subjected to a sinusoidal pressure loading.

Table 3
Nondimensionalized deﬂ ections of simply supported two-layer (45/ 45) square and
rectangular laminates under sinusoidal transverse load
a/h Source                                w
                   Square plate (a = b) Rectangular plate (b = 3a)
4  Ren [20]              1.4471                   3.9653
   HSDT                  1.0203                   3.1560
   FSDT                  1.1576                   3.3814
   RPT1                  1.0203                   3.0971
   RPT2                  1.0633                   3.1503
10 Ren [20]              0.6427                   2.3953
   HSDT                  0.5581                   2.2439
   FSDT                  0.5773                   2.2784
   RPT1                  0.5581                   2.2325
   RPT2                  0.5628                   2.2395
100 Ren [20]             0.4685                   2.0686
   HSDT                  0.4676                   2.0671
   FSDT                  0.4678                   2.0674
   RPT1                  0.4676                   2.0670
   RPT2                  0.4677                   2.0670
   CLPT                  0.4667                   2.0653
"""
module angle_ply_comp_pl_sin_examples
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.FEMMShellQ4RSCompModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function normalized_deflection(w, a, h, q, E2)
    return 100 * w * E2 * h^3 / (q * a^4)
end

const REFERENCE_SOLUTIONS = Dict(
        (3.0, 4.0) => 3.3814,
        (3.0, 10.0) => 2.2784,
        (3.0, 100.0) => 2.0674,
        (1.0, 4.0) => 1.1576,
        (1.0, 10.0) => 0.5773,
        (1.0, 100.0) => 0.4678
    )

function test_t3ff(b_over_a = 3.0, a_over_h = 4.0, nx = 8, ny = 8, visualize = true)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule

    a = 250 * phun("mm") # note the changed dimension
    b = b_over_a * a
    # high modulus orthotropic graphite-epoxy composite material
    E1 = 100 * phun("GPa")
    E2 = E1 / 40
    G12 = G13 = E2 * 0.5
    nu12 = 0.25
    G23 = E2 * 0.6
    q = 10.0 * phun("kilo*Pa")
    h = a / a_over_h
    
    tolerance = a / nx / 100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[
        CM.Ply("ply_+45", mater, h / 2, 45),
        CM.Ply("ply_-45", mater, h / 2, -45)
    ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("ren", plies, mcsys)

    fens, fes = T3block(a, b, nx, ny)
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)

    femm = formul.make(IntegDomain(fes, TriRule(1), h), layup)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    lcenter = selectnode(fens; box = Float64[a/2 a/2 b/2 b/2 0 0], inflate = tolerance)
    # Apply EBC's
    # The boundary conditions are SS-2: refer to the paper.
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf 0 0], inflate = tolerance)
    for i in [1, 4]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[a a -Inf Inf 0 0], inflate = tolerance)
    for i in [1, 4]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 0 0], inflate = tolerance)
    for i in [2, 5]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[-Inf Inf b b 0 0], inflate = tolerance)
    for i in [2, 5]
        setebc!(dchi, l1, true, i)
    end
    # Simple support
    l1 = connectednodes(meshboundary(fes))
    for i in [3]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, massem, geom0, u0, Rfield0, dchi)

    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64, 6, (forceout, XYZ, tangents, feid, qpid) -> begin
        forceout .= [0.0, 0.0, q * sin(XYZ[1] / a * π) * sin(XYZ[2] / b * π), 0.0, 0.0, 0.0]
        return forceout
    end)
    F = distribloads(lfemm, vassem, geom0, dchi, fi, 2)

    # Solve
    U = K \ F
    scattersysvec!(dchi, U[:])
    @info "b/a = $b_over_a, a/h = $a_over_h, nx = $nx, ny = $ny"
    @info "Center deflection = $(normalized_deflection(dchi.values[lcenter, 3][1], a, h, q, E2))"
    @info "Reference deflection = $(REFERENCE_SOLUTIONS[(b_over_a, a_over_h)])"

    # Visualization
    if visualize
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = layup.csys)
            push!(scalars, ("m$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = layup.csys)
            push!(scalars, ("em$nc", fld.values))
        end
        vtkwrite("angle_sin_t3ff-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = layup.csys)
            push!(scalars, ("n$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = layup.csys)
            push!(scalars, ("en$nc", fld.values))
        end
        vtkwrite("angle_sin_t3ff-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = layup.csys)
            push!(scalars, ("q$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = layup.csys)
            push!(scalars, ("eq$nc", fld.values))
        end
        vtkwrite("angle_sin_t3ff-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

    end
    true
end

function test_q4rs(b_over_a = 3.0, a_over_h = 4.0, nx = 8, ny = 8, visualize = true)
    formul = FEMMShellQ4RSCompModule
    CM = CompositeLayupModule

    a = 250 * phun("mm") # note the changed dimension
    b = b_over_a * a
    # high modulus orthotropic graphite-epoxy composite material
    E1 = 100 * phun("GPa")
    E2 = E1 / 40
    G12 = G13 = E2 * 0.5
    nu12 = 0.25
    G23 = E2 * 0.6
    q = 10.0 * phun("kilo*Pa")
    h = a / a_over_h
    
    tolerance = a / nx / 100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[
        CM.Ply("ply_+45", mater, h / 2, 45),
        CM.Ply("ply_-45", mater, h / 2, -45)
    ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("ren", plies, mcsys)

    fens, fes = Q4block(a, b, nx, ny)
    fens.xyz = xyz3(fens)
    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)

    femm = formul.make(IntegDomain(fes, GaussRule(2, 2), h), layup)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    lcenter = selectnode(fens; box = Float64[a/2 a/2 b/2 b/2 0 0], inflate = tolerance)
    # Apply EBC's
    # The boundary conditions are SS-2: refer to the paper.
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf 0 0], inflate = tolerance)
    for i in [1, 4]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[a a -Inf Inf 0 0], inflate = tolerance)
    for i in [1, 4]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 0 0], inflate = tolerance)
    for i in [2, 5]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[-Inf Inf b b 0 0], inflate = tolerance)
    for i in [2, 5]
        setebc!(dchi, l1, true, i)
    end
    # Simple support
    l1 = connectednodes(meshboundary(fes))
    for i in [3]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, massem, geom0, u0, Rfield0, dchi)

    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(Float64, 6, (forceout, XYZ, tangents, feid, qpid) -> begin
        forceout .= [0.0, 0.0, q * sin(XYZ[1] / a * π) * sin(XYZ[2] / b * π), 0.0, 0.0, 0.0]
        return forceout
    end)
    F = distribloads(lfemm, vassem, geom0, dchi, fi, 2)

    # Solve
    U = K \ F
    scattersysvec!(dchi, U[:])
    @info "b/a = $b_over_a, a/h = $a_over_h, nx = $nx, ny = $ny"
    @info "Center deflection = $(normalized_deflection(dchi.values[lcenter, 3][1], a, h, q, E2))"
    @info "Reference deflection = $(REFERENCE_SOLUTIONS[(b_over_a, a_over_h)])"

    # Visualization
    if visualize
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = layup.csys)
            push!(scalars, ("m$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = layup.csys)
            push!(scalars, ("em$nc", fld.values))
        end
        vtkwrite("angle_sin_t3ff-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = layup.csys)
            push!(scalars, ("n$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = layup.csys)
            push!(scalars, ("en$nc", fld.values))
        end
        vtkwrite("angle_sin_t3ff-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = layup.csys)
            push!(scalars, ("q$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = layup.csys)
            push!(scalars, ("eq$nc", fld.values))
        end
        vtkwrite("angle_sin_t3ff-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

    end
    true
end

function allrun()
    N = 64
    println("#####################################################")
    println("# test_t3ff ")
    test_t3ff(1.0, 4.0, N, N, false)
    test_t3ff(1.0, 10.0, N, N, false)
    test_t3ff(1.0, 100.0, N, N, false)
    test_t3ff(3.0, 4.0, N, N, false)
    test_t3ff(3.0, 10.0, N, N, false)
    test_t3ff(3.0, 100.0, N, N, false)
    println("#####################################################")
    println("# test_q4rs ")
    test_q4rs(1.0, 4.0, N, N, false)
    test_q4rs(1.0, 10.0, N, N, false)
    test_q4rs(1.0, 100.0, N, N, false)
    test_q4rs(3.0, 4.0, N, N, false)
    test_q4rs(3.0, 10.0, N, N, false)
    test_q4rs(3.0, 100.0, N, N, false)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
