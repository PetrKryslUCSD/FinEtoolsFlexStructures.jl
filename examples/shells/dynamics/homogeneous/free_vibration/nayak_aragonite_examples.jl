"""
Free vibration analysis of composite sandwich plates based on Reddy’s
higher-order theory
A.K. Nayak, S.S.J. Moy, R.A. Shenoi
Composites: Part B 33 (2002) 505–519
Section 4.3, Free vibration analysis of a simply supported square
plate (orthotropic)

From Table 2, nondimensional frequencies
(1,1)    0.0474
(1,2)    0.1033
(2,1)    0.1188
(2,2)    0.1694
(1,3)    0.1888
(3,1)    0.2180
(2,3)    0.2475
(3,2)    0.2624
(1,4)    0.2969
"""
module nayak_aragonite_examples

using Arpack
using LinearAlgebra: norm, Transpose, mul!, I, Symmetric
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test(nplies = 10, axes = (1, 2, 3))
@info "nplies = $nplies, axes = $axes"

    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule

    ax = ay = 100 * phun("mm")
    nx = ny = 90
    # aragonite crystals
    E1 = 143.52 * phun("GPa")
    E2 = 75.38 * phun("GPa")
    G12 = 42.03 * phun("GPa")
    G13 = 25.56 * phun("GPa")
    nu12 = 0.44
    nu23 = 0.23
    G23 = 42.65 * phun("GPa")
    rho = 1500.0 * phun("kg/m^3") # only a guess
    C11 = 159.85 * phun("GPa")
    thickness = ax / 10
    tolerance = ax / nx / 100
    

    mater = CM.lamina_material(rho, E1, E2, nu12, G12, G13, G23)
    plies = [CM.Ply("ply_$i", mater, thickness / nplies, 0) for i = 1:nplies]

    mcsys = CM.cartesian_csys(axes)
    layup = CM.CompositeLayup("Nayak 4.3", plies, mcsys)

    fens, fes = T3block(ax, ay, nx, ny)
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)

    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), layup)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # Pin one of the corners
    l1 = selectnode(fens; box = Float64[ax ax 0 0 -Inf Inf], inflate = tolerance)
    for i in [1, 2, 3]
        setebc!(dchi, l1, true, i)
    end
    # Roller at the other
    l1 = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], inflate = tolerance)
    for i in [2]
        setebc!(dchi, l1, true, i)
    end
    # Simple support
    l1 = connectednodes(meshboundary(fes))
    for i in [1, 2, 3]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, geom0, u0, Rfield0, dchi)
    M = formul.mass(femm, geom0, dchi)
    K_ff = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:ff]
    M_ff = matrix_blocked(M, nfreedofs(dchi), nfreedofs(dchi))[:ff]

    # Solve
    neigvs = 9
    d, v, nconv = eigs(
        Symmetric(K_ff),
        Symmetric(M_ff);
        nev = neigvs,
        which = :SM,
        explicittransform = :none,
    )
    @test nconv == neigvs
    fs = real(sqrt.(complex(d))) / (2 * pi)
    @info "Frequencies: $(fs)"
    @info "ND Frequencies: $(thickness * sqrt(rho / C11) .* (2 * pi * fs))"
    @info "Relative result: $(thickness * sqrt(rho / C11) .* (2 * pi * fs) ./ [0.0474, 0.1033, 0.1188, 0.1694, 0.1888, 0.2180, 0.2475, 0.2624, 0.2969])"

    vectors = []
    for i = 1:neigvs
        scattersysvec!(dchi, v[:, i])
        push!(vectors, ("mode_$i", deepcopy(dchi.values[:, 1:3])))
    end
    vtkwrite("plate-modes.vtu", fens, fes; vectors = vectors)

    true
end

function allrun()
    println("#####################################################")
    println("# test_convergence ")
    test()
    test(3, (2, -1, 3))
    test(5, (-1, -2, 3))
    test(4, (-2, 1, 3))
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
