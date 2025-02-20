module mbuckling1
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMCorotBeamModule
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
mass = FEMMCorotBeamModule.mass
geostiffness = FEMMCorotBeamModule.geostiffness
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
function test()
    # Parameters:
    E = 1000000.0
    nu = 0.3
    L = 30.0 # Length of the beam
    b = 0.5 # width
    h = 4.0 # height
    magn = -0.2056167583560 * 1e4 / L^2
    magn_scale = 1000
    neigvs = 4

    # Reference frequencies
    reffs = [48.5475, 124.1839]

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [-1.0, 0.0, 0.0])

    # Select the number of elements per leg.
    n = 32
    members = []
    push!(members, frame_member([0 0 0; 0 0 L], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000)

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    l1 = selectnode(fens; box = Float64[0 0 0 0 0 0], tolerance = L / 10000)
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the global discrete system
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)
    tipn = selectnode(fens; box = Float64[0 0 0 0 L L], tolerance = L / 10000)
    loadbdry = FESetP1(reshape(tipn, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity([0, -magn * b * h / magn_scale, 0, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 3)

    # Solve the static problem
    solve_blocked!(dchi, K, F)
    # @show  dchi.values[tipn, :]

    # Compute the geometric stiffness
    u0.values = dchi.values[:, 1:3]
    update_rotation_field!(Rfield0, dchi)
    Kg = geostiffness(femm, geom0, u0, Rfield0, dchi)

    K_ff = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:ff]
    Kg_ff = matrix_blocked(Kg, nfreedofs(dchi), nfreedofs(dchi))[:ff]

    # Solve the eigenvalue problem
    d, v, nev, nconv = eigs(-Kg_ff, K_ff; nev = neigvs, which = :LM)
    fs = 1.0 ./ (d) ./ magn_scale
    # println("Buckling factors: $fs [ND]")
    # println("Reference: $reffs [ND]")
    @test norm(reffs - abs.(fs[2:3])) ./ norm(reffs) <= 1.0e-3

    # # Visualize vibration modes
    # scattersysvec!(dchi, 10*v[:, 1])
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_space_box([[-L/2 -L/2 0]; [L/2 L/2 1.1*L]]),
    #  plot_nodes(fens),
    #     plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    #     dims = 1)
    # render(plots)
    # @show p.plot.data
    true
end
end
using .mbuckling1
mbuckling1.test()
