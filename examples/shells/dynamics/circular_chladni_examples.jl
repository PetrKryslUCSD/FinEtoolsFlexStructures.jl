"""
# Spherical cap, transient vibration.

The structure is loaded with a suddenly applied pressure.
The downward deflection at the centre is monitored.
One quarter of the entire shell is modeled.
"""
module circular_chladni_examples

using LinearAlgebra
using SparseArrays
using Arpack
using FinEtools
using FinEtools.MeshUtilModule: gradedspace
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.AssemblyModule
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using SymRCM
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection, vtkwrite

function biastowardsends(a, b, N)
    sort(vcat(gradedspace(a, (a + b) / 2, Int(round(N / 2)), 2), gradedspace(b, (a + b) / 2, Int(round(N / 2)), 2)))
end

const E = 70.0e3 * phun("MPa")
const nu = 0.34
const rho = 2700 * phun("kg/m^3")
const R = 120.0 * phun("mm")
const h = 1.0 * phun("mm")
const wmag = 600 * phun("micro*m")
const o2shift = (20 * 2 * pi)^2
const f = [3.28451e+01, 6.52558e+01, 9.00651e+01, 2.08699e+02, 3.61328e+02, 3.67857e+02, 3.83570e+02, 5.62524e+02, 6.04128e+02, 7.96744e+02, 9.04541e+02, 1.04571e+03, 1.06972e+03, 1.09128e+03, 1.25119e+03, 1.37448e+03, 1.37887e+03, 1.44462e+03, 1.44883e+03, 1.64328e+03, 1.64494e+03, 1.72060e+03, 1.72569e+03, 1.90974e+03, 1.91348e+03, 2.05683e+03, 2.07907e+03, 2.10875e+03, 2.13029e+03, 2.42437e+03, 2.52235e+03, 2.53267e+03, 2.55700e+03]
const N = 16
_frequencies = eltype(f)[]
for i in 1:length(f)-1
    global _frequencies = vcat(_frequencies, biastowardsends(f[i], f[i+1], N))
end
const frequencies = sort(_frequencies)
const color = "red"
const visualize = true
const a0 = 0.1
const a1 = 0.0001

function test(n=16, visualize=true)
    color = "red"
    tolerance = R / n / 100

    tolerance = R / n / 1000
    fens, fes = Q4circlen(R, n)
    fens, fes = Q4toT3(fens, fes)
    renumb = (c) -> c[[1, 3, 2]] # make sure the mirrored copies have consistent orientation
    fens1, fes1 = mirrormesh(fens, fes, [0.0, -1.0, 0.0], [0.0, 0.0, 0.0], renumb = renumb)
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, tolerance)
    fes = cat(fesa[1], fesa[2])
    fens1, fes1 = mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], renumb = renumb)
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, tolerance)
    fes = cat(fesa[1], fesa[2])
    

    bfes = meshboundary(fes)
    if visualize
        vtkwrite("circular_plate_$n-mesh.vtu", fens, fes)
        vtkwrite("circular_plate_$n-boundary-mesh.vtu", fens, bfes)
    end

    @info "Mesh $(count(fens)) nodes, $(count(fes)) elements"

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, TriRule(1)))
    # C = connectionmatrix(femm, count(fens))
    # perm = symrcm(C)

    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), h), mater)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    lcenter = selectnode(fens; box=Float64[0 0 0 0 0 0], inflate=tolerance)
    for i in 1:6
        setebc!(dchi, lcenter, true, i)
    end
    setebc!(dchi, lcenter, true, 3, wmag)
    applyebc!(dchi)
    numberdofs!(dchi)

    @info "$(nfreedofs(dchi)) free dofs"

    # Assemble the system matrix
    FEMMShellT3FFModule.associategeometry!(femm, geom0)

    K = FEMMShellT3FFModule.stiffness(femm, geom0, u0, Rfield0, dchi)
    M = FEMMShellT3FFModule.mass(femm, geom0, dchi)

    K_ff = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:ff]
    K_fd = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:fd]

    M_ff = matrix_blocked(M, nfreedofs(dchi), nfreedofs(dchi))[:ff]
    M_fd = matrix_blocked(M, nfreedofs(dchi), nfreedofs(dchi))[:fd]

    neigvs = 50
    evals, evecs, nconv = eigs(K_ff + o2shift * M_ff, M_ff; nev=neigvs, which=:SM, explicittransform=:none)
    evals[:] = evals .- o2shift
    fs = real(sqrt.(complex(evals)))/(2*pi)
    @info "Frequencies: $(round.(fs, digits=4))"

    C_ff = a0*M_ff + a1*K_ff
    C_fd = a0*M_fd + a1*K_fd

    # Solve
    U_d = gathersysvec(dchi, :d)
    F_f = zeros(Complex{Float64}, nfreedofs(dchi))
    
    U_f = zeros(Complex{Float64}, nfreedofs(dchi), length(frequencies))
    for k in eachindex(frequencies)
        frequency = frequencies[k]
        omega = 2 * pi * frequency
        F_f .= - (K_fd + 1im * omega * C_fd - (-omega^2 * M_fd)) * U_d
        U_f[:, k] = (-omega^2 * M_ff + 1im * omega * C_ff + K_ff) \ F_f
    end

    if visualize
        @info "Dumping visualization"
        times = Float64[]
        vectors = []
        for i in axes(U_f, 2)
            scattersysvec!(dchi, real(U_f[:, i]))
            u = dchi.values[:, 1:3]
            u[:, 3] .-= wmag
            push!(vectors, ("U", deepcopy(u)))
            push!(times, frequencies[i])
        end
        vtkwritecollection("circular_plate_$n", fens, fes, times; vectors=vectors)
    end

    # if visualize
    #     for ev in 1:neigvs
    #         U = evecs[:, ev]
    #         scattersysvec!(dchi, 1.0/maximum(abs.(U)).*U)
    #         vtkwrite("circular_plate_$n-mode-$(ev).vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
    #         # update_rotation_field!(Rfield0, dchi)
    #         # plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
    #         #     plot_nodes(fens),
    #         #     plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    #         #     dims = 1)
    #         # pl = render(plots; title="$(ev)")
    #     end
    # end

    return true
end

function allrun(ns=16)
    println("#####################################################")
    println("# test ")
    test(ns)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing

