"""
Vibration analysis of free-floating thin plate

Free-vibration problem is solved for a homogeneous free-floating
(unsupported) square plate. This is the NAFEMS Benchmark, Test No. FV12.

The plate is discretized with shell elements. Because no displacements are
prevented, the structure has six rigid body modes (six zero vibration
frequencies).

The mass matrix is computed without consideration for the rotations. Hence, the
frequencies are underestimated.    

The nonzero benchmark frequencies are (in hertz): 

1.622, 2.360, 2.922, 4.190, 4.190,  7.356, 7.356, 7.668.
"""
module FV12_vibration_examples

using LinearAlgebra
using Arpack
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3, local_frame!
using FinEtoolsFlexStructures.FEMMShellDSG3Module: FEMMShellDSG3, stiffness, mass
# using FinEtoolsFlexStructures.FEMMShellT3Module: FEMMShellT3, stiffness, mass
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function single_dsg3()
    E = 200e3*phun("MPa")
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    thickness = 0.05*phun("m");
    L = 10.0*phun("m");
    @show rho*L*L*thickness

# Mesh
    n = 64
    tolerance = L/n/1000
    fens, fes = T3block(L,L,n,n);
    fens.xyz[:, 1] .-= L/2
    fens.xyz[:, 2] .-= L/2
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellDSG3(IntegDomain(fes, TriRule(1), thickness), mater)

# Construct the requisite fields, geometry and displacement
# Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

# Apply EBC's
    l1 = collect(1:count(fens))
    # for i in [6]
    #     setebc!(dchi, l1, true, i)
    # end
    applyebc!(dchi)
    numberdofs!(dchi);

# Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    M = mass(femm, geom0, dchi);
    @show sum(sum(M, dims = 1))/3

# Solve
    OmegaShift = (0.5*2*pi)^2
    neigvs = 24
    evals, evecs, nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    @show nconv
    evals[:] = evals .- OmegaShift;
    fs = real(sqrt.(complex(evals)))/(2*pi)
    @show fs

    # sol = eigen(Matrix(K+OmegaShift*M), Matrix(M))
    # evals = sol.values
    # evals[:] = evals .- OmegaShift;
    # fs = real(sqrt.(complex(evals)))/(2*pi)
    #     @show fs[9]
    #     evecs = sol.vectors

    # for index_ in 1:neigvs
    #     om = evals[index_]
    #     v = evecs[:, index_]
    #     @show index_, norm(K*v-om*M*v)
    #     @show (v'*K*v)/om
    #     @show (v'*M*v)
    # end

    # for n in dchi.dofnums[:, 6]
    #     # @show sqrt(K[n, n] / M[n, n]) / 2/pi
    #     if n != 0
    #         @show K[n, :]
    #     end
    # end
        
# Visualization
    # for ev in 1:10
    #     U = evecs[:, ev]
    #     @show ev
    #     @show argmax(abs.(U))
    #     scattersysvec!(dchi, (0.5*L)/maximum(abs.(U)).*U)
    #     update_rotation_field!(Rfield0, dchi)
    #     plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
    #         plot_nodes(fens),
    #         plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    #         dims = 1)
    #     pl = render(plots; title="$(ev)")
    # end
end


function allrun()
    println("#####################################################")
    println("# single_dsg3 ")
    single_dsg3()
    return true
end # function allrun

end # module

using .FV12_vibration_examples
FV12_vibration_examples.allrun()