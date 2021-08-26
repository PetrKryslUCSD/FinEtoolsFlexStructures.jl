"""
Vibration analysis of free-floating thin plate

Free-vibration problem is solved for a homogeneous free-floating
(unsupported) square plate. This is the NAFEMS Benchmark, Test No. FV12.

The plate is discretized with shell elements. Because no displacements are
prevented, the structure has six rigid body modes (six zero vibration
frequencies).

The nonzero benchmark frequencies are (in hertz): 

1.622, 2.360, 2.922, 4.190, 4.190,  7.356, 7.356, 7.668.
"""
module FV12_vibration_dsg3_verification

using Test
using LinearAlgebra
using Arpack
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellDSG3Module
using FinEtoolsFlexStructures.FEMMShellDSG3IModule
using FinEtoolsFlexStructures.FEMMShellDSG3IFModule
# using FinEtoolsFlexStructures.FEMMShellT3Module
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function test_dsg3if(args...)
    return _execute_dsg_model(FEMMShellDSG3IFModule, args...)
end

function test_dsg3i(args...)
    return _execute_dsg_model(FEMMShellDSG3IModule, args...)
end

function test_dsg3(args...)
    return _execute_dsg_model(FEMMShellDSG3Module, args...)
end

function test_csdsg3(args...)
    return _execute_dsg_model(FEMMShellCSDSG3Module, args...)
end

function _execute_dsg_model(formul, n = 8, visualize = true)
    E = 200e3*phun("MPa")
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    thickness = 0.05*phun("m");
    L = 10.0*phun("m");
    
    # Report
    # @info "Mesh: $n elements per side"

    tolerance = L/n/1000
    fens, fes = T3block(L,L,n,n);
    fens.xyz[:, 1] .-= L/2
    fens.xyz[:, 2] .-= L/2
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    associate = formul.associategeometry!
    stiffness = formul.stiffness
    mass = formul.mass

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
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    M = mass(femm, geom0, dchi);
    # @show sum(sum(M, dims = 1))/3,.

    # Solve
    OmegaShift = (0.5*2*pi)^2
    neigvs = 24
    evals, evecs, nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    # @show nconv
    evals[:] = evals .- OmegaShift;
    fs = real(sqrt.(complex(evals)))/(2*pi)
    # @info "Frequencies: $(round.(fs[7:10], digits=4))"

    # Visualization
    # if !visualize
    #     return true
    # end
    # for ev in 1:10
    #     U = evecs[:, ev]
    #     scattersysvec!(dchi, (0.5*L)/maximum(abs.(U)).*U)
    #     update_rotation_field!(Rfield0, dchi)
    #     plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
    #         plot_nodes(fens),
    #         plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    #         dims = 1)
    #     pl = render(plots; title="$(ev)")
    # end
    return fs
end


function test_convergence(t)
    @info "FV12 free vibration, formulation=$(t)"
    for n in [2, 4, 8, 16, 32, 64]
        t(n, false)
    end
    return true
end

# 1.622, 2.360, 2.922, 4.190, 4.190,  7.356, 7.356, 7.668.
reffs = [0.0, 0.0, 2.604869127850317e-7, 4.698288861559094e-6, 6.749716652051837e-6, 9.829373450823236e-6, 1.572130183778014, 2.2424585076387427, 2.8079394352847316, 3.883763676656034, 4.039123204140305, 6.787320617260535, 6.920636670319986, 7.127888889722697] 
fs = test_dsg3if()
for j in eachindex(reffs)
    @test isapprox(fs[j], reffs[j], atol = 1.0e-8, rtol = 1.0e-6)
end

end # module
