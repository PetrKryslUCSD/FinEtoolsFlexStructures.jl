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
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellDSG3Module
# using FinEtoolsFlexStructures.FEMMShellT3Module
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function test_dsg3(n = 8, visualize = true)
    E = 200e3*phun("MPa")
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    thickness = 0.05*phun("m");
    L = 10.0*phun("m");
    
    formul = FEMMShellDSG3Module
       # Report
    @info "FV12 free vibration, t/L=$(thickness/L), formulation=$(formul)"
    @info "Mesh: $n elements per side"

    tolerance = L/n/1000
    fens, fes = T3block(L,L,n,n);
    fens.xyz[:, 1] .-= L/2
    fens.xyz[:, 2] .-= L/2
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.FEMMShellDSG3(IntegDomain(fes, TriRule(1), thickness), mater)
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
    @info "Frequencies: $(round.(fs[7:10], digits=4))"

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
    if !visualize
        return true
    end
    for ev in 1:10
        U = evecs[:, ev]
        scattersysvec!(dchi, (0.5*L)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
            plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots; title="$(ev)")
    end
    return true
end


function test_q4sri(n = 8, visualize = true)
    E = 200e3*phun("MPa")
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    thickness = 0.05*phun("m");
    L = 10.0*phun("m");
    
    formul = FEMMShellQ4SRIModule
       # Report
    @info "FV12 free vibration, t/L=$(thickness/L), formulation=$(formul)"
    @info "Mesh: $n elements per side"

    tolerance = L/n/1000
    fens, fes = Q4block(L,L,n,n);
    fens.xyz[:, 1] .-= L/2
    fens.xyz[:, 2] .-= L/2
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    
    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.FEMMShellQ4SRI(IntegDomain(fes, GaussRule(2, 2), thickness), mater)
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
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    M = mass(femm, geom0, dchi);
    # @show sum(sum(M, dims = 1))/3

    # Solve
    OmegaShift = (0.5*2*pi)^2
    neigvs = 24
    evals, evecs, nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    # @show nconv
    evals[:] = evals .- OmegaShift;
    fs = real(sqrt.(complex(evals)))/(2*pi)
    @info "Frequencies: $(round.(fs[7:10], digits=4))"

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
    if !visualize
        return true
    end
    for ev in 1:10
        U = evecs[:, ev]
        scattersysvec!(dchi, (0.5*L)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
            plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots; title="$(ev)")
    end
    return true
end

function single_q4sri()
    E = 200e3*phun("MPa")
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    thickness = 0.05*phun("m");
    L = 10.0*phun("m");
    
    lambdaij = sqrt.([13.49, 19.79, 24.43])
    fij(k) = lambdaij[k]^2/(2*pi*L^2)*sqrt(E*thickness^2/12/rho/(1-nu^2))

    @show fij.(1:3)
    # Mesh
    n = 8
    tolerance = L/n/1000
    fens, fes = Q4block(L,L,n,n);
    fens.xyz[:, 1] .-= L/2
    fens.xyz[:, 2] .-= L/2
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    
    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = FEMMShellQ4SRIModule.FEMMShellQ4SRI(IntegDomain(fes, GaussRule(2, 2), thickness), mater)
    stiffness = FEMMShellQ4SRIModule.stiffness
    mass = FEMMShellQ4SRIModule.mass

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
    # @show sum(sum(M, dims = 1))/3

    # Solve
    OmegaShift = (0.5*2*pi)^2
    neigvs = 10
    evals, evecs, nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    # @show nconv
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
    if !visualize
        return true
    end
    for ev in 7:neigvs
        U = evecs[:, ev]
        scattersysvec!(dchi, (0.5*L)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
            plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots; title="$(ev)")
    end
    return true
end

function test_dsg3_convergence()
    for n in [2, 4, 8, 16, 32, 64]
        test_dsg3(n, false)
    end
    return true
end

function test_q4sri_convergence()
    for n in [2, 4, 8, 16, 32, 64]
        test_q4sri(n, false)
    end
    return true
end

function allrun()
    println("#####################################################")
    println("# test_dsg3 ")
    test_dsg3()
    println("#####################################################")
    println("# test_q4sri  ")
    test_q4sri()
    println("#####################################################")
    println("# test_dsg3_convergence  ")
    test_dsg3_convergence()
    println("#####################################################")
    println("# test_q4sri_convergence  ")
    test_q4sri_convergence()
    return true
end # function allrun

end # module

using .FV12_vibration_examples
FV12_vibration_examples.test_dsg3_convergence()
FV12_vibration_examples.test_q4sri_convergence()