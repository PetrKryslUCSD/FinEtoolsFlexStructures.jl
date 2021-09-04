"""
Vibration analysis of simply-supported thick  (solid) plate

Free-vibration problem is solved for a homogeneous square plate,
simply-supported on the circumference.
This is the NAFEMS Benchmark, Test No. FV52.

The plate is discretized with shell elements. 
The in-plane displacements are prevented.

The mass matrix is computed without consideration for the rotations. Hence, the
frequencies are underestimated.    

The nonzero benchmark frequencies are (in hertz): 
45.897  109.44  109.44  167.89  
Note: these are only for the out of plane 
mode shapes. The in-plane modes are not included.
"""
module FV52_vibration_examples

using Arpack
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3, local_frame!
using FinEtoolsFlexStructures.FEMMShellT3DSGOModule
using FinEtoolsFlexStructures.FEMMShellT3DSGICModule
using FinEtoolsFlexStructures.FEMMShellT3DSGModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function _execute_dsg_model(formul, n, visualize)
    E = 200e3*phun("MPa")
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    thickness = 1.0*phun("m");
    L = 10.0*phun("m");

    # Report
    @info "Mesh: $n elements per side"

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
    # simple support
    l1 = selectnode(fens; box = Float64[L/2 L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [3,4,6]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[-L/2 -L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [3,4,6]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance)
    for i in [3,5,6]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[-Inf Inf -L/2 -L/2 -Inf Inf], inflate = tolerance)
    for i in [3,5,6]
        setebc!(dchi, l1, true, i)
    end
    # in-plane, rotations
    l1 = collect(1:count(fens))
    for i in [1, 2, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    M = mass(femm, geom0, dchi);

    # Solve
    OmegaShift = 0.1*2*pi
    neigvs = 4
    d, v, nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    d[:] = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    @show fs
        
    # Visualization
    U = v[:, 4]
    scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
end


function test_convergence()
    formul = FEMMShellT3DSGModule
    @info "FV12 free vibration, formulation=$(formul)"
    for n in [2, 4, 8, 16, 32, 64]
        _execute_dsg_model(formul, n, false)
    end
    return true
end

end # module

using .FV52_vibration_examples
m = FV52_vibration_examples
m.test_convergence()