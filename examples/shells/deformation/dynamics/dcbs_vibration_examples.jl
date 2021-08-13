"""
Vibration analysis of double-cell box structure
"""
module dbcs_vibration_examples

using Arpack
using FinEtools
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3, local_frame!
using FinEtoolsFlexStructures.FEMMShellDSG3Module
using FinEtoolsFlexStructures.FEMMShellDSG3IModule
using FinEtoolsFlexStructures.FEMMShellDSG3IFModule
# using FinEtoolsFlexStructures.FEMMShellT3Module: FEMMShellT3, stiffness
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

function _execute_dsg_model(formul, n, visualize)
    E = 207*phun("GPa")
    nu = 0.3;
    rho= 7850*phun("KG/M^3");
    L0 = 2.54 * phun("m")
    B0 = 1.27 * phun("m")
    thickness = 12.7*phun("mm");

    # Report
    @info "Mesh: $n elements per side"

    tolerance = B0/n/1000
    X = [0.0 0.0 0.0;
        B0 0.0 0.0;
        B0 B0 0.0;
        B0 2*B0 0.0;
        0.0 2*B0 0.0;
        0.0 B0 0.0;
    ]
    l2fens = FENodeSet(X);
    conn = [1 2; 6 3; 5 4; 2 3; 3 4; 5 6; 6 1]
    l2fes = FESetL2(conn);
    nLayers = 5
    ex = function (x, k)
        x[3] += k*L0/nLayers
        x
    end
    fens,fes = Q4extrudeL2(l2fens, l2fes, nLayers, ex); # Mesh
    fens, fes = Q4toT3(fens, fes)
    for r in 1:n
        fens, fes = T3refine(fens, fes)
    end
    @show count(fens)
    vtkwrite("dbcs.vtu", fens, fes)

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
    l1 = selectnode(fens; box = Float64[-Inf Inf -Inf Inf 0.0 0.0], inflate = tolerance)
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, l1, true, i)
    end
    
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    M = mass(femm, geom0, dchi);

    # Solve
    OmegaShift = 0.0*2*pi
    neigvs = 4
    d, v, nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    d[:] = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    @show fs
        
    # Visualization
    U = v[:, 4]
    scattersysvec!(dchi, (L0/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 0]; [L0 L0 L0]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
end


function test_convergence(t)
    @info "FV12 free vibration, formulation=$(t)"
    for n in [1, 2, 3, 4]
        t(n, false)
    end
    return true
end

function allrun()
    println("#####################################################")
    println("# test_dsg3 ")
    test_dsg3()
    return true
end # function allrun

end # module

using .dbcs_vibration_examples
m = dbcs_vibration_examples
m.test_convergence(m.test_dsg3if)