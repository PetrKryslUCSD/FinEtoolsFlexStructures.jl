"""
Vibration analysis of anisotropic layered square plates

Reference values: ? [Hz]
From: Table 3, page 76
Dynamic stiffness formulation for composite Mindlin plates for exact modal 
analysis of structures. Part II: Results and applications
M. Boscolo, J.R. Banerjee
Computers and Structures 96-97 (2012) 74–83
"""
module boscolo_examples

using Arpack
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function _execute(formul, aspect, n, refndom, visualize)
    CM = CompositeLayupModule
    # Material from section 2.1
    rho = 2000*phun("KG/M^3")
    E2 = 3*phun("GPa")
    E1 = E2 * 40
    G12 = E2 * 0.6
    G13 = G12
    G23 = E2 * 0.5
    nu12 = 0.25
    Material = (rho, E1, E2, nu12, G12, G13, G23)

       
    a = 1.0 * phun("m")
    thickness = a / aspect
    nondimensionalised_frequency = function(om) 
        rho = Material[1]
        E2 = Material[3]
        om * a^2 / thickness * sqrt(rho/E2) 
    end

    # Report
    @info "Aspect: $aspect"
    @info "Mesh: $n elements per side"

    tolerance = a/n/1000
   
    fens,fes = Q4block(a, a, n, n); # Mesh
    fens.xyz = xyz3(fens)
    fens, fes = Q4toT3(fens, fes)
    
    vtkwrite("boscolo_n=$n.vtu", fens, fes)

    mater = CM.lamina_material(Material...)
    plies = CM.Ply[
        CM.Ply("ply_0", mater, thickness / 3, 0),
        CM.Ply("ply_90", mater, thickness / 3, 90),
        CM.Ply("ply_0", mater, thickness / 3, 0),
    ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("boscolo_3-ply", plies, mcsys)

    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    # femm.transv_shear_formulation = formul.__TRANSV_SHEAR_FORMULATION_AVERAGE_K
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
    l1 = connectednodes(meshboundary(fes))
    for i in [1, 2, 3, ]
        setebc!(dchi, l1, true, i)
    end
    
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    M = mass(femm, geom0, dchi);

    K_ff = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:ff]
    M_ff = matrix_blocked(M, nfreedofs(dchi), nfreedofs(dchi))[:ff]

    # Solve
    OmegaShift = 0.0*2*pi
    neigvs = 1
    d, v, nconv = eigs(K_ff+OmegaShift*M_ff, M_ff; nev=neigvs, which=:SM, explicittransform=:none)
    d[:] = d .- OmegaShift;
    oms = real(sqrt.(complex(d)))
    fs = oms ./(2*pi)
    @info "Frequencies: $fs [Hz]"
    @info "Nondimensional angular frequency  $(nondimensionalised_frequency.(oms)) vs ref $(refndom)"
        
    # Visualization
    if visualize
        U = v[:, 4]
        scattersysvec!(dchi, (a*4)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 0]; [a a a]]),
                    #plot_nodes(fens),
                    plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
                    dims = 1)
        pl = render(plots)
    end
    nothing
end


function test_convergence()
    formul = FEMMShellT3FFModule
    @info "Boscolo-Banerjee plate vibration"
    refndoms = [5.500, 9.395, 10.854, 15.143, 17.660, 18.071]
    for (aspect, refndom) in zip([2, 4, 5, 10, 20, 25], refndoms)
        for n in [50]
            _execute(formul, aspect, n, refndom, false)
        end
    end
    return true
end


function allrun()
    println("#####################################################")
    println("# test_convergence ")
    test_convergence()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
