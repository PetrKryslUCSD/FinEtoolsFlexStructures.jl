"""
Vibration analysis of cross-ply layered square plates

Layup [0/90/0] (thicknesses t/4, t/2, t/4).
Material II.

Reference values: Table 6 (page 163)
From: J.N. REDDY and W.C. CHAO 
A COMPARISON OF CLOSED-FORM AND FINITE-ELEMENT SOLUTIONS 
OF THICK LAMINATED ANISOTROPIC RECTANGULAR PLATES,
Nuclear Engineering and Design 64 (1981) 153-167 153 

"""
module reddy_chao_examples

using Arpack
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function _execute(aspect, n, reference, visualize)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    # # Material I
    # rho = 2250*phun("KG/M^3")
    # E2 = 3*phun("GPa")
    # E1 = E2 * 25
    # G12 = E2 * 0.5
    # G13 = G12
    # G23 = E2 * 0.2
    # nu12 = 0.25
    # Material_I = (rho, E1, E2, nu12, G12, G13, G23)
    # Material II
    rho = 2250*phun("KG/M^3")
    E2 = 3*phun("GPa")
    E1 = E2 * 40
    G12 = E2 * 0.6
    G13 = G12
    G23 = E2 * 0.5
    nu12 = 0.25
    Material_II = (rho, E1, E2, nu12, G12, G13, G23)

    Material = Material_II
   
    L0 = 0.3 * phun("m")
    thickness = L0 / aspect
    nondimensionalised_frequency = function(om) 
        rho = Material[1]
        E2 = Material[3]
        om * L0^2 / thickness * sqrt(rho/E2) 
    end

    # Report
    @info "Aspect: $aspect"
    @info "Mesh: $n elements per side"

    tolerance = L0/n/1000
   
    fens,fes = Q4block(L0, L0, n, n); # Mesh
    fens.xyz = xyz3(fens)
    fens, fes = Q4toT3(fens, fes)
    
    # vtkwrite("reddy_chao_n=$n.vtu", fens, fes)

    mater = CM.lamina_material(Material...)
    plies = CM.Ply[
        CM.Ply("ply_0", mater, thickness / 4, 0),
        CM.Ply("ply_90", mater, thickness / 2, 90),
        CM.Ply("ply_0", mater, thickness / 4, 0),
    ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("Reddy_Chao_3-ply", plies, mcsys)

    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), layup)
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
    for i in  [3, 6] 
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens, box=Float64[0, 0, -Inf, Inf, 0, 0], inflate=tolerance)
    setebc!(dchi, l1, true, 2)
    setebc!(dchi, l1, true, 4)
    l1 = selectnode(fens, box=Float64[L0, L0, -Inf, Inf, 0, 0], inflate=tolerance)
    setebc!(dchi, l1, true, 2)
    setebc!(dchi, l1, true, 4)
    l1 = selectnode(fens, box=Float64[-Inf, Inf, 0, 0, 0, 0], inflate=tolerance)
    setebc!(dchi, l1, true, 1)
    setebc!(dchi, l1, true, 5)
    l1 = selectnode(fens, box=Float64[-Inf, Inf, L0, L0, 0, 0], inflate=tolerance)
    setebc!(dchi, l1, true, 1)
    setebc!(dchi, l1, true, 5)
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
    neigvs = 10
    d, v, nconv = eigs(K_ff+OmegaShift*M_ff, M_ff; nev=neigvs, which=:SM, explicittransform=:none)
    d[:] = d .- OmegaShift;
    oms = real(sqrt.(complex(d)))
    fs = oms ./(2*pi)
    # @info "Frequencies: $fs [Hz]"
    ndom = nondimensionalised_frequency.(oms[1])
    @info "Nondim freq: $(ndom) ($(ndom/reference*100) %)"
        
    # Visualization
    if visualize
        U = v[:, 4]
        scattersysvec!(dchi, (L0*4)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 0]; [L0 L0 L0]]),
                    #plot_nodes(fens),
                    plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
                    dims = 1)
        pl = render(plots)
    end
    nothing
end


function test_convergence()
    reference = Dict(
        "a/h=5" => 10.861,
        "a/h=10" => 15.156,
        "a/h=25" => 18.089,
        "a/h=50" => 18.698,
    )
    @info "Reddy-Chao square plate vibration"
    for aspect in [5, 10, 25, 50]
        for n in [10, 20, 40]
            _execute(aspect, n, reference["a/h=$aspect"], false)
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
