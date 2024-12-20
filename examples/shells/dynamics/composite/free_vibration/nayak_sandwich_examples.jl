"""
Free vibration analysis of composite sandwich plates based on Reddy’s
higher-order theory
A.K. Nayak, S.S.J. Moy, R.A. Shenoi
Composites: Part B 33 (2002) 505–519
Section 4.6, Free vibration analysis of a three layer sandwich plate 
with a soft core

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
module nayak_sandwich_examples

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

function _execute(n, visualize)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    # Material Face: aluminum
    rho = 2770*phun("KG/M^3")
    E2 = 68.9*phun("GPa")
    E1 = E2
    # Shear moduli not given in the article
    G12 = E2 / 2 / (1 + 0.34)
    G13 = G12
    G23 = G12
    nu12 = 0.3
    Material_Face = (rho, E1, E2, nu12, G12, G13, G23)
    # Material Core: aluminum honeycomb
    rho = 122*phun("KG/M^3")
    G12 = 0.134*phun("GPa")
    G13 = G12
    G23 = 0.052*phun("GPa")
    # Young's moduli, Poisson ratio not given
    E2 = 2*G12
    E1 = E2
    nu12 = 0.25
    Material_Core = (rho, E1, E2, nu12, G12, G13, G23)
    # Dimensions
    a = 1.83 * phun("m")
    b = 1.22 * phun("m")
    thickness_face = 0.406 * phun("mm")
    thickness_core = 6.4 * phun("mm")

    # Report
    @info "Mesh: $n elements per side"

    tolerance = a/n/1000
   
    fens,fes = Q4block(a, b, n, n); # Mesh
    fens.xyz = xyz3(fens)
    fens, fes = Q4toT3(fens, fes)
    
    vtkwrite("nayak_sandwich_n=$n.vtu", fens, fes)

    mater_face = CM.lamina_material(Material_Face...)
    mater_core = CM.lamina_material(Material_Core...)
    plies = CM.Ply[
        CM.Ply("ply_face", mater_face, thickness_face, 0),
        CM.Ply("ply_core", mater_core, thickness_core, 0),
        CM.Ply("ply_face", mater_face, thickness_face, 0),
    ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("nayak_sandwich_3-ply", plies, mcsys)

    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), CM.thickness(layup)), layup)
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
    neigvs = 10
    d, v, nconv = eigs(K_ff+OmegaShift*M_ff, M_ff; nev=neigvs, which=:SM, explicittransform=:none)
    d[:] = d .- OmegaShift;
    oms = real(sqrt.(complex(d)))
    fs = oms ./(2*pi)
    @info "Frequencies: $fs [Hz]"
            
    # Visualization
    if visualize
        L0 = max(a, b)
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
    @info "Nayak et al sandwich plate vibration"
        for n in [10, 100]
            _execute(n, true)
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
