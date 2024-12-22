"""
Vibration analysis of anisotropic layered square plates

Reference values: ? [Hz]
Table 5.6-14 
Effect of ply angle (``\\theta``) and number of layers (n) 
on dimensionless fundamental frequency, ZJ, of 
antisymmetric angle-play (``\\theta``/ ``-\\theta``/``\\theta``/ ... /``-\\theta``) 
square plates {a/h = 10).

SOMETHING WRONG: The reference values are not matching with the results. 
"""
module angle_ply_examples

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

function _execute(angle, n, nplies, refndom, visualize)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    # Material from section 2.1
    rho = 1600*phun("KG/M^3")
    E2 = 3*phun("GPa")
    E1 = E2 * 40
    G12 = E2 * 0.6
    G13 = G12
    G23 = E2 * 0.5
    nu12 = 0.25
    Material = (rho, E1, E2, nu12, G12, G13, G23)
    aspect = 10
       
    a = 1.0 * phun("m")
    thickness = a / aspect
    nondimensionalised_frequency = function(om) 
        rho = Material[1]
        E2 = Material[3]
        om * a^2 / thickness * sqrt(rho/E2) 
    end

    # Report
    @info "Angle: $angle"
    @info "Number of plies: $nplies"
    @info "Mesh: $n elements per side"

    tolerance = a/n/1000
   
    fens,fes = Q4block(a, a, n, n); # Mesh
    fens.xyz = xyz3(fens)
    fens, fes = Q4toT3(fens, fes)
    
    vtkwrite("angle_ply_n=$n.vtu", fens, fes)

    sign = -1
    plies = CM.Ply[
        CM.Ply("ply_$j", CM.lamina_material(Material...), thickness / nplies, sign^(j+1)*angle)
         for j = 1:nplies
    ]
    # @show [p.angle for p in plies]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("angle_ply_$nplies", plies, mcsys)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), layup)
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
    # simple support along two sides
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in  [3, ] 
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[a a -Inf Inf -Inf Inf], inflate = tolerance)
    for i in  [3, ] 
        setebc!(dchi, l1, true, i)
    end
    # one side is clamped
    l1 = selectnode(fens; box = Float64[-Inf Inf a a -Inf Inf], inflate = tolerance)
    for i in  [1, 2, 3, 4, 5, 6] 
        setebc!(dchi, l1, true, i)
    end
    # the remaining side is free
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    M = mass(femm, geom0, dchi);

    K_ff = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:ff]
    M_ff = matrix_blocked(M, nfreedofs(dchi), nfreedofs(dchi))[:ff]

    # Solve
    OmegaShift = (0.1*2*pi)^2
    neigvs = 9
    d, v, nconv = eigs(K_ff+OmegaShift*M_ff, M_ff; nev=neigvs, which=:SM, explicittransform=:none)
    d[:] = d .- OmegaShift;
    oms = real(sqrt.(complex(d)))
    fs = oms ./(2*pi)
    # @info "Frequencies: $fs [Hz]"
    @info "Nondimensional angular frequency  $(nondimensionalised_frequency.(oms)) vs ref $(refndom)"
        
    # Visualization
    if visualize
        # U = v[:, 4]
        # scattersysvec!(dchi, (a*4)/maximum(abs.(U)).*U)
        # update_rotation_field!(Rfield0, dchi)
        # plots = cat(plot_space_box([[0 0 0]; [a a a]]),
        #             #plot_nodes(fens),
        #             plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
        #             dims = 1)
        # pl = render(plots)
        vectors = []
        for ev in 1:neigvs
            U = v[:, ev]
            scattersysvec!(dchi, 1.0/maximum(abs.(U)).*U)
            push!(vectors, ("mode_$ev", deepcopy(dchi.values[:, 1:3])))
            # vtkwrite("angle_ply-mode-$(ev).vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        end
            vtkwrite("angle_ply-a=$angle-np=$(nplies)-modes.vtu", fens, fes; vectors = vectors)
            # update_rotation_field!(Rfield0, dchi)
            # plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
            #     plot_nodes(fens),
            #     plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            #     dims = 1)
            # pl = render(plots; title="$(ev)")
        # end
    end
    nothing
end


function test_convergence()
    @info "angle_ply plate vibration"
    @info "Reference from Miravete, Reddy Table 5.6-14"
    reference = (
    (angle = 30, nplies = 2, fundamental = 12.68),
    (angle = 30, nplies = 10, fundamental = 18.51),
    (angle = 45, nplies = 2, fundamental = 13.04),
    (angle = 45, nplies = 10, fundamental = 19.38),
    (angle = 60, nplies = 2, fundamental = 12.68),
    (angle = 60, nplies = 10, fundamental = 18.51),
    )
    n = 10
    for ref in reference
            _execute(ref.angle, n, ref.nplies, ref.fundamental, true)
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
