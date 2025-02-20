"""
Vibration analysis of antisymmetric angle ply square plates

Layup [th/-th] or [th/-th/th/-th/th/-th/th/-th], where th = 45.
Aspect Ratio variable.
Material: E1/ E2 = 40, G23 = 0*5E2, G12 = G13 = 0*6E, v = 0.25.

Reference values: Table 4
Journal of Sound and Vibration (1985) 98(2), 157-170 
STABILITY AND VIBRATION OF ISOTROPIC, ORTHOTROPIC 
AND LAMINATED PLATES ACCORDING TO A 
HIGHER-ORDER SHEAR DEFORMATION THEORY 
J. N. REDDY AND N. D. PHAN 

===========================================================
Can't get sufficiently close to these reference values!!!
===========================================================

"""
module reddy_phan_1985_examples

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

nondimensionalised_frequency = function(m, a, thickness, om) 
    # @show m, a, thickness
    rho = m[1]
    E2 = m[3]
    om * a^2 / thickness * sqrt(rho/E2) 
end

function _execute_angle(nplies, aspect, n, reference, visualize)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    angle = 45
    a = 1.0 * phun("m")
    thickness = a / aspect
    @assert nplies == 2  || nplies == 8
    rho = 1200 * phun("KG/M^3")
    E2 = 3 * phun("GPa")
    E1 = E2 * 40
    G12 = E2 * 0.6
    G13 = G12
    G23 = E2 * 0.5
    nu12 = 0.25
    material = (rho, E1, E2, nu12, G12, G13, G23)
       
    # Report
    @info "Mesh: $n elements per side"

    tolerance = a/n/1000
   
    fens,fes = Q4block(a, a, n, n); # Mesh
    fens.xyz = xyz3(fens)
    fens, fes = Q4toT3(fens, fes)
    
    vtkwrite("reddy_phan_1985_nplies=$(nplies)-aspect=$(aspect)-n=$n.vtu", fens, fes)

    plies = CM.Ply[
        CM.Ply("ply_$k", CM.lamina_material(material...), thickness / nplies, (-1)^(k+1)*angle)
        for k in 1:nplies
    ]
    
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("laminate_$nplies", plies, mcsys)
    
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
    # simple support
    l1 = connectednodes(meshboundary(fes))
    for i in  [1, 2, 3, ] 
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
    OmegaShift = (0.1*2*pi)^2
    neigvs = 1
    d, v, nconv = eigs(K_ff+OmegaShift*M_ff, M_ff; nev=neigvs, which=:SM, explicittransform=:none)
    d[:] = d .- OmegaShift;
    oms = real(sqrt.(complex(d)))
    fs = oms ./(2*pi)
    ndoms = nondimensionalised_frequency(material, a, thickness, oms)
    @info "Nondimensional frequency: $(ndoms)  vs Reference: $(reference) "
        
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
            push!(vectors, ("mode_$ev-$(round(ndoms[ev]; sigdigits=4))", deepcopy(dchi.values[:, 1:3])))
            # vtkwrite("reddy_1979-mode-$(ev).vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        end
        vtkwrite("reddy_1979-a=$angle-np=$(nplies)-modes.vtu", fens, fes; vectors=vectors)
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




function test_angle()
    @info "Antisymmetric angle-ply laminated plate vibration"
    @info "Reference from Reddy, Chao 1985.  "
    reference = Dict(
        "nplies=2, a/h=5" => 10.335,
        "nplies=2, a/h=10" => 13.044,
        "nplies=2, a/h=20" => 14.179,
        "nplies=2, a/h=100" => 14.618,
        "nplies=8, a/h=5" => 12.892,
        "nplies=8, a/h=10" => 19.289, 
        "nplies=8, a/h=20" => 23.259,
        "nplies=8, a/h=100" => 25.176
    )
    for nplies in [2, 8]
        @info "Number of plies: $nplies"
        for aspect in [5, 10, 20, 100]
            @info "Aspect: $aspect"
            for n in [10, 20, 40,]
                _execute_angle(nplies, aspect, n, reference["nplies=$nplies, a/h=$aspect"], true)
            end
        end
    end
    
    return true
end



function allrun()
    println("#####################################################")
    println("# test_angle ")
    test_angle()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
