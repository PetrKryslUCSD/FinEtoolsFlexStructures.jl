"""
Vibration analysis of anisotropic layered square plates

Reference values: From the chapter "Inﬂuence of Shear Correction Factors on
Eigenfrequencies of Layered Plates", Friedrich Gruttmann and Werner Wagner
in H. Altenbach et al. (eds.), Advances in Mechanics of Materials and
Structural Analysis, Advanced Structured Materials 80,
https://doi.org/10.1007/978-3-319-70563-7_6 © Springer International Publishing
AG 2018



"""
module gruttman_examples

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

function _execute_30(n, reference, visualize)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    # Material from Section 4.1
    rho = 1500*phun("KG/M^3")
    E1 = 172369*phun("MPa")
    E2 = 6895*phun("MPa")
    G12 = 3447*phun("MPa")
    G13 = G12
    G23 = 1379*phun("MPa")
    nu12 = 0.25
    Material = (rho, E1, E2, nu12, G12, G13, G23)
    angle = 30
    nplies = 4
    a = 1.0 * phun("m")
    thickness = 50.0 * phun("mm")
       
    nondimensionalised_frequency = function(om) 
        rho = Material[1]
        E2 = Material[3]
        om * a^2 / thickness * sqrt(rho/E2) 
    end

    # Report
    @info "Mesh: $n elements per side"

    tolerance = a/n/1000
   
    fens,fes = Q4block(a, a, n, n); # Mesh
    fens.xyz = xyz3(fens)
    fens, fes = Q4toT3(fens, fes)
    
    vtkwrite("gruttman_n=$n.vtu", fens, fes)

    plies = CM.Ply[
        CM.Ply("ply_1", CM.lamina_material(Material...), thickness / nplies, angle)
        CM.Ply("ply_2", CM.lamina_material(Material...), thickness / nplies, -angle)
        CM.Ply("ply_3", CM.lamina_material(Material...), thickness / nplies, -angle)
        CM.Ply("ply_4", CM.lamina_material(Material...), thickness / nplies, angle)
    ]
    # @show [p.angle for p in plies]
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
    neigvs = length(reference)
    d, v, nconv = eigs(K_ff+OmegaShift*M_ff, M_ff; nev=neigvs, which=:SM, explicittransform=:none)
    d[:] = d .- OmegaShift;
    oms = real(sqrt.(complex(d)))
    fs = oms ./(2*pi)
    @info "Frequencies: $(round.(fs, digits=3)) \n     vs Reference: $(reference)  [Hz]\n         Relative: $(round.(fs ./ reference, digits=4))"
        
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
            # vtkwrite("gruttman-mode-$(ev).vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        end
            vtkwrite("gruttman-a=$angle-np=$(nplies)-modes.vtu", fens, fes; vectors = vectors)
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


function _execute_45(n, reference, visualize)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    # Material from Section 4.2.2
    # glue laminated timber
    rho = 530*phun("KG/M^3")
    E1 = 11600*phun("MPa")
    E2 = 390*phun("MPa")
    G12 = 720*phun("MPa")
    G13 = G12
    G23 = 100*phun("MPa")
    nu12 = 0.03
    Material = (rho, E1, E2, nu12, G12, G13, G23)
    angle = 45
    nplies = 8
    a = 4.8 * phun("m")
    thickness = 120 * phun("mm")
       
    nondimensionalised_frequency = function(om) 
        rho = Material[1]
        E2 = Material[3]
        om * a^2 / thickness * sqrt(rho/E2) 
    end

    # Report
    @info "Mesh: $n elements per side"

    tolerance = a/n/1000
   
    fens,fes = Q4block(a, a, n, n); # Mesh
    fens.xyz = xyz3(fens)
    fens, fes = Q4toT3(fens, fes)
    
    vtkwrite("gruttman_n=$n.vtu", fens, fes)

    plies = CM.Ply[
        CM.Ply("ply_1", CM.lamina_material(Material...), thickness / nplies, angle)
        CM.Ply("ply_2", CM.lamina_material(Material...), thickness / nplies, -angle)
        CM.Ply("ply_3", CM.lamina_material(Material...), thickness / nplies, angle)
        CM.Ply("ply_4", CM.lamina_material(Material...), thickness / nplies, -angle)
    ]
    plies = vcat(plies, reverse(plies))
    # @show [p.angle for p in plies]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("glulam_$nplies", plies, mcsys)
    
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
    neigvs = length(reference)
    d, v, nconv = eigs(K_ff+OmegaShift*M_ff, M_ff; nev=neigvs, which=:SM, explicittransform=:none)
    d[:] = d .- OmegaShift;
    oms = real(sqrt.(complex(d)))
    fs = oms ./(2*pi)
    @info "Frequencies: $(round.(fs, digits=3)) \n     vs Reference: $(reference)  [Hz]\n         Relative: $(round.(fs ./ reference, digits=4))"
        
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
            # vtkwrite("gruttman-mode-$(ev).vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        end
            vtkwrite("gruttman-a=$angle-np=$(nplies)-modes.vtu", fens, fes; vectors = vectors)
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

function test_30()
    @info "Symmetric angle-ply laminated plate vibration"
    @info "Angle: 30"
    @info "Number of plies: 4"
    @info "Reference from Gruttmann, Wagner. Table 2 for the [30, -30, -30, 30] layup."
    
    # The reference solution is not well defined: it is produced for a single mesh of
    # solid shell elements. Nevertheless, the differences are probably less than one
    # percent.
    reference = (
    253.20, 444.31, 669.01, 710.02, 861.31, 1020.29, 1148.31, 1187.15, 1352.61, 1352.36
    )
    ns = [10, 20, 40, 80, 160]
    for n in ns
        _execute_30(n, reference, true)
    end
    
    return true
end

function test_45()
    @info "Symmetric angle-ply glulam plate vibration"
    @info "Angle: 45"
    @info "Number of plies: 8"
    @info "Reference from Gruttmann, Wagner. Table 3 for the [45, -45, 45, -45]s layup."
    
    # The reference solution is not well defined: it is produced for a single mesh of
    # solid shell elements. Nevertheless, the differences are probably less than one
    # percent.
    reference = (
    14.28, 31.32, 34.77, 52.23, 60.93, 63.09, 77.86, 85.82, 97.74, 98.99
    )
    ns = [10, 20, 40, 80, 160]
    for n in ns
        _execute_45(n, reference, true)
    end
    
    return true
end


function allrun()
    println("#####################################################")
    println("# test_30 ")
    test_30()
    println("#####################################################")
    println("# test_45 ")
    test_45()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
