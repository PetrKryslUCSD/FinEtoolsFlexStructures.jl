"""
Vibration analysis of antisymmetric angle ply square plates

Layup [th/-th]n, where th = 30, 45, and n = 1, 2, 5.
Aspect Ratio 10.

The boundary condition is assumed 

- simply supported (SS-2). Very good agreement.
- free floating plate. Atrocious disagreement.

Reference values: 
Joumal of Sound and Vibration (1988) 122(2),377-388

FREE VIBRATION OF ANTISYMMETRIC ANGLE-PLY 
LAMINATED PLATES INCLUDING VARIOUS 
BOUNDARY CONDITIONS

A. A. KHDEIR
"""
module khdeir_1988_examples

using Arpack
using LinearAlgebra
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
    rho = m[1]
    E2 = m[3]
    om * a^2 / thickness * sqrt(rho/E2) 
end

function _execute_ss(nplies, angle, n, reference, visualize)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    a = 1.0 * phun("m")
    thickness = a / 10
    materials = []
    # Material from from second paragraph of section 4
    rho = 1200 * phun("KG/M^3")
    E2 = 3 * phun("GPa")
    E1 = E2 * 40
    G12 = E2 * 0.6
    G13 = E2 * 0.6
    G23 = E2 * 0.5
    nu12 = 0.25
    material = (rho, E1, E2, nu12, G12, G13, G23)
       
    # Report
    @info "Mesh: $n elements per side"

    tolerance = a/n/1000
   
    fens,fes = Q4block(a, a, n, n); # Mesh
    fens.xyz = xyz3(fens)
    fens, fes = Q4toT3(fens, fes)
    
    appshape = deepcopy(fens.xyz)
    for j in axes(appshape, 1)
        appshape[j, 3] = sin(pi*appshape[j, 1]/a)*sin(pi*appshape[j, 2]/a)
        appshape[j, 1] = appshape[j, 2] = 0
    end
    vtkwrite("reddy_1979_n=$n.vtu", fens, fes, vectors = [("appshape", appshape)])

    plies = CM.Ply[
        CM.Ply("ply_$k", CM.lamina_material(material...), thickness / nplies, (-1)^(k+1)*angle)
        for k in 1:nplies
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
    # simple support (SS-2, meaning that displacements normal to the edge are suppressed, as is rotation about normal)
    l1 = connectednodes(meshboundary(fes))
    for i in  [3, 6] 
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens, box=Float64[0, 0, -Inf, Inf, 0, 0], inflate=tolerance)
    setebc!(dchi, l1, true, 1)
    setebc!(dchi, l1, true, 4)
    l1 = selectnode(fens, box=Float64[a, a, -Inf, Inf, 0, 0], inflate=tolerance)
    setebc!(dchi, l1, true, 1)
    setebc!(dchi, l1, true, 4)
    l1 = selectnode(fens, box=Float64[-Inf, Inf, 0, 0, 0, 0], inflate=tolerance)
    setebc!(dchi, l1, true, 2)
    setebc!(dchi, l1, true, 5)
    l1 = selectnode(fens, box=Float64[-Inf, Inf, a, a, 0, 0], inflate=tolerance)
    setebc!(dchi, l1, true, 2)
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
    OmegaShift = (0.1*2*pi)^2
    neigvs = 10
    d, v, nconv = eigs(K_ff+OmegaShift*M_ff, M_ff; nev=neigvs, which=:SM, explicittransform=:none)
    d[:] = d .- OmegaShift;
    oms = real(sqrt.(complex(d)))
    proj = 0.0
    bestfit = 0
    for ev in 1:neigvs
        U = v[:, ev]
        scattersysvec!(dchi, U)
        p = dot(dchi.values[:, 1:3][:], appshape[:])
        if abs(p) > proj
            proj = abs(p)
            bestfit = ev
        end
    end
    ndoms = nondimensionalised_frequency(material, a, thickness, oms)
    @info "Nondim. freq. ($bestfit): $(ndoms[bestfit])  vs ref: $(reference): $(ndoms[bestfit]/reference*100) %"
        
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
        vtkwrite("khdeir_1988-SS-a=$angle-np=$(nplies)-modes.vtu", fens, fes; vectors=vectors)
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




function test_ss()
    @info "Antisymmetric angle-ply laminated plate vibration"
    
    @info "Reference from KHDEIR 1988. Table 1, SS boundary condition "
    reference = Dict(
        "2/30" => 12.68,
        "4/30" => 17.63,
        "10/30" => 18.51,
        "2/45" => 13.04,
        "4/45" => 18.46,
        "10/45" => 19.38
    )
    for nplies in [2, 4, 10]
        @info "\nnplies: $nplies"
        for angle in [30, 45, ]
            @info "Angle: $angle"
            for n in [10, 20, 40, ]
                _execute_ss(nplies, angle, n, reference["$nplies/$angle"], true)
            end
        end
    end
    
    return true
end


function _execute_ff(nplies, angle, n, reference, visualize)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    a = 1.0 * phun("m")
    thickness = a / 10
    materials = []
    # Material from from second paragraph of section 4
    rho = 1200 * phun("KG/M^3")
    E2 = 3 * phun("GPa")
    E1 = E2 * 40
    G12 = E2 * 0.6
    G13 = E2 * 0.6
    G23 = E2 * 0.5
    nu12 = 0.25
    material = (rho, E1, E2, nu12, G12, G13, G23)
       
    # Report
    @info "Mesh: $n elements per side"

    tolerance = a/n/1000
   
    fens,fes = Q4block(a, a, n, n); # Mesh
    fens.xyz = xyz3(fens)
    fens, fes = Q4toT3(fens, fes)
    
    appshape = deepcopy(fens.xyz)
    for j in axes(appshape, 1)
        appshape[j, 3] = sin(pi*appshape[j, 1]/a)*sin(pi*appshape[j, 2]/a)
        appshape[j, 1] = appshape[j, 2] = 0
    end
    vtkwrite("reddy_1979_n=$n.vtu", fens, fes, vectors = [("appshape", appshape)])

    plies = CM.Ply[
        CM.Ply("ply_$k", CM.lamina_material(material...), thickness / nplies, (-1)^(k+1)*angle)
        for k in 1:nplies
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
    # FF - free floating plate
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    M = mass(femm, geom0, dchi);

    K_ff = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:ff]
    M_ff = matrix_blocked(M, nfreedofs(dchi), nfreedofs(dchi))[:ff]

    # Solve
    OmegaShift = (0.1*2*pi)^2
    neigvs = 10
    d, v, nconv = eigs(K_ff+OmegaShift*M_ff, M_ff; nev=neigvs, which=:SM, explicittransform=:none)
    d[:] = d .- OmegaShift;
    oms = real(sqrt.(complex(d)))
    proj = 0.0
    ndoms = nondimensionalised_frequency(material, a, thickness, oms)
    @info "Nondim. freq.: $(ndoms[7])  vs ref: $(reference): $(ndoms[7]/reference*100) %"
        
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
        vtkwrite("khdeir_1988-FF-a=$angle-np=$(nplies)-modes.vtu", fens, fes; vectors=vectors)
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

function test_ff()
    @info "Antisymmetric angle-ply laminated plate vibration"
    
    @info "Reference from KHDEIR 1988. Table 1, FF boundary condition "
    reference = Dict(
        "2/30" => 6.95,
        "4/30" => 9.61,
        "10/30" => 10.11,
        "2/45" => 4.76,
        "4/45" => 6.26,
        "10/45" => 6.57
    )
    for nplies in [2, 4, 10]
        @info "\nnplies: $nplies"
        for angle in [30, 45, ]
            @info "Angle: $angle"
            for n in [10, 20, 40, ]
                _execute_ff(nplies, angle, n, reference["$nplies/$angle"], true)
            end
        end
    end
    
    return true
end



function allrun()
    println("#####################################################")
    println("# test_ss ")
    test_ss()
    println("#####################################################")
    println("# test_ff ")
    test_ff()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
