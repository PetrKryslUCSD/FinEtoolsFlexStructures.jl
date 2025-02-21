"""
Vibration analysis of antisymmetric angle ply square plates

Layup [th/-th/th/-th], where th = 0, 30, 45.
Aspect Ratio 10.
Two different materials of the plies.

Decent agreement can be obtained with material II, but not with material I.

Reference values: Table 2 for all angles; also Table 6 (th=45), Table 7 (th=30)
Journal of Sound and Vibration (1979) 66(4), 565-576
FREE VIBRATION OF ANTISYMMETRIC,
ANGLE -PLY LAMINATED PLATES INCLUDING
TRANSVERSE SHEAR DEFORMATION BY THE
FINITE ELEMENT METHOD
J.N.REDDY

"""
module reddy_1979_examples

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

function _execute(mid, angle, n, reference, visualize)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    nplies = 4
    a = 1.0 * phun("m")
    thickness = a / 10
    materials = []
    # Material I
    rho = 1200 * phun("KG/M^3")
    E2 = 3 * phun("GPa")
    E1 = E2 * 40
    G12 = E2 * 0.6
    G13 = E2 * 0.6
    G23 = E2 * 0.5
    nu12 = 0.25
    push!(materials, (rho, E1, E2, nu12, G12, G13, G23))
    # Material II
    rho = 1200 * phun("KG/M^3")
    E2 = 3 * phun("GPa")
    E1 = E2 * 25
    G12 = E2 * 0.5
    G13 = E2 * 0.2
    G23 = E2 * 0.5
    nu12 = 0.25
    push!(materials, (rho, E1, E2, nu12, G12, G13, G23))
       
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
        CM.Ply("ply_$k", CM.lamina_material(materials[mid]...), thickness / nplies, (-1)^(k+1)*angle)
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
    # simple support
    l1 = connectednodes(meshboundary(fes))
    for i in  [3, ] 
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens, box=Float64[0, 0, -Inf, Inf, 0, 0], inflate=tolerance)
    setebc!(dchi, l1, true, 1)
    l1 = selectnode(fens, box=Float64[a, a, -Inf, Inf, 0, 0], inflate=tolerance)
    setebc!(dchi, l1, true, 1)
    l1 = selectnode(fens, box=Float64[-Inf, Inf, 0, 0, 0, 0], inflate=tolerance)
    setebc!(dchi, l1, true, 2)
    l1 = selectnode(fens, box=Float64[-Inf, Inf, a, a, 0, 0], inflate=tolerance)
    setebc!(dchi, l1, true, 2)
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
    ndoms = nondimensionalised_frequency(materials[mid], a, thickness, oms)
    @info "Nondimensional frequency ($bestfit): $(ndoms[bestfit])  vs Reference: $(reference) "
        
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
        vtkwrite("reddy_1979-mid=$mid-a=$angle-np=$(nplies)-modes.vtu", fens, fes; vectors=vectors)
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




function test()
    @info "Antisymmetric angle-ply laminated plate vibration"
    @info "Number of plies: 4"
    @info "Reference from Reddy 1979. Table 2 "
    reference = Dict(
        "1/0" => 14.193, "1/30" => 17.689, "1/45" => 18.609,
        "2/0" => 10.358, "2/30" => 12.739, "2/45" => 13.303,
    )
    for mid in [1, 2]
        @info "Material: $mid"
        for angle in [0, 30, 45, ]
            @info "Angle: $angle"
            for n in [10, 20, 40, ]
                _execute(mid, angle, n, reference["$mid/$angle"], true)
            end
        end
    end
    
    return true
end



function allrun()
    println("#####################################################")
    println("# test ")
    test()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
