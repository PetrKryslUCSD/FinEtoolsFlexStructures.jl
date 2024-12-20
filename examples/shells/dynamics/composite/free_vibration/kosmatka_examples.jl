"""
Vibration analysis of anisotropic layered square plates

Reference values: ? [Hz]

"""
module kosmatka_examples

using Arpack
using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, richextrapol
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function _execute_steel_triangle(aspect, n, refndom, visualize)
    formul = FEMMShellT3FFCompModule
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
   
    fens = FENodeSet([0 0; a 0; 0 a])
    fes = FESetT3([(1, 2, 3)])
    for r in 1:n
        fens, fes = T3refine(fens, fes)
    end
    fens.xyz = xyz3(fens)
    
    vtkwrite("kosmatka_n=$n.vtu", fens, fes)

    mater = CM.lamina_material(Material...)
    plies = CM.Ply[
        CM.Ply("ply_0", mater, thickness / 4, 0),
        CM.Ply("ply_90", mater, thickness / 4, 90),
        CM.Ply("ply_90", mater, thickness / 4, 90),
        CM.Ply("ply_0", mater, thickness / 4, 0),
    ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("kosmatka-0/90/90/0", plies, mcsys)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), layup)
    femm.transv_shear_formulation = formul.__TRANSV_SHEAR_FORMULATION_AVERAGE_K
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


function _execute_clamped_rectangular(n, orientation, reffs, visualize)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    # Geometry from figure 12
    a = 9.0 * phun("in")
    b = 3.0 * phun("in")
    h = 0.1 * phun("in")
    # Material from Table III
    rho = 0.05781 * phun("lbm/in^3")
    E1 = 22.43e6 * phun("psi")
    E2 = 1.45e6 * phun("psi")
    G12 = 0.917e6 * phun("psi")
    G13 = G12
    nu12 = 0.25
    nu22 = 0.4
    G23 = E2 / (2*(1+nu22))
    Material = (rho, E1, E2, nu12, G12, G13, G23)

     
    # Report
    @info "Mesh: $n elements per side"

    tolerance = a/n/1000
   
    fens, fes = T3block(a, b, 3*n, n, orientation); # Mesh
   
    fens.xyz = xyz3(fens)
    
    vtkwrite("kosmatka_n=$n.vtu", fens, fes)

    mater = CM.lamina_material(Material...)
    angles = [22.5, -22.5, 0, 0, -22.5, 0, 0, 22.5, 0, 0]
    # symmetric
    angles = vcat([a for a in angles], [a for a in reverse(angles)])
    
    plies = CM.Ply[]
    for (j, a) in enumerate(angles)
        push!(plies, CM.Ply("ply_$j", mater, h / length(angles), a))
    end
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("kosmatka-rect-plate", plies, mcsys)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), h), layup)
    # femm.transv_shear_formulation = formul.__TRANSV_SHEAR_FORMULATION_AVERAGE_K
    associategeometry! = formul.associategeometry!
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
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf 0 0], inflate = tolerance)
    for i in [1, 2, 3, 4, 5, 6]
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
    neigvs = length(reffs)
    d, v, nconv = eigs(K_ff+OmegaShift*M_ff, M_ff; nev=neigvs, which=:SM, explicittransform=:none)
    d[:] = d .- OmegaShift;
    oms = real(sqrt.(complex(d)))
    fs = oms ./(2*pi)
    @info "Angular: $oms [rad]"
    @info "Frequencies: $fs [Hz]"
        
    # Visualization
    if visualize
        for j in 1:neigvs
            U = v[:, j]
            scattersysvec!(dchi, (a * 4) / maximum(abs.(U)) .* U)
            update_rotation_field!(Rfield0, dchi)
            plots = cat(plot_space_box([[0 0 0]; [a a a]]),
                #plot_nodes(fens),
                plot_midsurface(fens, fes; x=geom0.values, u=dchi.values[:, 1:3], R=Rfield0.values);
                dims=1)
            pl = render(plots)
            sleep(1.5)
        end
    end
    return fs
end

function test_mesh()
    @info "Kosmatka plate vibration"
    reffs = [70.0, 271.0, 433.8, 888.3, 1196.5, 1370.9, 1730.3, 1737.1]
    for n in [4]
            _execute_clamped_rectangular(n, :a, reffs, false)
            _execute_clamped_rectangular(n, :b, reffs, false)
        end
    return true
end


function test_convergence()
    @info "Kosmatka plate vibration: approach to the limit"
    reffs = [70.0, 271.0, 433.8, 888.3, 1196.5, 1370.9, 1730.3, 1737.1]
    ns = 2 * [4, 8, 16]
    results = []
    for n in ns
        fs = _execute_clamped_rectangular(n, :b, reffs, false)
        push!(results, fs)
    end
    @info "Extrapolation:"
    q1, q2, q3 = results[1:3]
    for j in 1:length(reffs)
        e = richextrapol([q1[j], q2[j], q3[j]], 1.0 ./ Float64.(ns))
        @info "$j: $e"
    end
    @info "Reference (not in the limit): \n $reffs"
    return true
end

function allrun()
    println("#####################################################")
    println("# test_mesh ")
    test_mesh()
    println("#####################################################")
    println("# test_convergence ")
    test_convergence()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
