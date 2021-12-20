# Cook membrane problem, plane stress.
module cook_membrane_examples

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using Examples.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

function _execute(n = 8, visualize = true)
    E = 1.0;
    nu = 1.0/3;
    width = 48.0; 
    height = 44.0; 
    thickness  = 1.0;
    free_height  = 16;
    mid_edge  = [48.0, 52.0];# Location of tracked  deflection
    magn=1/free_height;# Magnitude of applied load
    convutip=23.97;
    formul = FEMMShellT3FFModule
    
    @info "Mesh: $n elements per side"
    tolerance = thickness/n/100
    fens, fes = T3block(width,height,n,n);
    for k in 1:count(fens)
        x = fens.xyz[k,1]
        y = fens.xyz[k,2] + (x/width)*(height -fens.xyz[k,2]/height*(height-free_height));
        fens.xyz[k, :] .= (x, y)
    end
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped edge
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    bfes = meshboundary(fes)
    nl = selectelem(fens, bfes; box = Float64[width width -Inf Inf -Inf Inf], tolerance = tolerance)
    loadbdry = subset(bfes, nl)
    lfemm = FEMMBase(IntegDomain(loadbdry, GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, magn, 0, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);

    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    
    l1 = selectnode(fens; box = [mid_edge[1] mid_edge[1] mid_edge[2] mid_edge[2] -Inf Inf], inflate = tolerance)
    resultpercent = dchi.values[l1, 2][1]/convutip*100
    @info "Solution: $(round(resultpercent, digits = 4))%"

    # Generate a graphical display of resultants
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc)
        push!(scalars, ("n$nc", fld.values))
    end
    vtkwrite("cook_membrane-n=$n-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    
    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (height/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 0]; [width height 0]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return resultpercent
end

function test_convergence()
    
    @info "Cook trapezoidal membrane"

       
    for n in [2, 4, 8, 16, 32, 64, 128]
        _execute(n, false)
    end
    return true
end

end # module

using .cook_membrane_examples
m = cook_membrane_examples
m.test_convergence()