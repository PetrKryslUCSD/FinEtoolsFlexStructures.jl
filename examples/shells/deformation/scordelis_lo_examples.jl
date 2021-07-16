# Pinched cylinder with diagphram supports and concentrated force
module scordelis_lo_examples

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3, local_frame!
using FinEtoolsFlexStructures.FEMMShellDSG3Module: FEMMShellDSG3, stiffness
# using FinEtoolsFlexStructures.FEMMShellT3Module: FEMMShellT3, stiffness
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function single_dsg3()
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    analyt_sol=-0.3024;
    # Parameters:
    E=4.32e8;
    nu=0.0;
    thickness = 0.25; # geometrical dimensions are in feet

# Mesh
    R = 25.0;
    L = 50.0;
    n = 8;
    tolerance = R/n/1000
    fens, fes = T3block(40/360*2*pi,L/2,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*cos(a))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellDSG3(IntegDomain(fes, TriRule(1), thickness), mater)

# Construct the requisite fields, geometry and displacement
# Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

# Apply EBC's
# rigid diaphragm
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,3,5]
        setebc!(dchi, l1, true, i)
    end
# plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
# plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

# Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

# Midpoint of the free edge
    nl = selectnode(fens; box = Float64[sin(40/360*2*pi)*25 sin(40/360*2*pi)*25 L/2 L/2 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt[0, 0, -90, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
# Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    @show dchi.values[nl, 3],  dchi.values[nl, 3]/analyt_sol*100

# Visualization
    scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
end


function allrun()
    println("#####################################################")
    println("# single_dsg3 ")
    single_dsg3()
    return true
end # function allrun

end # module

using .scordelis_lo_examples
scordelis_lo_examples.allrun()