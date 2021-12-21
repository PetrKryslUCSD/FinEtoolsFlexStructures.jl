"""
Slit cylinder under torsional loads 
"""
module slit_cylinder_examples

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using Examples.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function _execute(n = 2, visualize = true)
    E = 3e6;
    nu = 0.3;
    thickness = 3.0;
    # analytical solution for the vertical deflection under the load
    analyt_sol=-1.82488e-5;
    formul = FEMMShellT3FFModule

    @info "Mesh: $n elements per side"

    # Mesh
    R = 300.0;
    L = 6*R;
    tolerance = min(R, L) / n  / 100
    
    tolerance = R/n/1000
    fens, fes = T3block(360.0, L/2, n, n);
    bfes = meshboundary(fes)
    al0 = selectelem(fens, bfes, box = [0.0 0.0 -Inf Inf], inflate = tolerance)
    al360 = selectelem(fens, bfes, box = [360.0 360.0 -Inf Inf], inflate = tolerance)
    ll0 = selectelem(fens, bfes, box = [-Inf Inf 0.0 0.0], inflate = tolerance)
    llL2 = selectelem(fens, bfes, box = [-Inf Inf L/2 L/2], inflate = tolerance)
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a/180*pi), y, R*cos(a/180*pi))
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    femm.drilling_stiffness_scale = 0.1
    associategeometry! = formul.associategeometry!
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # No applied supports
    numberdofs!(dchi);

    gsl = connectednodes(subset(bfes, vcat(ll0, llL2)))
    gsfes = FESetP1(reshape(gsl, length(gsl), 1))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    # loadbdry = subset(bfes, al0)
    # lfemm = FEMMBase(IntegDomain(loadbdry, GaussRule(1, 2)))
    # fi = ForceIntensity(FFlt[0, 1/4, 0, 0, 0, 0]);
    # F = distribloads(lfemm, geom0, dchi, fi, 1);
    # loadbdry = subset(bfes, al360)
    # lfemm = FEMMBase(IntegDomain(loadbdry, GaussRule(1, 2)))
    # fi = ForceIntensity(FFlt[0, -1/4, 0, 0, 0, 0]);
    # F += distribloads(lfemm, geom0, dchi, fi, 1);
    loadbdry = subset(bfes, ll0)
    lfemm = FEMMBase(IntegDomain(loadbdry, GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, 0, 0, 0, 1/4, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 1);
    loadbdry = subset(bfes, llL2)
    lfemm = FEMMBase(IntegDomain(loadbdry, GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, 0, 0, 0, -1/4, 0]);
    F += distribloads(lfemm, geom0, dchi, fi, 1);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    # @show n,  dchi.values[nl, 3]/analyt_sol*100

    # Visualization
    if visualize
        scattersysvec!(dchi, (L/8)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end
    return true
end

function test_convergence()
    @info "Slit cylinder"
    for n in [18, ] # 3:64 #
        _execute(n, !false)
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

