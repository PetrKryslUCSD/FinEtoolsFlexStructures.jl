"""
Twist a ribbon
"""
module ribbon_twist_examples

using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json



function _execute_dsg_model(formul, n = 2, LW_ratio = 2.0, visualize = true)
    E = 1e7;
    nu = 0.25;
    W = 1.0;
    L = LW_ratio * W
    thickness = 0.05
    analyt_sol = 26e-3 * LW_ratio / 10
    
    tolerance = thickness/1000
    fens, fes = T3block(L,W,n,1,:b);
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    
    # Report
    @info "Ribbon twist, L/W=$(LW_ratio), formulation=$(formul)"
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    associate = formul.associategeometry!
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped end
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
   
    applyebc!(dchi)
    numberdofs!(dchi);

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    Kff = stiffness(femm, massem, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[L L 0 0 -Inf Inf], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(Float64[0, 0, -1, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 3);
    nl = selectnode(fens; box = Float64[L L W W -Inf Inf], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(Float64[0, 0, +1, 0, 0, 0]);
    Ff += distribloads(lfemm, vassem, geom0, dchi, fi, 3);

    # @infiltrate
    # Solve
     Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)
    
    targetu =  dchi.values[nl, 3][1]
    @info "Target: $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_convergence()
    formul = FEMMShellT3FFModule
    LW_ratio = 2.0
    @info "Twisted Ribbon, L/W=$(LW_ratio), formulation=$(formul)"
    for n in [8, 12, 16, 24]
        _execute_dsg_model(formul, n, LW_ratio, false)
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
