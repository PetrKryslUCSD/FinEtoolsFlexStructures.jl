"""
Pinched cylinder with free ends and a concentrated force.
Mostly inextensional deformation.

Asymptotic deformation modes of benchmark
problems suitable for evaluating shell elements
Demetres Briassoulis
Comput. Methods Appl. Mech. Engrg. 194 (2005) 2385–2405
"""
module pinched_cylinder_free_ends_examples

using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function _execute_dsg_model(formul, n = 2, visualize = true)
    E = 10.5e6;
    nu = 0.3125
    R = 4.953
    L = 10.35
    P = 100.0
   
    # [8] D. Briassoulis, Testing the asymptotic behaviour of shell elements.
    # Part I. The classical benchmark tests, Int. J. Numer. Meth. Engrg. 54(3)
    # (2002) 421–452.

    # Analytical: series solution; inextensible [8,17]
    thickness = R/1.0e6 
    analyt_sol = -1.907*P*R/E/thickness^3; 
    # thickness = R/1.0e4
    # analyt_sol = -1.912*P*R/E/thickness^3;
    # Finite element solution (25 · 35 mesh) [8]
    # thickness = R/1.0e2
    # analyt_sol = -1.997*P*R/E/thickness^3;

    @info "Mesh: $n elements per side"

    # Mesh
    
    tolerance = R/n/1000
    fens, fes = T3block(90/360*2*pi,L/2,n,n,:b);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*cos(a))
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    # e0
    associategeometry! = formul.associategeometry!
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # free end
    # l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    # for i in [1,3]
    #     setebc!(dchi, l1, true, i)
    # end
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
    # plane of symmetry perpendicular to Z
    l1 = selectnode(fens; box = Float64[-Inf Inf -Inf Inf 0.0 0.0], inflate = tolerance)
    for i in [3,4,5]
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
    nl = selectnode(fens; box = Float64[0 0 L/2 L/2 -Inf Inf], inflate = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(Float64[0, 0, -P/4, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 3);

    # Solve
     Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)
    
    @show n,  dchi.values[nl, 3]/analyt_sol*100

    # Visualization
    if !visualize
        return true
    end
    scattersysvec!(dchi, (L/8)/maximum(abs.(U)).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function test_convergence(formul = FEMMShellT3FFModule)
    @info "Pinched cylinder, formulation=$(formul)"
    for n in [2, 4, 8, 16, 32] # 3:64 #
        _execute_dsg_model(formul, n, false)
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

