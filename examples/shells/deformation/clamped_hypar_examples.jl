"""
From: CMES, vol.49, no.2, pp.81-110, 2009

The problem considered in this section is that of a hyperbolic paraboloid shell,
clamped along one side and free on three edges and loaded by self-weight
(Figure 16). This is a pure bending dominated problem and known to be a very
hard test for locking behaviour as suggested in References (Chapelle and Bathe,
1998; Bathe, Iosilevich, and Chapelle, 2000). 

The shell geometry is described by the
equation: z = x^2 −y^2 ; (x,y) ∈ − L/2 ; L/2
"""
module clamped_hypar_examples

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellDSG3Module
# using FinEtoolsFlexStructures.FEMMShellT3Module
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function test_dsg3()
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    # Parameters:
    E = 2.0e11;
    nu = 0.3;
    L = 1.0
    # Bathe, Iosilevich, and Chapelle (2000) with a refined mesh of
    # high-order element MITC16

    thickness = L/100; 
    analyt_sol=-9.3355e-5;
    g = 80*0.1^0

    thickness = L/1000; 
    analyt_sol=-6.3941e-5;
    g = 80*0.1^3#
    
    # Mesh
    n = 32

    formul = FEMMShellDSG3Module
    # Report
    @info "Clamped hypar, t/L=$(thickness/L), formulation=$(formul)"
    @info "Mesh: $n elements per side"

    tolerance = L/n/1000
    fens, fes = T3block(L,L,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        x=fens.xyz[i, 1]-L/2; y=fens.xyz[i, 2]-L/2;
        fens.xyz[i, :] .= (x, y, x^2 - y^2)
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.FEMMShellDSG3(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped edge
    l1 = selectnode(fens; box = Float64[-L/2 -L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[L/2 L/2 0 0 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    # computeforce!(forceout, XYZ, tangents, fe_label) = let
    #     n = cross(tangents[:, 1], tangents[:, 2])
    #     n = n / norm(n)
    #     forceout .= 0.0
    #     forceout[3] = -g * n[3]
    #     return forceout
    # end
    # fi = ForceIntensity(FFlt, 6, computeforce!)
    fi = ForceIntensity(FFlt[0, 0, -g, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    targetu = dchi.values[nl, 3][1]
    @info "Target: $targetu,  $(round(targetu/analyt_sol, digits = 4)*100)%"

        # Visualization
    scattersysvec!(dchi, (L/4)/abs(targetu).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end


function test_q4sri()
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    # Parameters:
    E = 2.0e11;
    nu = 0.3;
    L = 1.0
    # Bathe, Iosilevich, and Chapelle (2000) with a refined mesh of
    # high-order element MITC16

    thickness = L/100; 
    analyt_sol=-9.3355e-5;
    g = 80*0.1^0

    thickness = L/1000; 
    analyt_sol=-6.3941e-5;
    g = 80*0.1^3
    
    # Mesh
    n = 32

    formul = FEMMShellQ4SRIModule
    # Report
    @info "Clamped hypar, t/L=$(thickness/L), formulation=$(formul)"
    @info "Mesh: $n elements per side"

    tolerance = L/n/1000
    fens, fes = Q4block(L,L,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        x=fens.xyz[i, 1]-L/2; y=fens.xyz[i, 2]-L/2;
        fens.xyz[i, :] .= (x, y, x^2 - y^2)
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.FEMMShellQ4SRI(IntegDomain(fes, GaussRule(2, 2), thickness), mater)
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped edge
    l1 = selectnode(fens; box = Float64[-L/2 -L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[L/2 L/2 0 0 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    # computeforce!(forceout, XYZ, tangents, fe_label) = let
    #     n = cross(tangents[:, 1], tangents[:, 2])
    #     n = n / norm(n)
    #     forceout .= 0.0
    #     forceout[3] = -g * n[3]
    #     return forceout
    # end
    # fi = ForceIntensity(FFlt, 6, computeforce!)
    fi = ForceIntensity(FFlt[0, 0, -g, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    targetu = dchi.values[nl, 3][1]
    @info "Target: $targetu,  $(round(targetu/analyt_sol, digits = 4)*100)%"

    # Visualization
    scattersysvec!(dchi, (L/4)/abs(targetu).*U)
    update_rotation_field!(Rfield0, dchi)
    plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
        plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
        plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    dims = 1)
    pl = render(plots)
    return true
end

function allrun()
    println("#####################################################")
    println("# test_dsg3 ")
    test_dsg3()
    println("#####################################################")
    println("# test_q4sri ")
    test_q4sri()
    return true
end # function allrun

end # module

using .clamped_hypar_examples
clamped_hypar_examples.allrun()