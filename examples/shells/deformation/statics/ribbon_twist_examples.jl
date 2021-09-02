# Simply supported square plate with uniform distributed load
module ribbon_twist_examples

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3DSGOModule
using FinEtoolsFlexStructures.FEMMShellT3DSGICModule
using FinEtoolsFlexStructures.FEMMShellT3DSGModule
using FinEtoolsFlexStructures.FEMMShellCSDSG3Module
using FinEtoolsFlexStructures.FEMMShellIsoPModule
using FinEtoolsFlexStructures.FEMMShellQ4SRIModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

using Infiltrator


function test_st3dsg(args...)
    return _execute_dsg_model(FEMMShellT3DSGModule, args...)
end

function test_st3dsgic(args...)
  return _execute_dsg_model(FEMMShellT3DSGICModule, args...)
end

function test_dsg3(args...)
  return _execute_dsg_model(FEMMShellT3DSGOModule, args...)
end

function test_csdsg3(args...)
  return _execute_dsg_model(FEMMShellCSDSG3Module, args...)
end

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

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[L L 0 0 -Inf Inf], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, 0, -1, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    nl = selectnode(fens; box = Float64[L L W W -Inf Inf], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, 0, +1, 0, 0, 0]);
    F += distribloads(lfemm, geom0, dchi, fi, 3);

    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
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

function test_t6(n = 1, LW_ratio = 10.0, visualize = true)
    E = 1e7;
    nu = 0.25;
    W = 1.0;
    L = LW_ratio * W
    thickness = 0.05
    @show analyt_sol = 26e-3 * LW_ratio / 10
    
    tolerance = thickness/1000
    fens, fes = T6block(L,W,n,n,:b);
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    formul = FEMMShellIsoPModule
    # Report
    @info "Ribbon twist, L/W=$(LW_ratio), formulation=$(formul)"
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(3), thickness), mater)
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
    numberdofs!(dchi)
    @show dchi.nfreedofs

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[L L 0 0 -Inf Inf], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, 0, -1, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    nl = selectnode(fens; box = Float64[L L W W -Inf Inf], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, 0, +1, 0, 0, 0]);
    F += distribloads(lfemm, geom0, dchi, fi, 3);

    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
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

function test_q4sri(n = 80, LW_ratio = 10.0, visualize = true)
    E = 1e7;
    nu = 0.25;
    W = 1.0;
    L = LW_ratio * W
    thickness = 0.05
    analyt_sol = 26e-3 * LW_ratio / 10
    
    tolerance = thickness/1000
    fens, fes = Q4block(L,W,n,n);
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    formul = FEMMShellQ4SRIModule
    # Report
    @info "Ribbon twist, L/W=$(LW_ratio), formulation=$(formul)"
    
    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, GaussRule(2, 2), thickness), mater)
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

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[L L 0 0 -Inf Inf], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, 0, -1, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    nl = selectnode(fens; box = Float64[L L W W -Inf Inf], tolerance = tolerance)
    loadbdry = FESetP1(reshape(nl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, 0, +1, 0, 0, 0]);
    F += distribloads(lfemm, geom0, dchi, fi, 3);

    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
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

end # module

using .ribbon_twist_examples
m = ribbon_twist_examples
# ribbon_twist_examples.test_st3dsgic()
m.test_csdsg3()
m.test_st3dsg()
# ribbon_twist_examples.test_q4sri()