# Simply supported square plate with uniform distributed load
# Note: In order to get fast convergence, the hard SS needs 
# to be adopted.
module simply_supp_square_plate_udl_examples

using Arpack
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.FEMMShellQ4RSModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

const E = 30e6;
const nu = 0.3;
const L = 10.0;

loading(t) = 1.0e10*(t/L)^3

# From A Triangular Plate Bending Element Based on an Energy-Orthogonal Free, 
    # by Felippa and Bergan, 1986, Table 3. 
w_analyt_sol(p, L, D) = -4.06235e-3 * p * L^4 / D

function _cond(A)
    emax = eigs(A, nev = 1, which = :LM, tol = 1e-3)[1]
    emin = eigs(A, nev = 1, which = :SM, tol = 1e-3)[1]
    return (emax[1]) / (emin[1])
end

function   cartesian!(csmatout, XYZ, tangents, feid, qpid)
    csmatout[:, 1] .= (1.0, 0.0, 0.0)
    csmatout[:, 2] .= (0.0, 1.0, 0.0)
    csmatout[:, 3] .= (0.0, 0.0, 1.0)
    return csmatout
end

function _execute_t3ff_quarter_model(n = 2, tL_ratio = 0.01, visualize = true)
    formul = FEMMShellT3FFModule
    thickness = L * tL_ratio;
    D = E*thickness^3/12/(1-nu^2)
    p = loading(thickness)
    analyt_sol = w_analyt_sol(p, L, D)

    tolerance = L/n/1000
    fens, fes = T3block(L/2,L/2,n,n);
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    ocsys = CSys(3, 3, cartesian!)

    @info "Mesh: $n elements per side"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), 
        mater,
        (t, h) -> t^2 / (t^2 + 0.0 * h^2),
    )
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 0 L/2 -Inf Inf], inflate = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[0 L/2 0 0 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    # simple support
    l1 = selectnode(fens; box = Float64[L/2 L/2 0 L/2 -Inf Inf], inflate = tolerance)
    for i in [3,4,6]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 L/2 L/2 L/2 -Inf Inf], inflate = tolerance)
    for i in [3,5,6]
        setebc!(dchi, l1, true, i)
    end
    # in-plane, rotations
    l1 = selectnode(fens; box = Float64[0 L/2 0 L/2 -Inf Inf], inflate = tolerance)
    for i in [1, 2, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], tolerance = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64[0, 0, -p, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    targetu =  dchi.values[nl, 3][1]
    @info "w(center): $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

    # Visualization
    if visualize
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("q$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("eq$nc", fld.values))
        end
        vtkwrite("simply_supp_square_plate_udl-quarter-$(n)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

        scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end
    return true
end

function _execute_Q4RS_quarter_model(n = 2, tL_ratio = 0.01, visualize = true)
    formul = FEMMShellQ4RSModule
    thickness = L * tL_ratio
    D = E*thickness^3/12/(1-nu^2)
    p = loading(thickness)
    # analytical solution for the vertical deflection under the load
    # From A Triangular Plate Bending Element Based on an Energy-Orthogonal Free, 
    # by Felippa and Bergan, 1986, Table 3. 
    analyt_sol=-4.06235e-3*p*L^4/D;

    tolerance = L/n/1000
    fens, fes = T3block(L/2,L/2,n,n);
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    ocsys = CSys(3, 3, cartesian!)

    @info "Mesh: $n elements per side"

    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, GaussRule(2, 2), thickness), 
        mater, 
        (t, h) -> t^2 / (t^2 + 0.0 * h^2)
        )
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 0 L/2 -Inf Inf], inflate = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[0 L/2 0 0 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    # simple support
    l1 = selectnode(fens; box = Float64[L/2 L/2 0 L/2 -Inf Inf], inflate = tolerance)
    for i in [3,4,6]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 L/2 L/2 L/2 -Inf Inf], inflate = tolerance)
    for i in [3,5,6]
        setebc!(dchi, l1, true, i)
    end
    # in-plane, rotations
    l1 = selectnode(fens; box = Float64[0 L/2 0 L/2 -Inf Inf], inflate = tolerance)
    for i in [1, 2, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], tolerance = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(Float64[0, 0, -p, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    targetu =  dchi.values[nl, 3][1]
    @info "w(center): $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

    # Visualization
    if visualize

        # Generate a graphical display of resultants
        # scalars = []
        # for nc in 1:3
        #     fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
        #     push!(scalars, ("m$nc", fld.values))
        #     fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
        #     push!(scalars, ("em$nc", fld.values))
        # end
        # vtkwrite("scordelis_lo_examples-$(n)-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        # scalars = []
        # for nc in 1:3
        #     fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
        #     push!(scalars, ("n$nc", fld.values))
        #     fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
        #     push!(scalars, ("en$nc", fld.values))
        # end
        # vtkwrite("scordelis_lo_examples-$(n)-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("q$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("eq$nc", fld.values))
        end
        vtkwrite("simply_supp_square_plate_udl-quarter-$(n)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

        # scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
        # update_rotation_field!(Rfield0, dchi)
        # plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        # #plot_nodes(fens),
        #     plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
        #     dims = 1)
        # pl = render(plots)
    end
    return true
end

function test_convergence_quarter()
    formul = FEMMShellT3FFModule
    tL_ratio = 0.1
    @info "Simply supported square plate with uniform load, T3FF"
    @info "thickness/length = $tL_ratio"
    # for n in [32,  ]
    for n in [2, 4, 8, 16, 32, 64]
        _execute_t3ff_quarter_model(n, tL_ratio, true)
    end
    formul = FEMMShellQ4RSModule
    tL_ratio = 0.1
    @info "Simply supported square plate with uniform load, Q4RS"
    @info "thickness/length = $tL_ratio"
    # for n in [32,  ]
    for n in [2, 4, 8, 16, 32, 64]
        _execute_Q4RS_quarter_model(n, tL_ratio, true)
    end
    return true
end

function _execute_t3ff_full_model(
    n=2,
    tL_ratio=0.01,
    simple_support=:hard,
    stab_alpha=0.1,
    tweak_mesh=false,
    visualize=true
    )
    formul = FEMMShellT3FFModule
    thickness = L * tL_ratio;
    D = E/12/(1-nu^2)*thickness^3
    p = loading(thickness)
    analyt_sol = w_analyt_sol(p, L, D)

    tolerance = L/n/1000
    fens, fes = T3block(L,L,n,n);
    fens.xyz = xyz3(fens)
    fens.xyz[:, 1] .-= L/2
    fens.xyz[:, 2] .-= L/2
    nm = selectnode(fens; box = Float64[-L/2+L/n -L/2+L/n -L/2+L/n -L/2+L/n -Inf Inf], inflate = tolerance)
    if tweak_mesh
        shift = L/n/4
        fens.xyz[nm, 1] .+= shift
        fens.xyz[nm, 2] .+= shift
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    ocsys = CSys(3, 3, cartesian!)

    @info "Mesh: $n elements per side"

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), 
        mater,
        (t, h) -> t^2 / (t^2 + stab_alpha * h^2),
        )
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # edges parallel to Y, soft support: rotations in plane are free; hard simple support: rotation orthogonal to edge = 0
    dof = [1,2,3,6]
    simple_support == :hard && (dof = [1,2,3,4,6])
    l1 = selectnode(fens; box = Float64[-L/2 -L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[+L/2 +L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    # edges parallel to X, soft support: rotations in plane are free; hard simple support: rotation orthogonal to edge = 0
    dof = [1,2,3,6]
    simple_support == :hard && (dof = [1,2,3,5,6])
    l1 = selectnode(fens; box = Float64[-Inf Inf -L/2 -L/2 -Inf Inf], inflate = tolerance)
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[-Inf Inf +L/2 +L/2 -Inf Inf], inflate = tolerance)
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nl = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], tolerance = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64[0, 0, -p, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    targetu =  dchi.values[nl, 3][1]
    @info "w(center): $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"
    epsrel = abs(targetu - analyt_sol) / abs(analyt_sol)
    @info "Digits of accuracy: $(-log10(epsrel))"
    # Visualization
    if visualize
        vtkwrite("ss_sqpl_udl-t3ff-full-$(simple_support)-r=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        # ocsys = CSys(3)
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys=ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
            push!(scalars, ("m$nc", fld.values))
        end
        vtkwrite("ss_sqpl_udl-t3ff-full-$(simple_support)-r=$(tL_ratio)-s=$(stab_alpha)-m.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys=ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("q$nc", fld.values))
        end
        vtkwrite("ss_sqpl_udl-t3ff-full-$(simple_support)-r=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-q.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    end
    return true
end

function _execute_q4rs_full_model(
    n=2,
    tL_ratio=0.01,
    simple_support=:hard,
    stab_alpha=0.1,
    tweak_mesh=false,
    visualize=true
    )
    formul = FEMMShellQ4RSModule
    thickness = L * tL_ratio;
    D = E/12/(1-nu^2)*thickness^3
    p = loading(thickness)
    analyt_sol = w_analyt_sol(p, L, D)
   
    tolerance = L/n/1000
    fens, fes = Q4block(L,L,n,n);
    fens.xyz = xyz3(fens)
    fens.xyz[:, 1] .-= L/2
    fens.xyz[:, 2] .-= L/2
    nm = selectnode(fens; box = Float64[-L/2+L/n -L/2+L/n -L/2+L/n -L/2+L/n -Inf Inf], inflate = tolerance)
    if tweak_mesh
        shift = L/n/4
        fens.xyz[nm, 1] .+= shift
        fens.xyz[nm, 2] .+= shift
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    ocsys = CSys(3, 3, cartesian!)

    @info "Mesh: $n elements per side"

    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, GaussRule(2, 2), thickness), 
        mater,
        (t, h) -> t^2 / (t^2 + stab_alpha * h^2),
        )
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # edges parallel to Y, soft support: rotations in plane are free; hard simple support: rotation orthogonal to edge = 0
    dof = [1,2,3,6]
    simple_support == :hard && (dof = [1,2,3,4,6])
    l1 = selectnode(fens; box = Float64[-L/2 -L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[+L/2 +L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    # edges parallel to X, soft support: rotations in plane are free; hard simple support: rotation orthogonal to edge = 0
    dof = [1,2,3,6]
    simple_support == :hard && (dof = [1,2,3,5,6])
    l1 = selectnode(fens; box = Float64[-Inf Inf -L/2 -L/2 -Inf Inf], inflate = tolerance)
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[-Inf Inf +L/2 +L/2 -Inf Inf], inflate = tolerance)
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    setebc!(dchi, 1:count(fens), true, 6)
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    # @info "stab_alpha = $stab_alpha, Condition number: $(_cond(matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:ff]))"

    # Load
    nl = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], tolerance = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(Float64[0, 0, -p, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);

    # @infiltrate
    # Solve
    solve_blocked!(dchi, K, F)
    targetu =  dchi.values[nl, 3][1]
    @info "w(center): $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"
    epsrel = abs(targetu - analyt_sol) / abs(analyt_sol)
    @info "Digits of accuracy: $(-log10(epsrel))"
    # Visualization
    if visualize
        vtkwrite("ss_sqpl_udl-q4rs-full-$(simple_support)-r=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
        # ocsys = CSys(3)
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys=ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
            push!(scalars, ("m$nc", fld.values))
        end
        vtkwrite("ss_sqpl_udl-q4rs-full-$(simple_support)-r=$(tL_ratio)-s=$(stab_alpha)-m.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys=ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("q$nc", fld.values))
        end
        vtkwrite("ss_sqpl_udl-q4rs-full-$(simple_support)-r=$(tL_ratio)-s=$(stab_alpha)-n=$(n)-q.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    end
    return true
end

function test_convergence_full()
    formul = FEMMShellT3FFModule
    stab_alpha = 0.1
    visualize = false
    support = :soft
    tweak_mesh = false
    tL_ratio = 0.00001
    tL_ratio = 0.001
    @info "Simply supported square plate with uniform load, T3FF"
    @info "thickness/length = $tL_ratio formulation=$(formul)"
    # for n in [32,  ]
    for n in [2, 4, 8, 16, 32, 64, 128, 256, 512]
        @info "---------------------------------------------------------------------"
        _execute_t3ff_full_model(n, tL_ratio, support, stab_alpha, visualize)
        formul = FEMMShellT3FFModule
    end
    formul = FEMMShellQ4RSModule
    tL_ratio = 0.00001
    tL_ratio = 0.001
    @info "Simply supported square plate with uniform load, Q4RS"
    @info "thickness/length = $tL_ratio formulation=$(formul)"
    # for n in [32,  ]
    for n in [2, 4, 8, 16, 32, 64, 128, 256, 512]
        @info "---------------------------------------------------------------------"
        _execute_q4rs_full_model(n, tL_ratio, support, stab_alpha, tweak_mesh, visualize)
    end
    return true
end

function test_convergence_full_thickness(tL_ratios = [1.0e-1, 1.0e-2, 1.0e-4, 1.0e-8])
    stab_alpha = 0.001
    visualize = !false
    support = :soft
    tweak_mesh = false
    ns = [2, 4, 8, 16, 32, 64]
    formul = FEMMShellT3FFModule
    @info "Simply supported square plate with uniform load, T3FF --------------------"
    for tL_ratio in tL_ratios
        @info "thickness/length = $tL_ratio"
        for n in ns
            _execute_t3ff_full_model(n, tL_ratio, support, stab_alpha, visualize)
            formul = FEMMShellT3FFModule
        end
    end
    formul = FEMMShellQ4RSModule
    @info "Simply supported square plate with uniform load, Q4RS --------------------"
    for tL_ratio in tL_ratios
        @info "thickness/length = $tL_ratio"
        for n in ns
            _execute_q4rs_full_model(n, tL_ratio, support, stab_alpha, tweak_mesh, visualize)
        end
    end
    return true
end

function test_stabilization_q4rs(tL_ratios = [1.0e-2, 1.0e-4, 1.0e-8])
    stab_alpha = 0.001
    visualize = true
    support = :hard
    tweak_mesh = true
    n = 16
    formul = FEMMShellQ4RSModule
    @info "Simply supported square plate with uniform load, Q4RS --------------------"
    for tL_ratio in tL_ratios
        @info ">>>>>>>>>>>>> thickness/length = $tL_ratio"
        for stab_alpha in [0.0, 0.000001, 0.001, 0.1]
            @info "  stab_alpha = $stab_alpha"
            _execute_q4rs_full_model(n, tL_ratio, support, stab_alpha, tweak_mesh, visualize)
        end
    end
    return true
end

function test_stabilization_t3ff(tL_ratios = [1.0e-2, 1.0e-4, 1.0e-8])
    stab_alpha = 0.001
    visualize = true
    support = :hard
    tweak_mesh = true
    n = 16
    formul = FEMMShellT3FFModule
    @info "Simply supported square plate with uniform load, T3FF --------------------"
    for tL_ratio in tL_ratios
        @info ">>>>>>>>>>>>> thickness/length = $tL_ratio"
        for stab_alpha in [0.0, 0.000001, 0.001, 0.1]
            @info "  stab_alpha = $stab_alpha"
            _execute_t3ff_full_model(n, tL_ratio, support, stab_alpha, tweak_mesh, visualize)
        end
    end
    return true
end

function allrun()
    println("#####################################################")
    println("# test_stabilization_t3ff ")
    test_stabilization_t3ff()
    println("#####################################################")
    println("# test_stabilization_q4rs ")
    test_stabilization_q4rs()
    println("#####################################################")
    println("# test_convergence_full ")
    # test_convergence_full()
    println("#####################################################")
    println("# test_convergence_full_thickness ")
    # test_convergence_full_thickness()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
