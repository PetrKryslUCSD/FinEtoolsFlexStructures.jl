"""
Square plate with uniform distributed load.
Two opposite sides have one condition, the others another.

"""
module square_plate_udl_various_bcs_examples

using Arpack
using LinearAlgebra: norm
using Random: rand, Xoshiro
using CliqueTrees
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.FEMMShellQ4RSModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule: VTKWrite
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using FinEtools.MeshExportModule.CSV: savecsv
using DelimitedFiles
using PGFPlotsX

const E = 30e6
const nu = 0.3
const L = 10.0
const tL_ratio = 1/10

loading(tL_ratio) = 1.0e6 * (tL_ratio)^3
const p = loading(tL_ratio)

function _execute_q4rs_quarter_model(
    n=2,
    horizontal_support=:hard,
    vertical_support=:hard,
    mesh=:uniform,
    visualize=true
    )
    formul = FEMMShellQ4RSModule
    thickness = L * tL_ratio
    D = E / 12 / (1 - nu^2) * thickness^3
    
    if mesh == :uniform
        fens, fes = Q4block(L/2, L/2, n, n);
    elseif mesh == :biased
        xs = biasedspace(0.0, L/2, n+1, 1000/tL_ratio)
        ys = biasedspace(0.0, L/2, n+1, 1000/tL_ratio)
        fens, fes = Q4blockx(xs, ys);
    elseif mesh == :graded
        xs = gradedspace(0.0, L/2, n, 3)
        ys = gradedspace(0.0, L/2, n, 3)
        fens, fes = Q4blockx(xs, ys);
    elseif mesh == :striped
        xs = L/2 .* vcat(linearspace(0.0, tL_ratio, Int(n/2)), linearspace(tL_ratio, 1.0, Int(n/2))[2:end])
        ys = L/2 .* vcat(linearspace(0.0, tL_ratio, Int(n/2)), linearspace(tL_ratio, 1.0, Int(n/2))[2:end])
        fens, fes = Q4blockx(xs, ys);
    else
        @error "Unknown mesh"
    end
    bfes = meshboundary(fes)
    elleft = selectelem(fens, bfes; facing = true, direction = Float64[-1, 0])
    nlleft = connectednodes(subset(bfes, elleft))
    nlleft = nlleft[sortperm(fens.xyz[nlleft, 2])]
    nllefts = fens.xyz[nlleft, 2] ./ (L/2) # .+ eps(1.0)
    elright = selectelem(fens, bfes; facing = true, direction = Float64[+1, 0])
    elbott = selectelem(fens, bfes; facing = true, direction = Float64[0, -1])
    nlbott = connectednodes(subset(bfes, elbott))
    nlbott = nlbott[sortperm(fens.xyz[nlbott, 1])]
    nlbotts = fens.xyz[nlbott, 1] ./ (L/2) # .+ eps(1.0)
    eltop = selectelem(fens, bfes; facing = true, direction = Float64[0, +1])
    ncorner = selectnode(fens; box=Float64[(0.0) (0.0) (0.0) (0.0)], tolerance=eps(1.0))
    nmidfreeedge = selectnode(fens; box=Float64[(0.0) (0.0) (L/2) (L/2)], tolerance=eps(1.0))
    fens.xyz = xyz3(fens)

    ecorner = selectelem(fens, fes; withnodes = ncorner, allin=false)[1]
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    function _cartesian!(csmatout, XYZ, tangents, feid, qpid)
        csmatout[:, 1] .= (1.0, 0.0, 0.0)
        csmatout[:, 2] .= (0.0, 1.0, 0.0)
        csmatout[:, 3] .= (0.0, 0.0, 1.0)
        return csmatout
    end

    ocsys = CSys(3, 3, _cartesian!)

    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    ir = Simpson13Rule(2) # 
    # ir = GaussRule(2, 2)
    femm = formul.make(IntegDomain(fes, ir, thickness),
        mater,
    )
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    dof = [1, 2, 6] # leave only plate degrees of freedom
    for i in dof
        setebc!(dchi, 1:count(fens), true, i)
    end
    # right -- symmetry
    dof = [1, 2, 5, 6]
    for i in dof
        setebc!(dchi, connectednodes(subset(bfes, elright)), true, i)
    end
    # top -- symmetry
    dof = [1, 2, 4, 6]
    for i in dof
        setebc!(dchi, connectednodes(subset(bfes, eltop)), true, i)
    end
    # bott (horizontal)
    dof = []    # free
    if horizontal_support == :clamp
        dof = [1, 2, 3, 4, 5, 6]
    elseif horizontal_support == :hard
        dof = [1, 2, 3, 5, 6]
    elseif horizontal_support == :soft
        dof = [1, 2, 3, 6]    
    end
    for i in dof
        setebc!(dchi, connectednodes(subset(bfes, elbott)), true, i)
    end
    # left (vertical)
    dof = []    # free
    if vertical_support == :clamp
        dof = [1, 2, 3, 4, 5, 6]
    elseif vertical_support == :hard
        dof = [1, 2, 3, 4, 6]
    elseif vertical_support == :soft
        dof = [1, 2, 3, 6]    
    end
    for i in dof
        setebc!(dchi, connectednodes(subset(bfes, elleft)), true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Load
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(Float64[0, 0, -p, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 2)

    # Solve
    Kff = matrix_blocked_ff(K, nfreedofs(dchi))
    Ff = vector_blocked_f(F, nfreedofs(dchi))
    # factor = ilu0(Kff)
    # @show mKd = minimum(diag(Kff))
    # # factor = ilu(K, τ = 1e-3 * mKd) # This may work for incompressible materials
    # opM = LinearOperator(Float64, nfreedofs(dchi), nfreedofs(dchi), false, false, (y, v) -> ldiv!(y, factor, v))
    # (U, stats) = Krylov.cg(Kff, Ff; M = opM, itmax = Int(round(nfreedofs(dchi) / 10)), verbose = 0)
    # @info "$stats"
    fact = CliqueTrees.cholesky(Kff)
    U = fact \ Ff
    scattersysvec!(dchi, U)
    # solve_blocked!(dchi, K, F)
    targetu = dchi.values[nmidfreeedge, 3][1]
    @info "w(nmidfreeedge): $(round(targetu, digits=8) / (p*L^4/D/100))"

    # Visualization
    if visualize
        u = deepcopy(dchi.values[:, 1:3])
        ur = deepcopy(dchi.values[:, 4:6])
        vtkwrite("sqpl_udl-q4rs-$(horizontal_support)-$(vertical_support)-$(mesh)-tL=$(tL_ratio)-n=$(n)-uur.vtu", fens, fes; 
            vectors=[("u", u), ("ur", ur)])  
        # ocsys = CSys(3)
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys=ocsys, nodevalmethod=:averaging)
            push!(scalars, ("m$nc", fld.values))
            savecsv("sqpl_udl-left-$(horizontal_support)-$(vertical_support)-$(mesh)-$(n)-m$(nc).csv", s=nllefts, v=fld.values[nlleft])
            savecsv("sqpl_udl-bott-$(horizontal_support)-$(vertical_support)-$(mesh)-$(n)-m$(nc).csv", s=nlbotts, v=fld.values[nlbott])
        end
        vtkwrite("sqpl_udl-q4rs-$(horizontal_support)-$(vertical_support)-$(mesh)-tL=$(tL_ratio)-n=$(n)-m.vtu", fens, fes; 
            scalars=scalars, 
            vectors=[("u", u), ("ur", ur)])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys=ocsys, nodevalmethod=:averaging)
            push!(scalars, ("q$nc", fld.values))
            savecsv("sqpl_udl-left-$(horizontal_support)-$(vertical_support)-$(mesh)-$(n)-q$(nc).csv", s=nllefts, v=fld.values[nlleft])
            savecsv("sqpl_udl-bott-$(horizontal_support)-$(vertical_support)-$(mesh)-$(n)-q$(nc).csv", s=nlbotts, v=fld.values[nlbott])
           @info "q$nc Range: $(minimum(fld.values) / (p * L )) to $(maximum(fld.values) / (p * L ))"
        end
        vtkwrite("sqpl_udl-q4rs-$(horizontal_support)-$(vertical_support)-$(mesh)-tL=$(tL_ratio)-n=$(n)-q.vtu", fens, fes; 
            scalars=scalars,
            vectors=[("u", u), ("ur", ur)])

    end
    return targetu
end

function _execute_q4rs_model(
    n=2,
    horizontal_support=:hard,
    vertical_support=:hard,
    mesh=:uniform,
    visualize=true
    )
    formul = FEMMShellQ4RSModule
    thickness = L * tL_ratio
    D = E / 12 / (1 - nu^2) * thickness^3
    
    tol = if mesh == :uniform
        fens, fes = Q4block(L/2, L/2, n, n);
        L / n / 1000
    elseif mesh == :biased
        xs = biasedspace(0.0, L/2, n+1, 1000/tL_ratio)
        ys = biasedspace(0.0, L/2, n+1, 1000/tL_ratio)
        fens, fes = Q4blockx(xs, ys);
        minimum(diff(xs)) / 10l
    elseif mesh == :graded 
        xs = gradedspace(0.0, L/2, n, 3)
        ys = gradedspace(0.0, L/2, n, 3)
        fens, fes = Q4blockx(xs, ys);
        minimum(diff(xs)) / 10
    elseif mesh == :striped
        xs = L/2 .* vcat(linearspace(0.0, tL_ratio, Int(n/2)), linearspace(tL_ratio, 1.0, Int(n/2))[2:end])
        ys = L/2 .* vcat(linearspace(0.0, tL_ratio, Int(n/2)), linearspace(tL_ratio, 1.0, Int(n/2))[2:end])
        fens, fes = Q4blockx(xs, ys);
        minimum(diff(xs)) / 10
    else
        @error "Unknown mesh"
    end
    fens1, fes1 = mirrormesh(fens, fes, vec([1.0, 0.0]), vec([L/2, L/2]), renumb = (c) -> c[[2, 1, 4, 3]])
    fens, fes1, fes2 = mergemeshes(fens, fes, fens1, fes1, tol)
    fes = cat(fes1, fes2)
    fens1, fes1 = mirrormesh(fens, fes, vec([0.0, 1.0]), vec([L/2, L/2]), renumb = (c) -> c[[2, 1, 4, 3]])
    fens, fes1, fes2 = mergemeshes(fens, fes, fens1, fes1, tol)
    fes = cat(fes1, fes2)
    bfes = meshboundary(fes)
    elleft = selectelem(fens, bfes; facing = true, direction = Float64[-1, 0])
    nlleft = connectednodes(subset(bfes, elleft))
    nlleft = nlleft[sortperm(fens.xyz[nlleft, 2])]
    nllefts = fens.xyz[nlleft, 2] ./ (L) # .+ eps(1.0)
    elright = selectelem(fens, bfes; facing = true, direction = Float64[+1, 0])
    elbott = selectelem(fens, bfes; facing = true, direction = Float64[0, -1])
    nlbott = connectednodes(subset(bfes, elbott))
    nlbott = nlbott[sortperm(fens.xyz[nlbott, 1])]
    nlbotts = fens.xyz[nlbott, 1] ./ (L) # .+ eps(1.0)
    eltop = selectelem(fens, bfes; facing = true, direction = Float64[0, +1])
    ncorner = selectnode(fens; box=Float64[(0.0) (0.0) (0.0) (0.0)], tolerance=eps(1.0))
    nmidfreeedge = selectnode(fens; box=Float64[(0.0) (0.0) (L/2) (L/2)], tolerance=eps(1.0))
    fens.xyz = xyz3(fens)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    function _cartesian!(csmatout, XYZ, tangents, feid, qpid)
        csmatout[:, 1] .= (1.0, 0.0, 0.0)
        csmatout[:, 2] .= (0.0, 1.0, 0.0)
        csmatout[:, 3] .= (0.0, 0.0, 1.0)
        return csmatout
    end

    ocsys = CSys(3, 3, _cartesian!)

    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    ir = Simpson13Rule(2) # 
    # ir = GaussRule(2, 2)
    femm = formul.make(IntegDomain(fes, ir, thickness),
        mater,
    )
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    dof = [1, 2, 6] # leave only plate degrees of freedom
    for i in dof
        setebc!(dchi, 1:count(fens), true, i)
    end
    # top and bott (horizontal)
    dof = []    # free
    if horizontal_support == :clamp
        dof = [1, 2, 3, 4, 5, 6]
    elseif horizontal_support == :hard
        dof = [1, 2, 3, 5, 6]
    elseif horizontal_support == :soft
        dof = [1, 2, 3, 6]    
    end
    for i in dof
        setebc!(dchi, connectednodes(subset(bfes, elbott)), true, i)
        setebc!(dchi, connectednodes(subset(bfes, eltop)), true, i)
    end
    # left and right (vertical)
    dof = []    # free
    if vertical_support == :clamp
        dof = [1, 2, 3, 4, 5, 6]
    elseif vertical_support == :hard
        dof = [1, 2, 3, 4, 6]
    elseif vertical_support == :soft
        dof = [1, 2, 3, 6]    
    end
    for i in dof
        setebc!(dchi, connectednodes(subset(bfes, elleft)), true, i)
        setebc!(dchi, connectednodes(subset(bfes, elright)), true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Load
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(Float64[0, 0, -p, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 2)

    # Solve
    Kff = matrix_blocked_ff(K, nfreedofs(dchi))
    Ff = vector_blocked_f(F, nfreedofs(dchi))
    fact = CliqueTrees.cholesky(Kff)
    U = fact \ Ff
    scattersysvec!(dchi, U)
    # solve_blocked!(dchi, K, F)
    targetu = dchi.values[nmidfreeedge, 3][1]
    @info "w(nmidfreeedge): $(round(targetu, digits=8) / (p*L^4/D/100))"

    # Visualization
    if visualize
        u = deepcopy(dchi.values[:, 1:3])
        ur = deepcopy(dchi.values[:, 4:6])
        vtkwrite("sqpl_udl-q4rs-$(horizontal_support)-$(vertical_support)-$(mesh)-tL=$(tL_ratio)-n=$(n)-uur.vtu", fens, fes; 
            vectors=[("u", u), ("ur", ur)])  
        # ocsys = CSys(3)
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys=ocsys, nodevalmethod=:averaging)
            push!(scalars, ("m$nc", fld.values))
            savecsv("sqpl_udl-left-$(horizontal_support)-$(vertical_support)-$(mesh)-$(n)-m$(nc).csv", s=nllefts, v=fld.values[nlleft])
            savecsv("sqpl_udl-bott-$(horizontal_support)-$(vertical_support)-$(mesh)-$(n)-m$(nc).csv", s=nlbotts, v=fld.values[nlbott])
        end
        vtkwrite("sqpl_udl-q4rs-$(horizontal_support)-$(vertical_support)-$(mesh)-tL=$(tL_ratio)-n=$(n)-m.vtu", fens, fes; 
            scalars=scalars, 
            vectors=[("u", u), ("ur", ur)])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys=ocsys, nodevalmethod=:averaging)
            push!(scalars, ("q$nc", fld.values))
            savecsv("sqpl_udl-left-$(horizontal_support)-$(vertical_support)-$(mesh)-$(n)-q$(nc).csv", s=nllefts, v=fld.values[nlleft])
            savecsv("sqpl_udl-bott-$(horizontal_support)-$(vertical_support)-$(mesh)-$(n)-q$(nc).csv", s=nlbotts, v=fld.values[nlbott])
           @info "q$nc Range: $(minimum(fld.values) / (p * L )) to $(maximum(fld.values) / (p * L ))"
        end
        vtkwrite("sqpl_udl-q4rs-$(horizontal_support)-$(vertical_support)-$(mesh)-tL=$(tL_ratio)-n=$(n)-q.vtu", fens, fes; 
            scalars=scalars,
            vectors=[("u", u), ("ur", ur)])

    end
    return targetu
end

const VISUALIZE = true
const NS = [32]
const STAB_ALPHA = 0.1
const SUPPORTS = [
    (:hard, :hard), 
    (:soft, :soft), 
    (:hard, :free), 
    (:soft, :free), 
    (:clamp, :hard),
    (:clamp, :soft),
    (:clamp, :free),
    ]

function test_q4rs()
    @info "thickness/length = $tL_ratio"
    for support in SUPPORTS
        @info "Support $support --------------------------------------------------"
        for mesh in [:graded] # :biased
            @info "Mesh distortion: $mesh"
            for n in NS
                @info "n = $n"
                _execute_q4rs_model(n, support[1], support[2], mesh, VISUALIZE)
            end
        end
    end
    return true
end

const COLORS = ["black", "red", "green", "blue", "cyan", "magenta", "yellow", "gray"]

function plot_curve(objects, support, A, color)
    @pgf o = PGFPlotsX.Plot(
        {
        color = color,
        line_width  = 1.0
        },
        Coordinates([v for v in  zip(A[:,1], A[:,2] ./ (p * L ))])
        )
    push!(objects, o)
    push!(objects, LegendEntry("$(support[1])-$(support[2])"))
end

function plot()
    N = NS[1]
    for edge in ["bott", "left"]
        for nc in [1, 2]
            objects = []
            for (support, color) in zip(SUPPORTS, COLORS)
                f = "sqpl_udl-$(edge)-$(support[1])-$(support[2])-graded-$(N)-q$(nc).csv"
                A = readdlm(f, ',', Float64; skipstart=1)
                plot_curve(objects, support, A, color)
            end
            @pgf ax = Axis(
                {
                    title = "$(edge)",
                    xlabel = "Distance from corner",
                    ylabel = "q$(nc)",
                    xmode = "linear",
                    ymode = "linear",
                    yminorgrids = "true",
                    grid = "both",
                    legend_style = {
                        at = Coordinate(1.005, 0.5),
                        anchor = "west",
                    },
                },
                objects...
            )
            display(ax)
            pgfsave("sqpl_udl-$(edge)-q$(nc).pdf", ax)
        end
    end
    return true
end

function allrun()
    println("#####################################################")
    println("# test_q4rs ")
    test_q4rs()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing

