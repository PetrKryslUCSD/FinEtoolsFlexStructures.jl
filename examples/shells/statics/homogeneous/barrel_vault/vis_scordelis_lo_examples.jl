"""
The barrel vault (Scordelis-Lo) roof is one of the benchmarks for linear elastic
analysis of shells. 

The candidate element's usefulness in irregular geometries (and most practical
cases involve a high degree of geometric irregularity) is tested. As would be
expected,the irregular mesh results are not as good as those provided by a
regular mesh with the same number of variables. 

Problem description

The physical basis of the problem is a deeply arched roof supported only
by diaphragms at its curved edges (an aircraft hangar), deforming under its own
weight. It is interesting to observe that the geometry is such that the
centerpoint of the roof moves upward under the self-weight(downwardly directed)
load. Perhaps this is one reason why the problem is not straightforward
numerically.

Analytical solution for the vertical deflection and the midpoint of the
free edge is often cited as 0.3024. However, this number is suspect. 
"""
module vis_scordelis_lo_examples

using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtools.AlgoBaseModule: richextrapol
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
using RichardsonExtrapolationUQ
using DelimitedFiles
using PGFPlotsX

# Analytical solution?
const ANALYT_SOL = Dict(:soft => -0.3020326, :hard => -0.3014486)

# Parameters:
const E = 4.32e8;
const nu = 0.0;
const thickness = 0.25; # geometrical dimensions are in feet
const R = 25.0;
const L = 50.0;

cylindrical!(csmatout, XYZ, tangents, feid, qpid) = begin
    r = vec(deepcopy(XYZ));
    r[2] = 0.0;
    r[3] += R # see def of Z below
    csmatout[:, 3] .= vec(r)/norm(vec(r))
    csmatout[:, 2] .= (0.0, 1.0, 0.0) #  this is along the axis
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

function _execute_t3ff(formul, n=8, visualize=true)
    tolerance = R/n/10
    # fens, fes = T3blockrand(40/360*2*pi,L/2,n,n);
    xs = reverse(gradedspace(0.0, 40/360*2*pi, n, 4))
    ys = reverse(gradedspace(0.0, L/2, n, 4))
    fens, fes = T3blockx(xs, ys);
    # fens, fes = T3block(40/360*2*pi,L/2,n,n,:b);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1];
        y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*(cos(a)-1))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # rigid diaphragm
    l1 = selectnode(fens; box=Float64[-Inf Inf 0 0 -Inf Inf], inflate=tolerance)
    for i in [1, 3, 5]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box=Float64[-Inf Inf L/2 L/2 -Inf Inf], inflate=tolerance)
    for i in [2, 4, 6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box=Float64[0 0 -Inf Inf -Inf Inf], inflate=tolerance)
    for i in [1, 5, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    # vtkwrite("debug-normals-t3ff-$(n).vtu", fens, fes; vectors = [("normals", deepcopy(femm._normals[:, 1:3]))])
    Kff = stiffness(femm, massem, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box=Float64[sin(40/360*2*pi)*25 sin(40/360*2*pi)*25 L/2 L/2 -Inf Inf], inflate=tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64[0, 0, -90, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 3);

    # Solve
    Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)

    result = dchi.values[nl, 3][1]
    @info "Solution: $(result), $(round(result/ANALYT_SOL*100, digits = 4))%"

    # Visualization
    if visualize
        # Generate a graphical display of resultants
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys=ocsys)
            push!(scalars, ("m$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys=ocsys)
            push!(scalars, ("em$nc", fld.values))
        end
        vtkwrite("scolo_t3ff-$(n)-m.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys=ocsys)
            push!(scalars, ("n$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys=ocsys)
            push!(scalars, ("en$nc", fld.values))
        end
        vtkwrite("scolo_t3ff-$(n)-n.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys=ocsys)
            push!(scalars, ("q$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys=ocsys)
            push!(scalars, ("eq$nc", fld.values))
        end
        vtkwrite("scolo_t3ff-$(n)-q.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])

        # vtkwrite("scolo_t3ff-$(n)-uur.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])

        # scattersysvec!(dchi, (L/8)/maximum(abs.(U)).*U)
        # update_rotation_field!(Rfield0, dchi)
        # plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #     #plot_nodes(fens),
        #     plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
        #     plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
        #     dims = 1)
        # pl = render(plots)
    end

    result
end


function _closest(loc, x, tol)
    mind = Inf
    c = 0 # no closest point (i.e. no point within tolerance) 
    for i in axes(x, 1)
        d = norm(vec(x[i, :]) - vec(loc))
        if d < mind
            mind = d
            c = i
        end
    end
    return c
end

struct Idat
    nc::Int
    v::Vector{Float64}
    tol::Float64
end

function _inspector(idat, j, conn, x, out, loc)
    c = _closest(loc, x, idat.tol)
    if c > 0
        idat.v[c] = out[idat.nc]
    end
    return idat
end

function _elementwise_arrays(femm, geom0, dchi, quantity, nc, tol)
    fes = finite_elements(femm)
    connectivity = zeros(Int, count(fes), 4)
    points = zeros(4 * count(fes), 3)
    values = zeros(4 * count(fes))
    idat = Idat(nc, zeros(4), tol)
    ecoords = zeros(4, 3)
    for i in eachindex(fes)
        gathervalues_asmat!(geom0, ecoords, fes.conn[i])
        inspectintegpoints(femm, geom0, dchi, [i], _inspector, idat, quantity)
        r = (i - 1) * 4 .+ (1:4)
        connectivity[i, :] .= r
        points[r, :] = ecoords
        values[r] = idat.v
        idat.v .= 0
    end
    return connectivity, points, values
end

function _execute_q4rs(mesh=:uniform, n=8, support = :soft)
    formul = FEMMShellQ4RSModule
    @info "Mesh: $mesh; $n elements per side"
    bias = 100
    if mesh == :uniform
        fens, fes = Q4block(40/360*2*pi, L/2, n, 2*n);
        tolerance = L/2/n/100
    elseif mesh == :biased
        xs = 40/360*2*pi .- reverse(biasedspace(0.0, 40/360*2*pi, n+1, bias))
        ys = biasedspace(0.0, L/2, 2*n+1, bias)
        fens, fes = Q4blockx(xs, ys);
        tolerance = (minimum(abs.(diff(xs))) + minimum(abs.(diff(ys)))) / 100
    else 
        @error "Unknown mesh"
    end
    bfes = meshboundary(fes)
    ela0 = selectelem(fens, bfes; facing=true, direction=Float64[-1, 0])
    ela1 = selectelem(fens, bfes; facing=true, direction=Float64[+1, 0])
    ell0 = selectelem(fens, bfes; facing=true, direction=Float64[0, -1])
    ell1 = selectelem(fens, bfes; facing=true, direction=Float64[0, +1])
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1];
        y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*(cos(a)-1))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, GaussRule(2, 2), thickness), ocsys, mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # rigid diaphragm: SOFT simple support
    l1 = connectednodes(subset(bfes, ell0))
    dof = support == :soft ? [1, 3] : [1, 3, 5]
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = connectednodes(subset(bfes, ell1))
    for i in [2, 4, 6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = connectednodes(subset(bfes, ela0))
    for i in [1, 5, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    # vtkwrite("debug-normals-q4rs-$(n).vtu", fens, fes; vectors = [("normals", deepcopy(femm._normals[:, 1:3]))])
    Kff = stiffness(femm, massem, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box=Float64[sin(40/360*2*pi)*25 sin(40/360*2*pi)*25 L/2 L/2 -Inf Inf], inflate=tolerance)
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(Float64[0, 0, -90, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 2);

    # Solve
    Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)

    result = dchi.values[nl, 3][1]
    @info "Solution: $(result), $(round(result/ANALYT_SOL[support]*100, digits = 4))%"

    midsection_nodes = connectednodes(subset(bfes, ell1))
    midsection_nodes_x = sortperm(fens.xyz[midsection_nodes, 1])
    midsection_nodes_ordered = midsection_nodes[midsection_nodes_x]
    midsection_angles = asin.(fens.xyz[midsection_nodes_ordered, 1] ./ R) * 180 / pi / 40

    diaphragm_nodes = connectednodes(subset(bfes, ell0))
    diaphragm_nodes_x = sortperm(fens.xyz[diaphragm_nodes, 1])
    diaphragm_nodes_ordered = diaphragm_nodes[diaphragm_nodes_x]
    diaphragm_angles = asin.(fens.xyz[diaphragm_nodes_ordered, 1] ./ R) * 180 / pi / 40

    peak_nodes = connectednodes(subset(bfes, ela0))
    peak_nodes_y = sortperm(fens.xyz[peak_nodes, 2])
    peak_nodes_ordered = peak_nodes[peak_nodes_y]
    peak_dists = fens.xyz[peak_nodes_ordered, 2] ./ L

    free_nodes = connectednodes(subset(bfes, ela1))
    free_nodes_y = sortperm(fens.xyz[free_nodes, 2])
    free_nodes_ordered = free_nodes[free_nodes_y]
    free_dists = fens.xyz[free_nodes_ordered, 2] ./ L

    # Visualization
    basef = "scolo_q4rs-$(support)-$(mesh)-$(n)"
    # Generate a graphical display of resultants
    # Here we generate an ad hoc machine: based on a single integration point, it improves the post processing of the stresses, especially the shear membrane.
    pfemm = formul.make(IntegDomain(fes, GaussRule(2, 1), thickness), ocsys, mater)
    associategeometry!(pfemm, geom0)
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(pfemm, geom0, dchi, :moment, nc, outputcsys=ocsys)
        push!(scalars, ("m$nc", fld.values))
        @info "m$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
        savecsv("$(basef)-midsection-m$(nc).csv", a=midsection_angles, v=fld.values[midsection_nodes_ordered])
        savecsv("$(basef)-diaphragm-m$(nc).csv", a=diaphragm_angles, v=fld.values[diaphragm_nodes_ordered])
        savecsv("$(basef)-peak-m$(nc).csv", a=peak_dists, v=fld.values[peak_nodes_ordered])
        savecsv("$(basef)-free-m$(nc).csv", a=free_dists, v=fld.values[free_nodes_ordered])
        fld = elemfieldfromintegpoints(pfemm, geom0, dchi, :moment, nc, outputcsys=ocsys)
        push!(scalars, ("em$nc", fld.values))
        @info "em$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
    end
    vtkwrite("$(basef)-m.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(pfemm, geom0, dchi, :membrane, nc, outputcsys=ocsys)
        push!(scalars, ("n$nc", fld.values))
        @info "n$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
        savecsv("$(basef)-midsection-n$(nc).csv", a=midsection_angles, v=fld.values[midsection_nodes_ordered])
        savecsv("$(basef)-diaphragm-n$(nc).csv", a=diaphragm_angles, v=fld.values[diaphragm_nodes_ordered])
        savecsv("$(basef)-peak-n$(nc).csv", a=peak_dists, v=fld.values[peak_nodes_ordered])
        savecsv("$(basef)-free-n$(nc).csv", a=free_dists, v=fld.values[free_nodes_ordered])
        fld = elemfieldfromintegpoints(pfemm, geom0, dchi, :membrane, nc, outputcsys=ocsys)
        push!(scalars, ("en$nc", fld.values))
        @info "en$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
    end
    vtkwrite("$(basef)-n.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    # for nc in 1:3
    #     scalars = []
    #     connectivity, points, values = _elementwise_arrays(pfemm, geom0, dchi, :membrane, nc, tolerance)
    #     push!(scalars, ("dn$nc", deepcopy(values)))
    #     vtkwrite("$(basef)-dn$nc.vtu", connectivity, points, VTKWrite.Q4; scalars=scalars)
    # end
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(pfemm, geom0, dchi, :shear, nc, outputcsys=ocsys)
        push!(scalars, ("q$nc", fld.values))
        @info "q$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
        savecsv("$(basef)-midsection-q$(nc).csv", a=midsection_angles, v=fld.values[midsection_nodes_ordered])
        savecsv("$(basef)-diaphragm-q$(nc).csv", a=diaphragm_angles, v=fld.values[diaphragm_nodes_ordered])
        savecsv("$(basef)-peak-q$(nc).csv", a=peak_dists, v=fld.values[peak_nodes_ordered])
        savecsv("$(basef)-free-q$(nc).csv", a=free_dists, v=fld.values[free_nodes_ordered])
        fld = elemfieldfromintegpoints(pfemm, geom0, dchi, :shear, nc, outputcsys=ocsys)
        push!(scalars, ("eq$nc", fld.values))
        @info "eq$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
    end
    vtkwrite("$(basef)-q.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    
    result
end


const COLORS = ["black", "red", "green", "blue", "cyan", "magenta", "yellow", "gray"]
const MARKERS = [
    "square",
    "triangle",
    "diamond",
    "triangle*",
    "x",
    "square*",
    "diamond*",
]
const MARK_REPEAT = [13, 15, 17, 19, 23, 27, 29] .+ 73

function plot(basef = "", 
              res = "q")
    ncs = res == "q" ? (1:2) : (1:3)
    for edge in ["midsection", "diaphragm", "peak", "free"]
        objects = []
        for nc in ncs
            f = "$(basef)-$(edge)-$(res)$(nc).csv"
            A = readdlm(f, ',', Float64; skipstart=1)
            @pgf o = PGFPlotsX.Plot(
            {
            color = COLORS[nc],
            mark=MARKERS[nc], mark_size=2.5, mark_repeat=MARK_REPEAT[nc],
            line_width  = 1.0
            },
            Coordinates([v for v in  zip(A[:,1], A[:,2])])
            )
            push!(objects, o)
            # push!(objects, LegendEntry("$(res)$(nc)"))
        end
        @pgf ax = Axis(
            {
                title = "$(edge)",
                xlabel = "Normalized distance",
                ylabel = "Resultant",
                # xmin = -0.01,
                # xmax = 1.01,
                xmode = "linear",
                ymode = "linear",
                yminorgrids = "true",
                grid = "both",
                # legend_style = {
                #     at = Coordinate(1.005, 0.5),
                #     anchor = "west",
                # },
            },
            objects...
        )
        display(ax)
        pgfsave("$(basef)-$(edge)-$(res).pdf", ax)
    end
    return true
end

const CORNERS = [
    (:A, "midsection", :b),
    (:B, "midsection", :e),
    (:C, "diaphragm", :e),
    (:D, "diaphragm", :b),
    ]

function extrapolate(basef = "", ns = [16, 32, 64, 128], res = "q")
    ncs = res == "q" ? (1:2) : (1:3)
    results = [] 
    for corner in CORNERS
        edge = corner[2]
        for nc in ncs
            r = Float64[]
            for n in ns
                f = "$(basef)-$(n)-$(edge)-$(res)$(nc).csv"
                A = readdlm(f, ',', Float64; skipstart=1)
                push!(r, corner[3] == :b ? A[1, 2] : A[end, 2])
            end
            extrapolation = nothing
            try
                extrapolation = RichardsonExtrapolationUQ.richextrapol_uq(r, 1.0 ./ ns)# richextrapol(r, 1.0 ./ ns)
            catch
                extrapolation = (data = r,)
            end
            push!(results, (corner = corner[1], nc = nc, extrapolation = extrapolation))
        end
    end
    return results
end

const MESH = :biased
const SUPPORT = :hard

function test_t3ff(ns=[4,], visualize=true)
    formul = FEMMShellT3FFModule
    @info "Scordelis-Lo shell, formulation=$(formul)"
    results = []
    for n in ns
        v = _execute_t3ff(formul, n, visualize)
        push!(results, v/(ANALYT_SOL)*100)
    end
    return ns, results
end

function test_q4rs(ns=[16, 32, 64, 128, 256, 512])
    @info "Scordelis-Lo shell, formulation=Q4RS"
    results = []
    for n in ns
        v = _execute_q4rs(MESH, n, SUPPORT)
        push!(results, v)
    end
    return ns, results
end


function allrun()
    # println("#####################################################")
    # println("# test_t3ff ")
    # test_t3ff()
    println("#####################################################")
    println("# test_q4rs ")
    test_q4rs()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
