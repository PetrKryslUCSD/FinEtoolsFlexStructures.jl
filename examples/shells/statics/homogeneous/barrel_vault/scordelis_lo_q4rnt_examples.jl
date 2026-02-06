module scordelis_lo_q4rnt

using Test
using LinearAlgebra
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellQ4RNTModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

function _execute(n = 8, visualize = true)
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    analyt_sol = -0.3024
    # Parameters:
    E = 4.32e8
    nu = 0.0
    thickness = 0.25 # geometrical dimensions are in feet
    R = 25.0
    L = 50.0
    formul = FEMMShellQ4RNTModule

    tolerance = R / n / 1000
    fens, fes = Q4block(40 / 360 * 2 * pi, L / 2, n, n)
    fens.xyz = xyz3(fens)
    for i = 1:count(fens)
        a = fens.xyz[i, 1]
        y = fens.xyz[i, 2]
        fens.xyz[i, :] .= (R * sin(a), y, R * (cos(a) - 1))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, CompositeRule(GaussRule(2, 2), GaussRule(2, 1)), thickness), mater)
    # femm.mult_el_size = 0.2
    femm.drilling_stiffness_scale = 1.0
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
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1, 3, 5]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf L / 2 L / 2 -Inf Inf], inflate = tolerance)
    for i in [2, 4, 6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1, 5, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Midpoint of the free edge
    nl = selectnode(
        fens;
        box = Float64[sin(40 / 360 * 2 * pi) * 25 sin(40 / 360 * 2 * pi) * 25 L / 2 L / 2 -Inf Inf],
        inflate = tolerance,
    )
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(Float64[0, 0, -90, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 2)
    @show sum(F)

    # Solve
    solve_blocked!(dchi, K, F)
    @show resultpercent = dchi.values[nl, 3][1] / analyt_sol * 100

    # Visualization
    if visualize
        update_rotation_field!(Rfield0, dchi)
        vtkwrite("scordelis_lo_q4rnt-$(n)-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])

        # # Generate a graphical display of resultants
        # scalars = []
        # for nc in 1:3
        #     fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
        #     push!(scalars, ("m$nc", fld.values))
        #     fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
        #     push!(scalars, ("em$nc", fld.values))
        # end
        # vtkwrite("scordelis_lo_q4rnt-$(n)-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        # scalars = []
        # for nc in 1:3
        #     fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
        #     push!(scalars, ("n$nc", fld.values))
        #     fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
        #     push!(scalars, ("en$nc", fld.values))
        # end
        # vtkwrite("scordelis_lo_q4rnt-$(n)-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        # scalars = []
        # for nc in 1:2
        #     fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
        #     push!(scalars, ("q$nc", fld.values))
        #     fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
        #     push!(scalars, ("eq$nc", fld.values))
        # end
        # vtkwrite("scordelis_lo_q4rnt-o-$(n)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

        

        
    end

    return resultpercent
end

function test_convergence()
    formul = FEMMShellQ4RNTModule
    results = [
        66.54771615057949,
        85.54615143134853,
        89.85075281481419,
        92.50616661644985,
        95.40469210310079,
        97.65376880486126,
    ]
    for (n, res) in zip([4, 8, 10, 12, 16, 24], results)
        v = _execute(n, true)
        # @test isapprox(res, v, rtol = 1.0e-4)
        # @show v
    end
    return true
end

test_convergence()

end # module