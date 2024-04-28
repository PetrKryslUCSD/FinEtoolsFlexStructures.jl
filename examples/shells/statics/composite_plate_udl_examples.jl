"""
Composite plate [0/90/90/0] with uniform distributed load.
Square shape.
"""
module composite_plate_udl_examples

using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    # E-glass/epoxy
    E1 = 60000 * phun("psi")
    E2 = 15000 * phun("psi")
    G13 = G12 = 6200 * phun("psi")
    nu12 = 0.28  
    G23 = 5000 * phun("psi")
    ax = 12.0 * phun("in")
    ay = 12.0 * phun("in")
    thickness = 0.2 * phun("in")
    q = 0.05 * phun("psi")
    nx = ny = 100
    tolerance = ax / nx / 100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[
        CM.Ply("p", mater, thickness/4, 0),
        CM.Ply("p", mater, thickness/4, 90),
        CM.Ply("p", mater, thickness/4, 90),
        CM.Ply("p", mater, thickness/4, 0),
        ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("composite_plate", plies, mcsys)

    fens, fes = T3block(ax, ay, nx, ny)
    bfes = meshboundary(fes)
    lx1 = selectelem(fens, bfes; facing = true, direction = [-1.0, 0.0])
    lx2 = selectelem(fens, bfes; facing = true, direction = [+1.0, 0.0])
    ly1 = selectelem(fens, bfes; facing = true, direction = [0.0, -1.0])
    ly2 = selectelem(fens, bfes; facing = true, direction = [0.0, +1.0])
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)

    femm = formul.make(IntegDomain(fes, TriRule(1), CM.thickness(layup)), layup)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # Simply supported Condition
    for el in [lx1, lx2, ly1, ly2]
        nl = connectednodes(subset(bfes, el))
        for d in [1, 2, 3]
            setebc!(dchi, nl, true, d)
        end
    end
    numberdofs!(dchi)

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, geom0, u0, Rfield0, dchi)

    # Edge load
    lfemm = FEMMBase(IntegDomain(fes, TriRule(1)))
    fi = ForceIntensity(Float64[0, 0, -q, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 2)

    # Solve
    solve_blocked!(dchi, K, F)
    for i = 1:3
        @show maximum(dchi.values[:, i]) / phun("in")
        @show minimum(dchi.values[:, i]) / phun("in")
    end
    
    vtkwrite("composite_plate-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])
    true
end

function allrun()
    println("#####################################################")
    println("# test ")
    test()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
