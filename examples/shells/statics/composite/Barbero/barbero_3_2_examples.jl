"""
From Barbero's Finite Element Analysis using Abaqus ... book Example 3.2

Use Abaqus to model a simply supported rectangular plate with dimensions
ax = 4, 000 mm, ay = 2, 000 mm, and plate thickness t = 10 mm.
Apply a uniform
transverse load q0 = 0.12×10−3 MPa. The material is a unidirectional AS4D/9310 lamina
(Table 3.1), with the fibers oriented in the x-direction. Determine the deflection of the center
point of the plate. Use (1.102–1.104) to calculate the orthotropic stiffness matrix and enter
the material as Type:
Orthotropic. This example is continued in Example 3.8, p. 151.
The maximum deflection is 17.43 mm at the center of the plate.
"""
module barbero_3_2_examples

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
    # ASFD/9310
    E1 = 133860 * phun("MPa")
    E2 = 7706 * phun("MPa")
    G12 = 4306 * phun("MPa")
    G13 = G12
    nu12 = 0.301
    nu23 = 0.396
    G23 = 2760 * phun("MPa")
    ax = 4.0 * phun("m")
    ay = 2.0 * phun("m")
    aspect_ratio = 0.001
    thickness = 10.0 * phun("mm")
    q = 0.12e-3 * phun("mega*Pa")
    # For 14 elements we should get transverse deflection -1.75776e+01
    nx = ny = 14
    ply_thickness = thickness
    tolerance = ax / nx / 100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[CM.Ply("p", mater, ply_thickness, 0),]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("example_3.2", plies, mcsys)

    fens, fes = T3block(ax, ay, nx, ny)
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
    l1 = connectednodes(meshboundary(fes))
    for i in [1, 2, 3, ]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
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
        @show maximum(dchi.values[:, i]) / phun("mm")
        @show minimum(dchi.values[:, i]) / phun("mm")
    end
    
    vtkwrite("example_3.2-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])
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
