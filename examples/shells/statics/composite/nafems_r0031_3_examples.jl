
module nafems_r0031_3_examples
# From Barbero's Finite Element Analysis using Abaqus ... book Example 3.7
# Here we actually solve original benchmark:
# R0031(3): Three-layer sandwich shell under normal pressure loading
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test(N = 18)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule

    a = 2 * 127 * phun("mm")
    nx = ny = N

    E1 = 68947.57 * phun("MPa")
    E2 = 27579.03 * phun("MPa")
    G12 = G13 = 12927.67 * phun("MPa")
    nu12 = 0.3
    G23 = 12927.67 * phun("MPa")
    thickness_face = 0.7112 * phun("mm")
    face = CM.lamina_material(E1, E2, nu12, G12, G13, G23)

    E1 = 0.068948 * phun("MPa")
    E2 = 0.068948 * phun("MPa")
    nu12 = 0.0
    G12 = 0.068948 * phun("MPa")
    G13 = 209.6 * phun("MPa")
    G23 = 82.737 * phun("MPa")
    thickness_core = 19.05 * phun("mm")
    core = CM.lamina_material(E1, E2, nu12, G12, G13, G23)

    q = 0.689476 * phun("MPa")

    tolerance = a / nx / 100

    plies = CM.Ply[
        CM.Ply("ply_face_1", face, thickness_face, 0),
        CM.Ply("ply_core", core, thickness_core, 0),
        CM.Ply("ply_face_2", face, thickness_face, 0),
    ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("Barbero 3.7", plies, mcsys)

    fens, fes = T3block(a, a, nx, ny)
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
    # Simple support
    l1 = connectednodes(meshboundary(fes))
    for i in [1, 2, 3]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, geom0, u0, Rfield0, dchi)

    lfemm = FEMMBase(IntegDomain(fes, TriRule(1)))
    fi = ForceIntensity(Float64[0, 0, q, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 2)

    # Solve
    solve_blocked!(dchi, K, F)

    # R0031(3): Three-layer sandwich shell under normal pressure loading
    # lists 0.123'' = 3.1242 mm
    @show maxdisp = maximum(dchi.values[:, 3]) ./ phun("mm") 
    @show maxdisp = maximum(dchi.values[:, 3]) ./ phun("in")
    
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

