module mcompshell10
# A mixed-enhanced finite-element for the analysis of laminated composite plates
# FERDINANDO AURICCHIO, ELIO SACCO

# Laminated [0/90] cross-ply square simply supported plate. Uniform Distributed
# Loading. Results in Table IV.

# Note: The scaling of the horizontal displacement of the midpoint of the edge
# is wrong in the paper. The power of the aspect ratio should also be 4.
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test(visualize = true)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    
    a = 250*phun("mm") # note the changed dimension
    nx = ny = 28
    # high modulus orthotropic graphite=epoxy composite material
    E1 = 100*phun("GPa")
    E2 = E1/25
    G12 = G13 = E2*0.5
    nu12 = 0.25;
    nu23 = 0.25
    G23 = E2*0.2
    q = 10.0*phun("kilo*Pa")
    thickness = a/10;
    # With these inputs, the Abaqus-verified solution is -1176.8346659735575

    tolerance = a/nx/100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[
        CM.Ply("ply_0", mater, thickness/2, 0),
        CM.Ply("ply_90", mater, thickness/2, 90)
        ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("sacco", plies, mcsys)

    fens, fes = T3block(a,a,nx,ny);
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)

    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), layup)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Simple support
    # The boundary conditions are a bit peculiar: refer to the paper
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf 0 0], inflate = tolerance)
    for i in [2, 4]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[a a -Inf Inf 0 0], inflate = tolerance)
    for i in [2, 4]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 0 0], inflate = tolerance)
    for i in [1, 5]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[-Inf Inf a a 0 0], inflate = tolerance)
    for i in [1, 5]
        setebc!(dchi, l1, true, i)
    end
    # Simple support
    l1 = connectednodes(meshboundary(fes))
    for i in [3]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, geom0, u0, Rfield0, dchi);

    lfemm = FEMMBase(IntegDomain(fes, TriRule(1)))
    fi = ForceIntensity(FFlt[0, 0, q, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])

    
    # Visualization
    if visualize
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = layup.csys)
            push!(scalars, ("m$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = layup.csys)
            push!(scalars, ("em$nc", fld.values))
        end
        vtkwrite("sacco_comp_pl_examples-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = layup.csys)
            push!(scalars, ("n$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = layup.csys)
            push!(scalars, ("en$nc", fld.values))
        end
        vtkwrite("sacco_comp_pl_examples-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = layup.csys)
            push!(scalars, ("q$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = layup.csys)
            push!(scalars, ("eq$nc", fld.values))
        end
        vtkwrite("sacco_comp_pl_examples-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

    end
    true
end

function allrun()
    println("#####################################################")
    println("# test_convergence ")
    test()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
