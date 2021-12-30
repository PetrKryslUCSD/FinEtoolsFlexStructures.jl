module mcompshell0
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using Test

function test()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    # From Barbero's Finite Element Analysis using Abaqus ... book Example 3.1
    ax = ay = 2000*phun("mm")
    nx = ny = 3
    # ASFD/9310
    E1 = 133860*phun("MPa")
    E2 = 7706*phun("MPa")
    G12 = 4306*phun("MPa")
    G13 = G12
    nu12 = 0.301;
    nu23 = 0.396
    G23 = 2760*phun("MPa")
    npairs = 2
    thickness = 10*phun("mm");
    tolerance = ax/nx/100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[]
    for p in npairs
        push!(plies, CM.Ply("ply_0_$p", mater, thickness, 0))
        push!(plies, CM.Ply("ply_90_$p", mater, thickness, 90))
    end
    mcsys = CM.cartesian_csys((1, 2, 3))
    cl = CM.CompositeLayup("example_3.1", plies, mcsys)

    fens, fes = T3block(ax,ay,nx,ny);
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)

    layups = CM.CompositeLayups([cl])
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), layups)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Pin one of the corners
    l1 = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,2, 3,]
        setebc!(dchi, l1, true, i)
    end
    # Roller at the other
    l1 = selectnode(fens; box = Float64[ax ax 0 0 -Inf Inf], inflate = tolerance)
    for i in [2,]
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

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[sin(40/360*2*pi)*25 sin(40/360*2*pi)*25 L/2 L/2 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt[0, 0, -90, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    resultpercent = dchi.values[nl, 3][1]/analyt_sol*100

    true
end
end
using .mcompshell0
mcompshell0.test()
