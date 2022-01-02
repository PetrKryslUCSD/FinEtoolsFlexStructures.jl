module mcompshell0
# From Barbero's Finite Element Analysis using Abaqus ... book Example 3.1
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    
    ax = ay = 2000*phun("mm")
    nx = ny = 8
    # ASFD/9310
    E1 = 133860*phun("MPa")
    E2 = 7706*phun("MPa")
    G12 = 4306*phun("MPa")
    G13 = G12
    nu12 = 0.301;
    nu23 = 0.396
    G23 = 2760*phun("MPa")
    npairs = 1
    thickness = 10*phun("mm");
    tolerance = ax/nx/100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[]
    for p in 1:npairs
        push!(plies, CM.Ply("ply_0_$p", mater, thickness/npairs/2, 0))
        push!(plies, CM.Ply("ply_90_$p", mater, thickness/npairs/2, 90))
    end
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("example_3.1", plies, mcsys)

    fens, fes = T3block(ax,ay,nx,ny);
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
    # Pin one of the corners
    l1 = selectnode(fens; box = Float64[ax ax 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,2, 3,]
        setebc!(dchi, l1, true, i)
    end
    # Roller at the other
    l1 = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], inflate = tolerance)
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

    # Edge load
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, box = Float64[-Inf Inf 0 0 0 0 ], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, 0.1*phun("MPa")*thickness, 0, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 1);
    l1 = selectelem(fens, bfes, box = Float64[-Inf Inf ay ay 0 0 ], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, -0.1*phun("MPa")*thickness, 0, 0, 0, 0]);
    F += distribloads(lfemm, geom0, dchi, fi, 1);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    for i in 1:3
        # @show  maximum(dchi.values[:, i]) ./phun("mm")
        # @show minimum(dchi.values[:, i]) ./phun("mm")
    end
    # vtkwrite("plate-$npairs-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])
    @test maximum(dchi.values[:, 3]) ./phun("mm") ≈ 0.22004349767718365
    true
end
end
using .mcompshell0
mcompshell0.test()

module mcompshell1
# From Barbero's Finite Element Analysis using Abaqus ... book Example 3.1
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    
    ax = ay = 2000*phun("mm")
    nx = ny = 8
    # ASFD/9310
    E1 = 133860*phun("MPa")
    E2 = 7706*phun("MPa")
    G12 = 4306*phun("MPa")
    G13 = G12
    nu12 = 0.301;
    nu23 = 0.396
    G23 = 2760*phun("MPa")
    npairs = 10
    thickness = 10*phun("mm");
    tolerance = ax/nx/100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[]
    for p in 1:npairs
        push!(plies, CM.Ply("ply_0_$p", mater, thickness/npairs/2, 0))
        push!(plies, CM.Ply("ply_90_$p", mater, thickness/npairs/2, 90))
    end
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("example_3.1", plies, mcsys)
    
    fens, fes = T3block(ax,ay,nx,ny);
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
    # Pin one of the corners
    l1 = selectnode(fens; box = Float64[ax ax 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,2, 3,]
        setebc!(dchi, l1, true, i)
    end
    # Roller at the other
    l1 = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], inflate = tolerance)
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

    # Edge load
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, box = Float64[-Inf Inf 0 0 0 0 ], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, 0.1*phun("MPa")*thickness, 0, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 1);
    l1 = selectelem(fens, bfes, box = Float64[-Inf Inf ay ay 0 0 ], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, -0.1*phun("MPa")*thickness, 0, 0, 0, 0]);
    F += distribloads(lfemm, geom0, dchi, fi, 1);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    for i in 1:3
        # @show  maximum(dchi.values[:, i]) ./phun("mm")
        # @show minimum(dchi.values[:, i]) ./phun("mm")
    end
    # vtkwrite("plate-$npairs-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])
    @test maximum(dchi.values[:, 3]) ./phun("mm") ≈ 0.010276258430770539
    true
end
end
using .mcompshell1
mcompshell1.test()

module mcompshell2
# From Barbero's Finite Element Analysis using Abaqus ... book Example 3.3
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    
    ax = ay = 2000*phun("mm")
    nx = ny = 8
    # ASFD/9310
    E1 = 37880*phun("MPa")
    E2 = 9407*phun("MPa")
    G12 = 3405*phun("MPa")
    G13 = G12
    nu12 = 0.299;
    nu23 = 0.422
    G23 = 3308*phun("MPa")
    ply_thickness = 1*phun("mm");
    tolerance = ax/nx/100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[
    CM.Ply("p", mater, ply_thickness, +45),
    CM.Ply("p", mater, ply_thickness, -45),
    CM.Ply("p", mater, ply_thickness, 0),
    CM.Ply("p", mater, ply_thickness, 0),
    CM.Ply("p", mater, ply_thickness, -45),
    CM.Ply("p", mater, ply_thickness, +45),
    ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("example_3.3", plies, mcsys)

    fens, fes = T3block(ax,ay,nx,ny);
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)

    femm = formul.make(IntegDomain(fes, TriRule(1), length(plies)*ply_thickness), layup)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Pin one of the corners
    l1 = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], inflate = tolerance)
    for i in [1, 2, 3,]
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

    # Edge load
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, box = Float64[0 0 -Inf Inf 0 0 ], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[-200*phun("N/mm"), 0, 0, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 1);
    l1 = selectelem(fens, bfes, box = Float64[ax ax -Inf Inf 0 0 ], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[200*phun("N/mm"), 0, 0, 0, 0, 0]);
    F += distribloads(lfemm, geom0, dchi, fi, 1);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    for i in 1:3
        # @show maximum(dchi.values[:, i]) ./phun("mm")
        # @show minimum(dchi.values[:, i]) ./phun("mm")
    end
    # vtkwrite("plate-3.3-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])

     # The book: The target value is 1.667 mm at either of the loaded edges; so
     # the total displacement is approximately 2*1.667=3.334. The Abaqus model
     # on the other hand gives 1.658*2 = 3.316.
    @test maximum(dchi.values[:, 1]) ./phun("mm") ≈ 3.3161098653904544
    true
end
end
using .mcompshell2
mcompshell2.test()

module mcompshell3
# Example  from Section 6.17.1
# Lazarus Teneketzis Tenek, John Argyris, Finite Element Analysis for Composite Structures, 1998
# Clamped isotropic plate under uniform transverse pressure
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test_composite()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    E = 200*phun("GPa");
    nu = 0.3;
    L = 1.0*phun("m");
    t_radius_ratio = 0.001
    thickness = L * t_radius_ratio;
    q = 1.0*phun("kilo*Pa");
    D = E*thickness^3/12/(1-nu^2)
    wcref = 0.00126*q*L^4/D
    nx = ny = 4
    ply_thickness = thickness;
    tolerance = L/nx/100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E, nu)
    plies = CM.Ply[
    CM.Ply("p", mater, ply_thickness, 0),
    ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("Tenek, Argyris 6.17.1", plies, mcsys)

    fens, fes = T3block(L,L,nx,ny);
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)

    femm = formul.make(IntegDomain(fes, TriRule(1), CM.thickness(layup)), layup)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped Condition
    l1 = connectednodes(meshboundary(fes))
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, geom0, u0, Rfield0, dchi);

    # Edge load
    lfemm = FEMMBase(IntegDomain(fes, TriRule(1)))
    fi = ForceIntensity(FFlt[0, 0, -q, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    for i in 1:3
        # @show maximum(dchi.values[:, i])
        # @show minimum(dchi.values[:, i])
    end
    # vtkwrite("plate-6.17.1-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])

     
    @test isapprox(-minimum(dchi.values[:, 3]), wcref, rtol=0.02)
    true
end,
function test_homogeneous()
    formul = FEMMShellT3FFModule
    E = 200*phun("GPa");
    nu = 0.3;
    L = 1.0*phun("m");
    t_radius_ratio = 0.001
    thickness = L * t_radius_ratio;
    q = 1.0*phun("kilo*Pa");
    D = E*thickness^3/12/(1-nu^2)
    wcref = 0.00126*q*L^4/D
    nx = ny = 4
    ply_thickness = thickness;
    tolerance = L/nx/100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E, nu)

    fens, fes = T3block(L,L,nx,ny);
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)

    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped Condition
    l1 = connectednodes(meshboundary(fes))
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, geom0, u0, Rfield0, dchi);

    # Edge load
    lfemm = FEMMBase(IntegDomain(fes, TriRule(1)))
    fi = ForceIntensity(FFlt[0, 0, -q, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    for i in 1:3
        # @show maximum(dchi.values[:, i])
        # @show minimum(dchi.values[:, i])
    end
    # vtkwrite("plate-6.17.1-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])

     
    @test isapprox(-minimum(dchi.values[:, 3]), wcref, rtol=0.02)
    true
end
end
using .mcompshell3
mcompshell3.test_composite()
mcompshell3.test_homogeneous()

module mcompshell4
# Example  from Section 6.17.1
# Lazarus Teneketzis Tenek, John Argyris, Finite Element Analysis for Composite Structures, 1998
# Clamped isotropic plate under uniform transverse pressure
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test_composite()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    E = 200*phun("GPa");
    nu = 0.3;
    L = 1.0*phun("m");
    t_radius_ratio = 0.001
    thickness = L * t_radius_ratio;
    q = 1.0*phun("kilo*Pa");
    D = E*thickness^3/12/(1-nu^2)
    wcref = 0.00126*q*L^4/D
    nx = ny = 4
    ply_thickness = thickness/3;
    tolerance = L/nx/100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E, nu)
    plies = CM.Ply[
    CM.Ply("p", mater, ply_thickness, 0),
    CM.Ply("p", mater, ply_thickness, 45),
    CM.Ply("p", mater, ply_thickness, 90),
    ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("Tenek, Argyris 6.17.1", plies, mcsys)

    fens, fes = T3block(L,L,nx,ny);
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)

    femm = formul.make(IntegDomain(fes, TriRule(1), CM.thickness(layup)), layup)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped Condition
    l1 = connectednodes(meshboundary(fes))
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, geom0, u0, Rfield0, dchi);

    # Edge load
    lfemm = FEMMBase(IntegDomain(fes, TriRule(1)))
    fi = ForceIntensity(FFlt[0, 0, -q, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    for i in 1:3
        # @show maximum(dchi.values[:, i])
        # @show minimum(dchi.values[:, i])
    end
    # vtkwrite("plate-6.17.1-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])

     
    @test isapprox(-minimum(dchi.values[:, 3]), wcref, rtol=0.02)
    true
end

end
using .mcompshell4
mcompshell4.test_composite()


module mcompshell5
# Example  from Section 6.17.2
# Lazarus Teneketzis Tenek, John Argyris, Finite Element Analysis for Composite Structures, 1998
# Simply supported sandwich plate.
# The deflection provided in the textbook of 31.45 mm is probably wrong.
# The deflection of 2.026 mm was obtained with Abaqus.
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
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
    Ef = 6.8*phun("GPa");
    nuf = 0.3;
    Ec = 480*phun("MPa");
    nuc = 0.3
    L = 0.1384*phun("m");
    tf = 0.001*phun("m");
    tc = 0.015*phun("m");
    q = 1.0*phun("MPa");
    nx = ny = 140
    tolerance = L/nx/100
    CM = CompositeLayupModule

    plies = CM.Ply[
    CM.Ply("face_1", CM.lamina_material(Ef, nuf), tf, 0),
    CM.Ply("core", CM.lamina_material(Ec, nuc), tc, 0),
    CM.Ply("face_2", CM.lamina_material(Ef, nuf), tf, 0),
    ]
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("Tenek, Argyris 6.17.2", plies, mcsys)

    fens, fes = T3block(L,L,nx,ny);
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)

    femm = formul.make(IntegDomain(fes, TriRule(1), CM.thickness(layup)), layup)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Pin one of the corners
    l1 = selectnode(fens; box = Float64[L L 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,2, 3,]
        setebc!(dchi, l1, true, i)
    end
        # Roller at the other
    l1 = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], inflate = tolerance)
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

    # Face load
    lfemm = FEMMBase(IntegDomain(fes, TriRule(1)))
    fi = ForceIntensity(FFlt[0, 0, -q, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    # for i in 1:3
    #     @show maximum(dchi.values[:, i]) ./phun("mm")
    #     @show minimum(dchi.values[:, i]) ./phun("mm")
    # end
    # vtkwrite("plate-6.17.2-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])

     
    @test isapprox(-minimum(dchi.values[:, 3]), 2.026e-3, rtol=0.02)
    true
end

end
using .mcompshell5
mcompshell5.test()

module mcompshell6
# Similar definition of the plate as in:
# From Barbero's Finite Element Analysis using Abaqus ... book Example 3.1
# but the loading is transverse, pressure of magnitude 1.0
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    
    ax = ay = 1000*phun("mm") # note the changed dimension
    nx = ny = 50
    # ASFD/9310
    E1 = 133860*phun("MPa")
    E2 = 7706*phun("MPa")
    G12 = 4306*phun("MPa")
    G13 = G12
    nu12 = 0.301;
    nu23 = 0.396
    G23 = 2760*phun("MPa")
    q = 1.0*phun("MPa")
    npairs = 1
    thickness = 10*phun("mm");
    # With these inputs, the Abaqus-verified solution is -1502.0324040447633

    tolerance = ax/nx/100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[]
    for p in 1:npairs
        push!(plies, CM.Ply("ply_0_$p", mater, thickness/npairs/2, 0))
        push!(plies, CM.Ply("ply_90_$p", mater, thickness/npairs/2, 90))
    end
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("example_3.1", plies, mcsys)

    fens, fes = T3block(ax,ay,nx,ny);
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
    l1 = connectednodes(meshboundary(fes))
    for i in [1, 2, 3]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, geom0, u0, Rfield0, dchi);

    lfemm = FEMMBase(IntegDomain(fes, TriRule(1)))
    fi = ForceIntensity(FFlt[0, 0, -q, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    for i in 1:3
        # @show  maximum(dchi.values[:, i]) ./phun("mm")
        # @show minimum(dchi.values[:, i]) ./phun("mm")
    end
    # vtkwrite("plate-$npairs-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])

    @test minimum(dchi.values[:, 3]) ./phun("mm") ≈ -1502.0324040447633
    true
end
end
using .mcompshell6
mcompshell6.test()


module mcompshell7
# Similar definition of the plate as in mcompshell6 but the layup is [-45/45]
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    
    ax = ay = 1000*phun("mm") # note the changed dimension
    nx = ny = 50
    # ASFD/9310
    E1 = 133860*phun("MPa")
    E2 = 7706*phun("MPa")
    G12 = 4306*phun("MPa")
    G13 = G12
    nu12 = 0.301;
    nu23 = 0.396
    G23 = 2760*phun("MPa")
    q = 1.0*phun("MPa")
    npairs = 1
    thickness = 10*phun("mm");
    # With these inputs, the Abaqus-verified solution is -1176.8346659735575

    tolerance = ax/nx/100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[]
    for p in 1:npairs
        push!(plies, CM.Ply("ply_-45_$p", mater, thickness/npairs/2, -45))
        push!(plies, CM.Ply("ply_+45_$p", mater, thickness/npairs/2, +45))
    end
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("example_3.1", plies, mcsys)

    fens, fes = T3block(ax,ay,nx,ny);
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
    l1 = connectednodes(meshboundary(fes))
    for i in [1, 2, 3]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, geom0, u0, Rfield0, dchi);

    lfemm = FEMMBase(IntegDomain(fes, TriRule(1)))
    fi = ForceIntensity(FFlt[0, 0, -q, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 2);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    for i in 1:3
        # @show  maximum(dchi.values[:, i]) ./phun("mm")
        # @show minimum(dchi.values[:, i]) ./phun("mm")
    end
    # vtkwrite("plate-$npairs-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])

    @test minimum(dchi.values[:, 3]) ./phun("mm") ≈ -1176.8346659735575
    true
end
end
using .mcompshell7
mcompshell7.test()

module mcompshell8
# From Barbero's Finite Element Analysis using Abaqus ... book Example 3.1
# Variation which tests the extension-bending coupling: unsymmetric cross-ply
# layup should bend one way under compression, and another way under tension.
# Two plies, 90 degrees at the bottom. Applied tensile load in the zero
# direction (the reference direction of the composite layup) should bend the
# plate downwards (in the negative z direction). And vice versa.

using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test(load_multiplier)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    
    ax = ay = 2000*phun("mm")
    nx = ny = 8
    # ASFD/9310
    E1 = 133860*phun("MPa")
    E2 = 7706*phun("MPa")
    G12 = 4306*phun("MPa")
    G13 = G12
    nu12 = 0.301;
    nu23 = 0.396
    G23 = 2760*phun("MPa")
    thickness = 10*phun("mm");
    tolerance = ax/nx/100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[]
    push!(plies, CM.Ply("ply_90", mater, thickness/2, 90))
    push!(plies, CM.Ply("ply_0", mater, thickness/2, 0))
    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("example_3.1", plies, mcsys)

    fens, fes = T3block(ax,ay,nx,ny);
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
    # Pin one of the corners
    l1 = selectnode(fens; box = Float64[ax ax 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,2, 3,]
        setebc!(dchi, l1, true, i)
    end
    # Roller at the other
    l1 = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], inflate = tolerance)
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

    # Edge load
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, box = Float64[0 0 -Inf Inf 0 0 ], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0.1*phun("MPa")*load_multiplier*thickness, 0, 0, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 1);
    l1 = selectelem(fens, bfes, box = Float64[ax ax -Inf Inf 0 0 ], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[-0.1*phun("MPa")*load_multiplier*thickness, 0, 0, 0, 0, 0]);
    F += distribloads(lfemm, geom0, dchi, fi, 1);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    for i in 1:3
        # @show maximum(dchi.values[:, i]) ./phun("mm")
        # @show minimum(dchi.values[:, i]) ./phun("mm")
    end
    # vtkwrite("plate-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])
    
    if load_multiplier < 0
        @test minimum(dchi.values[:, 3]) ./phun("mm") ≈ -0.22004349767718365
    else
        @test maximum(dchi.values[:, 3]) ./phun("mm") ≈ 0.22004349767718365
    end
    true
end
end
using .mcompshell8
mcompshell8.test(+1)
mcompshell8.test(-1)

module mcompshell9
# As mcompshell8, but the plate is oriented in the yz plane, with the loading in
# the z direction.

using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test(load_multiplier)
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    
    ax = ay = 2000*phun("mm")
    nx = ny = 8
    # ASFD/9310
    E1 = 133860*phun("MPa")
    E2 = 7706*phun("MPa")
    G12 = 4306*phun("MPa")
    G13 = G12
    nu12 = 0.301;
    nu23 = 0.396
    G23 = 2760*phun("MPa")
    thickness = 10*phun("mm");
    tolerance = ax/nx/100
    CM = CompositeLayupModule

    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[]
    push!(plies, CM.Ply("ply_90", mater, thickness/2, 90))
    push!(plies, CM.Ply("ply_0", mater, thickness/2, 0))
    mcsys = CSys(Float64[
        0 0 1; 
        1 0 0; 
        0 1 0])
    layup = CM.CompositeLayup("example_3.1", plies, mcsys)

    fens, fes = T3block(ax,ay,nx,ny);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        x, y, z = fens.xyz[i, :]
        fens.xyz[i, :] .= (z, x, y)
    end
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
    # Pin one of the corners
    l1 = selectnode(fens; box = Float64[0 0 0 0 0 0], inflate = tolerance)
    for i in [1,2, 3,]
        setebc!(dchi, l1, true, i)
    end
    # Roller at the other
    l1 = selectnode(fens; box = Float64[0 0 0 0 ax ax], inflate = tolerance)
    for i in [2,]
        setebc!(dchi, l1, true, i)
    end
    # Simple support
    l1 = connectednodes(meshboundary(fes))
    for i in [1]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, geom0, u0, Rfield0, dchi);

    # Edge load
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, box = Float64[0 0 -Inf Inf ax ax], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, 0, 0.1*phun("MPa")*load_multiplier*thickness, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 1);
    l1 = selectelem(fens, bfes, box = Float64[0 0 -Inf Inf 0 0], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, 0, -0.1*phun("MPa")*load_multiplier*thickness, 0, 0, 0]);
    F += distribloads(lfemm, geom0, dchi, fi, 1);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    for i in 1:3
        # @show maximum(dchi.values[:, i]) ./phun("mm")
        # @show minimum(dchi.values[:, i]) ./phun("mm")
    end
    # vtkwrite("plate-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])
    
    if load_multiplier < 0
        @test minimum(dchi.values[:, 1]) ./phun("mm") ≈ -0.22004349767718365
    else
        @test maximum(dchi.values[:, 1]) ./phun("mm") ≈ 0.22004349767718365
    end
    true
end
end
using .mcompshell9
mcompshell9.test(+1)
mcompshell9.test(-1)

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

function test()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    
    a = 250*phun("mm") # note the changed dimension
    nx = ny = 8
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

    # @show maximum(dchi.values[:, 3]) * (100*E2) / (q*thickness*(a/thickness)^4)
    # @show maximum(dchi.values[:, 1]) * (100*E2) / (q*thickness*(a/thickness)^4)
    @test maximum(dchi.values[:, 3]) * (100*E2) / (q*thickness*(a/thickness)^4) ≈ 1.9529295928169146
    @test maximum(dchi.values[:, 1]) * (100*E2) / (q*thickness*(a/thickness)^4) ≈ 0.12375916912652352
    
    # vtkwrite("plate--uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])

    # @test minimum(dchi.values[:, 3]) ./phun("mm") ≈ -1176.8346659735575
    true
end
end
using .mcompshell10
mcompshell10.test()