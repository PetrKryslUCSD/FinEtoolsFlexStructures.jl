module mcompshelldyn0

# Introduction to Finite Element Vibration Analysis 2nd Edition
# by Maurice Petyt 
#  Publisher: Cambridge University Press; 2nd edition (August 1, 2010)
# ISBN-10: 0521191602
# ISBN-13: 978-0521191609 

# Table 9.2. Non-dimensional natural frequencies of a simply supported
# square laminated plate (a/h = 10) [45/−45/45/−45]
# Natural frequencies, non dimensional
# 18.46, 34.87, 54.27, 67.17, 75.58, 82.84

using Arpack
using LinearAlgebra: norm, Transpose, mul!, I, Symmetric
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using Test

function test_homogeneous()
    formul = FEMMShellT3FFModule
    
    ax = ay = 1000*phun("mm")
    nx = ny = 8
    E = 1000.0*phun("MPa")
    nu = 0.396
    rho = 2.0*phun("kg/m^3")
    thickness = ax/1000;
    tolerance = ax/nx/100
    D = E*thickness^3/12/(1-nu^2)
    
    fens, fes = T3block(ax,ay,nx,ny);
    fens.xyz = xyz3(fens)
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)

    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # No in plane displacements
    # l1 = collect(1:count(fens))
    # for i in [1, 2,]
    #     setebc!(dchi, l1, true, i)
    # end
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
    M = formul.mass(femm, geom0, dchi);

    # Solve
    neigvs = 8
    d, v, nconv = eigs(Symmetric(K), Symmetric(M); nev=neigvs, which=:SM, explicittransform=:none)
    # @show nconv
    fs = real(sqrt.(complex(d))) / (2 * pi)
    # From Blevins, Table 5.3, square plate SSSS
    # for m in 1:3
    #     for n in 1:3
    #         @show pi^2 * (m^2 + n^2) * sqrt(D/rho/thickness) / (2 * pi * ax^2)
    #     end
    # end
        
    # vectors = []
    # for i in 1:neigvs
    #     scattersysvec!(dchi, v[:, i])
    #     push!(vectors, ("mode_$i", deepcopy(dchi.values[:, 1:3])))
    # end
    # vtkwrite("plate-modes.vtu", fens, fes; vectors = vectors)
    # @test maximum(dchi.values[:, 3]) ./phun("mm") ≈ 0.22004349767718365
    @test norm(fs - [21.822774909379287, 54.52203720717488, 54.60475772077202, 86.5480437899711, 108.65349048567792, 108.78050954575491, 138.29506871044424, 140.89370172063866]) < 1.0e-3 * norm(fs)
    true
end


function test_composite()
    formul = FEMMShellT3FFCompModule
    CM = CompositeLayupModule
    
    ax = ay = 1000*phun("mm")
    nx = ny = 8
    E = 1000.0*phun("MPa")
    nu = 0.396
    rho = 2.0*phun("kg/m^3")
    thickness = ax/1000;
    tolerance = ax/nx/100
    D = E*thickness^3/12/(1-nu^2)

    mater = CM.lamina_material(rho, E, nu)
    plies = CM.Ply[]
    push!(plies, CM.Ply("ply_1", mater, thickness/4, 45))
    push!(plies, CM.Ply("ply_2", mater, thickness/4, -45))
    push!(plies, CM.Ply("ply_3", mater, thickness/4, 45))
    push!(plies, CM.Ply("ply_4", mater, thickness/4, -45))

    mcsys = CM.cartesian_csys((1, 2, 3))
    layup = CM.CompositeLayup("petyt_9.2", plies, mcsys)

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
    for i in [1, 2, 3]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    formul.associategeometry!(femm, geom0)
    K = formul.stiffness(femm, geom0, u0, Rfield0, dchi);
    M = formul.mass(femm, geom0, dchi);

    # Solve
    neigvs = 8
    d, v, nconv = eigs(Symmetric(K), Symmetric(M); nev=neigvs, which=:SM, explicittransform=:none)
    # @show nconv
    fs = real(sqrt.(complex(d))) / (2 * pi)
    # oms = real(sqrt.(complex(d)))
    # @show (oms*ax^2*sqrt(rho/E/thickness^2))
    # vectors = []
    # for i in 1:neigvs
    #     scattersysvec!(dchi, v[:, i])
    #     push!(vectors, ("mode_$i", deepcopy(dchi.values[:, 1:3])))
    # end
    # vtkwrite("plate-modes.vtu", fens, fes; vectors = vectors)
    @test norm(fs - [21.822774909379287, 54.52203720717488, 54.60475772077202, 86.5480437899711, 108.65349048567792, 108.78050954575491, 138.29506871044424, 140.89370172063866]) < 1.0e-3 * norm(fs)
    true
end
end
using .mcompshelldyn0
mcompshelldyn0.test_homogeneous()
mcompshelldyn0.test_composite()

# module mcompshelldyn1

# # Introduction to Finite Element Vibration Analysis 2nd Edition
# # by Maurice Petyt 
# #  Publisher: Cambridge University Press; 2nd edition (August 1, 2010)
# # ISBN-10: 0521191602
# # ISBN-13: 978-0521191609 

# # Table 9.2. Non-dimensional natural frequencies of a simply supported
# # square laminated plate (a/h = 10) [45/−45/45/−45]
# # Natural frequencies, non dimensional
# # 18.46, 34.87, 54.27, 67.17, 75.58, 82.84

# using Arpack
# using LinearAlgebra: norm, Transpose, mul!, I, Symmetric
# using FinEtools
# using FinEtoolsDeforLinear
# using FinEtoolsFlexStructures.CompositeLayupModule
# using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
# using FinEtoolsFlexStructures.FEMMShellT3FFCompModule
# using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
# using FinEtools.MeshExportModule.VTKWrite: vtkwrite
# using Test

# function test()
#     formul = FEMMShellT3FFCompModule
#     CM = CompositeLayupModule
    
#     ax = ay = 1000*phun("mm")
#     nx = ny = 50
#     E1 = 1000.0*phun("MPa")
#     E2 = E1/40
#     G12 = E2*0.6
#     G13 = E2*0.5
#     nu12 = 0.25;
#     nu23 = 0.396
#     G23 = E2*0.5
#     rho = 10.0*phun("kg/m^3")
#     thickness = ax/10;
#     tolerance = ax/nx/100
#     CM = CompositeLayupModule

#     mater = CM.lamina_material(rho, E1, E2, nu12, G12, G13, G23)
#     plies = CM.Ply[]
#     push!(plies, CM.Ply("ply_1", mater, thickness/4, 45))
#     push!(plies, CM.Ply("ply_2", mater, thickness/4, -45))
#     push!(plies, CM.Ply("ply_3", mater, thickness/4, 45))
#     push!(plies, CM.Ply("ply_4", mater, thickness/4, -45))

#     mcsys = CM.cartesian_csys((1, 2, 3))
#     layup = CM.CompositeLayup("petyt_9.2", plies, mcsys)

#     fens, fes = T3block(ax,ay,nx,ny);
#     fens.xyz = xyz3(fens)
#     sfes = FESetShellT3()
#     accepttodelegate(fes, sfes)

#     femm = formul.make(IntegDomain(fes, TriRule(1), thickness), layup)

#     # Construct the requisite fields, geometry and displacement
#     # Initialize configuration variables
#     geom0 = NodalField(fens.xyz)
#     u0 = NodalField(zeros(size(fens.xyz,1), 3))
#     Rfield0 = initial_Rfield(fens)
#     dchi = NodalField(zeros(size(fens.xyz,1), 6))

#     # Apply EBC's
#     # Pin one of the corners
#     l1 = selectnode(fens; box = Float64[ax ax 0 0 -Inf Inf], inflate = tolerance)
#     for i in [1,2, 3,]
#         setebc!(dchi, l1, true, i)
#     end
#     # Roller at the other
#     l1 = selectnode(fens; box = Float64[0 0 0 0 -Inf Inf], inflate = tolerance)
#     for i in [2,]
#         setebc!(dchi, l1, true, i)
#     end
#     # Simple support
#     l1 = connectednodes(meshboundary(fes))
#     for i in [1, 2, 3]
#         setebc!(dchi, l1, true, i)
#     end
#     applyebc!(dchi)
#     numberdofs!(dchi);

#     # Assemble the system matrix
#     formul.associategeometry!(femm, geom0)
#     K = formul.stiffness(femm, geom0, u0, Rfield0, dchi);
#     M = formul.mass(femm, geom0, dchi);

#     # Solve
#     neigvs = 8
#     d, v, nconv = eigs(Symmetric(K), Symmetric(M); nev=neigvs, which=:SM, explicittransform=:none)
#     @show nconv
#     oms = real(sqrt.(complex(d)))
#     @show (oms*ax^2*sqrt(rho/E2/thickness^2))
#     vectors = []
#     for i in 1:neigvs
#         scattersysvec!(dchi, v[:, i])
#         push!(vectors, ("mode_$i", deepcopy(dchi.values[:, 1:3])))
#     end
#     vtkwrite("plate-modes.vtu", fens, fes; vectors = vectors)
#     # @test maximum(dchi.values[:, 3]) ./phun("mm") ≈ 0.22004349767718365
#     true
# end
# end
# using .mcompshelldyn1
# mcompshelldyn1.test()
