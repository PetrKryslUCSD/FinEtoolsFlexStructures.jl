module mssbeam1
#  Simply supported beam with distributed load. Soft axis bending

using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMCorotBeamModule
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
distribloads_global = FEMMCorotBeamModule.distribloads_global
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test

function test()
    E=30002*1000.0;#30002ksi
    b=2;#in
    h=18;#in
    L=240;# in
    q=1;# 1200 lbf/ft

    # Cross-sectional properties
    A=b*h;#in^2
    I2=b*h^3/12;#cm^4
    I3=b^3*h/12;#cm^4
    nu=0.0;

    ##
    # Exact deformation under the load

    deflex=norm(q)*L^4/(384*E*I3);

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [-1.0, 0.0, 0.0])

    # Select the number of elements per leg.
    n = 4;
    members = []
    push!(members, frame_member([0 -L/2 0; 0 L/2 0], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    tipl = selectnode(fens; box = Float64[0 0 0 0 0 0], tolerance = L/10000)
    l1 = selectnode(fens; box = Float64[0 0 -L/2 -L/2 0 0], tolerance = L/10000)
    for i in [1,2,3,5,6] # cylindrical joint
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 0 +L/2 +L/2 0 0], tolerance = L/10000)
    for i in [1,2,3,5,6] # cylindrical joint
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    fi = ForceIntensity(FFlt[q, 0, 0]);
    F = distribloads_global(femm, geom0, u0, Rfield0, dchi, fi);

    # Solve the static problem
    solve!(dchi, K, F)
    @test norm(dchi.values[tipl, 1] .- deflex) / deflex < 1.0e-5

    true
end
test()
end # module

module mssbeam2
#  Simply supported beam with distributed load. Stiff axis bending

using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMCorotBeamModule
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
distribloads_global = FEMMCorotBeamModule.distribloads_global
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test

function test()
    E=30002*1000.0;#30002ksi
    b=2;#in
    h=18;#in
    L=240;# in
    q=1;# 1200 lbf/ft

    # Cross-sectional properties
    A=b*h;#in^2
    I2=b*h^3/12;#cm^4
    I3=b^3*h/12;#cm^4
    nu=0.0;

    ##
    # Exact deformation under the load

    deflex=norm(q)*L^4/(384*E*I2);

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0])

    # Select the number of elements per leg.
    n = 4;
    members = []
    push!(members, frame_member([0 -L/2 0; 0 L/2 0], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    tipl = selectnode(fens; box = Float64[0 0 0 0 0 0], tolerance = L/10000)
    l1 = selectnode(fens; box = Float64[0 0 -L/2 -L/2 0 0], tolerance = L/10000)
    for i in [1,2,3,5,6] # cylindrical joint
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 0 +L/2 +L/2 0 0], tolerance = L/10000)
    for i in [1,2,3,5,6] # cylindrical joint
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    fi = ForceIntensity(FFlt[q, 0, 0]);
    F = distribloads_global(femm, geom0, u0, Rfield0, dchi, fi);

    # Solve the static problem
    solve!(dchi, K, F)
    @test norm(dchi.values[tipl, 1] .- deflex) / deflex < 1.0e-5

    true
end
test()
end # module

nothing
