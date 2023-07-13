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

    deflex=5*norm(q)*L^4/(384*E*I3);

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
    for i in [1,2,3,4,5] # cylindrical joint
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 0 +L/2 +L/2 0 0], tolerance = L/10000)
    for i in [1,2,3,4,5] # cylindrical joint
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
    # @show dchi.values[tipl, 1], deflex
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

    deflex=5*norm(q)*L^4/(384*E*I2);

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

    fi = ForceIntensity(FFlt[0, 0, q]);
    F = distribloads_global(femm, geom0, u0, Rfield0, dchi, fi);

    # Solve the static problem
    solve!(dchi, K, F)
    # @show dchi.values[tipl, 3], deflex
    @test norm(dchi.values[tipl, 3] .- deflex) / deflex < 1.0e-5

    true
end
test()
end # module

module mccbeam1
#  Clamped-clamped beam with distributed load. Soft axis bending

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
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 0 +L/2 +L/2 0 0], tolerance = L/10000)
    for i in [1,2,3,4,5,6] # clamped
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
    # @show dchi.values[tipl, 1], deflex
    @test norm(dchi.values[tipl, 1] .- deflex) / deflex < 1.0e-5


    true
end
test()
end # module

module mccbeam2
#  Clamped-clamped beam with distributed load. Stiff axis bending

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
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 0 +L/2 +L/2 0 0], tolerance = L/10000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    fi = ForceIntensity(FFlt[0, 0, q]);
    F = distribloads_global(femm, geom0, u0, Rfield0, dchi, fi);

    # Solve the static problem
    solve!(dchi, K, F)
    # @show dchi.values[tipl, 3], deflex
    @test norm(dchi.values[tipl, 3] .- deflex) / deflex < 1.0e-5

    true
end
test()
end # module


module margyrislframe1
##  Argyris L-frame, in-plane load

##
# References
#
# J.H. Argyris, P.C. Dunne and D.W. Scharpf, On large
# displacement-small strain analysis of structures with
# rotation degree of freedom. Comput. Methods Appl. Mech.
# Engrg. 14 (1978) 401-451; 15 (1978) 99-135.
#
# J.H. Argyris, 0. Hilpert, GA. Malejannakis and D.W. Scharpf, On the
# geometrical stiffness of a beam in space-A consistent V.W.
# Approach, Comput. Methods Appl. Mech. Engrg. 20 (1979) 105-131.
#
#
# Verified with continuum element solution with high order hexahedron.
# Displacement of the point of application of the load
#  1.0e-005 *[-0.1831  -0.0714  -0.0000]
# Note the orientation of the structure is different here. The obtained
# displacements correlate in the obvious way.

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
    E=71240.0;# MPa
    nu=0.31;# Poisson ratio
    b=0.6; h=30; L=240;# cross-sectional dimensions and length of each leg in millimeters
    magn=+1e-5;# Magnitude of the total applied force, Newton
    deflex =  [-1.91840e-06 0.00000e+00 -7.18697e-07]

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 1.0, 0.0])

    # Select the number of elements per leg.
    n = 8;
    members = []
    push!(members, frame_member([0 0 L; L 0 L], n, cs))
    push!(members, frame_member([L 0 L; L 0 0], n, cs))
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
    clampedn = selectnode(fens; box = Float64[0 0 0 0 L L], tolerance = L/10000)
    tipn = selectnode(fens; box = Float64[L L 0 0 0 0], tolerance = L/10000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, clampedn, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    loadbdry = FESetP1(reshape(tipn, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[-magn, 0, 0, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)

    # @show dchi.values[tipn, 1:3], deflex
    @test norm(vec(dchi.values[tipn, 1:3]) - vec(deflex)) / norm(deflex) < 1.0e-5

    true
end
test()
end # module

module margyrislframe2
##  Argyris L-frame, out-of-plane load

##
# References
#
# J.H. Argyris, P.C. Dunne and D.W. Scharpf, On large
# displacement-small strain analysis of structures with
# rotation degree of freedom. Comput. Methods Appl. Mech.
# Engrg. 14 (1978) 401-451; 15 (1978) 99-135.
#
# J.H. Argyris, 0. Hilpert, GA. Malejannakis and D.W. Scharpf, On the
# geometrical stiffness of a beam in space-A consistent V.W.
# Approach, Comput. Methods Appl. Mech. Engrg. 20 (1979) 105-131.
#
#
# Verified with continuum element solution with high order hexahedron.
# Displacement of the point of application of the load in the direction of the
# force 0.0047, 0.0 in the remaining two directions. Note the orientation of
# the structure is different here. The obtained displacements correlate in the
# obvious way.

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
    E=71240.0;# MPa
    nu=0.31;# Poisson ratio
    b=0.6; h=30; L=240;# cross-sectional dimensions and length of each leg in millimeters
    magn=+1e-5;# Magnitude of the total applied force, Newton
    deflex =  [0 4.78846e-03 0]

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 1.0, 0.0])

    # Select the number of elements per leg.
    n = 8;
    members = []
    push!(members, frame_member([0 0 L; L 0 L], n, cs))
    push!(members, frame_member([L 0 L; L 0 0], n, cs))
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
    clampedn = selectnode(fens; box = Float64[0 0 0 0 L L], tolerance = L/10000)
    tipn = selectnode(fens; box = Float64[L L 0 0 0 0], tolerance = L/10000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, clampedn, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    loadbdry = FESetP1(reshape(tipn, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, magn, 0, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)

    # @show dchi.values[tipn, 1:3], deflex
    @test norm(vec(dchi.values[tipn, 1:3]) - vec(deflex)) / norm(deflex) < 1.0e-5

    true
end
test()
end # module

nothing

