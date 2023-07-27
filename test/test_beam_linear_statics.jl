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

"""
CURVED BEAM WITH STATIC LOADS
PROBLEM DESCRIPTION
In this example a curved cantilever beam, modeled with shell elements, is
subjected to unit forces at the tip in the in-plane and out-of-plane directions, that
is, the Y and Z directions, respectively. The in-plane and out-of-plane loads are
applied in different load cases. The tip displacements in the direction of the load
are compared with independent hand calculated results.
The geometry, properties and loading are as suggested in MacNeal and Harder
1985. The cantilever beam is bent into a 90° arc. It has a 4.12 inch inner radius
and a 4.32 inch outer radius. Thus it is 0.2 inch wide and approximately 6.63 inch
long at its centerline. The beam is 0.1 inch thick in the Y direction. For modeling
in SAP2000, the curved beam is meshed into six area objects, each subtending a
15° arc.

Material Properties
E = 10,000,000 lb/in2
ν = 0.25
G = 4,000,000 lb/in2
Section Properties
Thickness = 0. 1 in
"""
module mcurvedbeam1

using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMRITBeamModule
using FinEtoolsFlexStructures.FEMMRITBeamModule: FEMMRITBeam
stiffness = FEMMRITBeamModule.stiffness
distribloads_global = FEMMRITBeamModule.distribloads_global
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
# # using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function test()
    E = 10000000.0 #  lb/in2
    nu = 0.25
    # Section Properties
    b = Thickness = 0.1 # in
    h = Depth = 0.2 # in
    radius = (4.32 + 4.12) / 2 # in
    k_s = 5/6 # shear correction factor
    force = 1.0 # lb
    # Select the number of elements per leg.
    nel = 4
    direction = 2
    deflex = 0.0886
    # direction = 3
    # deflex = 0.5

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s)

    ang=90
    members = []
    xyz=fill(0.0, 2, 3)
    for i in 1:nel
        xyz[1, :] .= (radius*cos((i-1)/nel*ang/360*2*pi), radius*sin((i-1)/nel*ang/360*2*pi), 0)
        xyz[2, :] .= (radius*cos((i)/nel*ang/360*2*pi), radius*sin((i)/nel*ang/360*2*pi), 0)
        push!(members, frame_member(xyz, 1, cs))
    end
    fens, fes = merge_members(members; tolerance = radius / 1000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    tipl = selectnode(fens; box = Float64[0 0 radius radius 0 0], inflate = radius/1000)
    basel = selectnode(fens; box = Float64[radius radius 0 0 0 0], inflate = radius/1000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, basel, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMRITBeam(IntegDomain(fes, GaussRule(1, 1)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    loadbdry = FESetP1(reshape(tipl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = fill(0.0, 6); q[direction] = force
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)
    # @show dchi.values[tipl, direction], deflex
    @test norm(dchi.values[tipl, direction] .- deflex) / deflex < 0.05

    # scaling = 1e0
    # dchi.values .*= scaling
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_space_box([[-radius -radius -radius]; [radius radius radius]]),
    #     plot_nodes(fens),
    #     plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    # pl = render(plots)

    true
end
test()
end # module


"""
CURVED BEAM WITH STATIC LOADS
PROBLEM DESCRIPTION
In this example a curved cantilever beam, modeled with shell elements, is
subjected to unit forces at the tip in the in-plane and out-of-plane directions, that
is, the Y and Z directions, respectively. The in-plane and out-of-plane loads are
applied in different load cases. The tip displacements in the direction of the load
are compared with independent hand calculated results.
The geometry, properties and loading are as suggested in MacNeal and Harder
1985. The cantilever beam is bent into a 90° arc. It has a 4.12 inch inner radius
and a 4.32 inch outer radius. Thus it is 0.2 inch wide and approximately 6.63 inch
long at its centerline. The beam is 0.1 inch thick in the Y direction. For modeling
in SAP2000, the curved beam is meshed into six area objects, each subtending a
15° arc.

Material Properties
E = 10,000,000 lb/in2
ν = 0.25
G = 4,000,000 lb/in2
Section Properties
Thickness = 0. 1 in
"""
module mcurvedbeam2

using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMRITBeamModule
using FinEtoolsFlexStructures.FEMMRITBeamModule: FEMMRITBeam
stiffness = FEMMRITBeamModule.stiffness
distribloads_global = FEMMRITBeamModule.distribloads_global
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
# # using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function test()
    E = 10000000.0 #  lb/in2
    nu = 0.25
    # Section Properties
    b = Thickness = 0.1 # in
    h = Depth = 0.2 # in
    radius = (4.32 + 4.12) / 2 # in
    k_s = 5/6 # shear correction factor
    force = 1.0 # lb
    # Select the number of elements per leg.
    nel = 4
    # direction = 2
    # deflex = 0.0886
    direction = 3
    deflex = 0.5

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s)

    ang=90
    members = []
    xyz=fill(0.0, 2, 3)
    for i in 1:nel
        xyz[1, :] .= (radius*cos((i-1)/nel*ang/360*2*pi), radius*sin((i-1)/nel*ang/360*2*pi), 0)
        xyz[2, :] .= (radius*cos((i)/nel*ang/360*2*pi), radius*sin((i)/nel*ang/360*2*pi), 0)
        push!(members, frame_member(xyz, 1, cs))
    end
    fens, fes = merge_members(members; tolerance = radius / 1000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    tipl = selectnode(fens; box = Float64[0 0 radius radius 0 0], inflate = radius/1000)
    basel = selectnode(fens; box = Float64[radius radius 0 0 0 0], inflate = radius/1000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, basel, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMRITBeam(IntegDomain(fes, GaussRule(1, 1)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    loadbdry = FESetP1(reshape(tipl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = fill(0.0, 6); q[direction] = force
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)
    # @show dchi.values[tipl, direction], deflex
    @test norm(dchi.values[tipl, direction] .- deflex) / deflex < 0.05

    # scaling = 1e0
    # dchi.values .*= scaling
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_space_box([[-radius -radius -radius]; [radius radius radius]]),
    #     plot_nodes(fens),
    #     plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    # pl = render(plots)

    true
end
test()
end # module

"""
CURVED BEAM WITH STATIC LOADS
PROBLEM DESCRIPTION
In this example a curved cantilever beam, modeled with shell elements, is
subjected to unit forces at the tip in the in-plane and out-of-plane directions, that
is, the Y and Z directions, respectively. The in-plane and out-of-plane loads are
applied in different load cases. The tip displacements in the direction of the load
are compared with independent hand calculated results.
The geometry, properties and loading are as suggested in MacNeal and Harder
1985. The cantilever beam is bent into a 90° arc. It has a 4.12 inch inner radius
and a 4.32 inch outer radius. Thus it is 0.2 inch wide and approximately 6.63 inch
long at its centerline. The beam is 0.1 inch thick in the Y direction. For modeling
in SAP2000, the curved beam is meshed into six area objects, each subtending a
15° arc.

Material Properties
E = 10,000,000 lb/in2
ν = 0.25
G = 4,000,000 lb/in2
Section Properties
Thickness = 0. 1 in
"""
module mcurvedbeam3

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
# # using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function test()
    E = 10000000.0 #  lb/in2
    nu = 0.25
    # Section Properties
    b = Thickness = 0.1 # in
    h = Depth = 0.2 # in
    radius = (4.32 + 4.12) / 2 # in
    k_s = 5/6 # shear correction factor
    force = 1.0 # lb
    # Select the number of elements per leg.
    nel = 4
    direction = 2
    deflex = 0.0886
    # direction = 3
    # deflex = 0.5

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s)

    ang=90
    members = []
    xyz=fill(0.0, 2, 3)
    for i in 1:nel
        xyz[1, :] .= (radius*cos((i-1)/nel*ang/360*2*pi), radius*sin((i-1)/nel*ang/360*2*pi), 0)
        xyz[2, :] .= (radius*cos((i)/nel*ang/360*2*pi), radius*sin((i)/nel*ang/360*2*pi), 0)
        push!(members, frame_member(xyz, 1, cs))
    end
    fens, fes = merge_members(members; tolerance = radius / 1000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    tipl = selectnode(fens; box = Float64[0 0 radius radius 0 0], inflate = radius/1000)
    basel = selectnode(fens; box = Float64[radius radius 0 0 0 0], inflate = radius/1000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, basel, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 1)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    loadbdry = FESetP1(reshape(tipl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = fill(0.0, 6); q[direction] = force
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)
    # @show dchi.values[tipl, direction], deflex
    @test norm(dchi.values[tipl, direction] .- deflex) / deflex < 0.05

    # scaling = 1e0
    # dchi.values .*= scaling
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_space_box([[-radius -radius -radius]; [radius radius radius]]),
    #     plot_nodes(fens),
    #     plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    # pl = render(plots)

    true
end
test()
end # module

"""
CURVED BEAM WITH STATIC LOADS
PROBLEM DESCRIPTION
In this example a curved cantilever beam, modeled with shell elements, is
subjected to unit forces at the tip in the in-plane and out-of-plane directions, that
is, the Y and Z directions, respectively. The in-plane and out-of-plane loads are
applied in different load cases. The tip displacements in the direction of the load
are compared with independent hand calculated results.
The geometry, properties and loading are as suggested in MacNeal and Harder
1985. The cantilever beam is bent into a 90° arc. It has a 4.12 inch inner radius
and a 4.32 inch outer radius. Thus it is 0.2 inch wide and approximately 6.63 inch
long at its centerline. The beam is 0.1 inch thick in the Y direction. For modeling
in SAP2000, the curved beam is meshed into six area objects, each subtending a
15° arc.

Material Properties
E = 10,000,000 lb/in2
ν = 0.25
G = 4,000,000 lb/in2
Section Properties
Thickness = 0. 1 in
"""
module mcurvedbeam4

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
# # using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function test()
    E = 10000000.0 #  lb/in2
    nu = 0.25
    # Section Properties
    b = Thickness = 0.1 # in
    h = Depth = 0.2 # in
    radius = (4.32 + 4.12) / 2 # in
    k_s = 5/6 # shear correction factor
    force = 1.0 # lb
    # Select the number of elements per leg.
    nel = 4
    # direction = 2
    # deflex = 0.0886
    direction = 3
    deflex = 0.5

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s)

    ang=90
    members = []
    xyz=fill(0.0, 2, 3)
    for i in 1:nel
        xyz[1, :] .= (radius*cos((i-1)/nel*ang/360*2*pi), radius*sin((i-1)/nel*ang/360*2*pi), 0)
        xyz[2, :] .= (radius*cos((i)/nel*ang/360*2*pi), radius*sin((i)/nel*ang/360*2*pi), 0)
        push!(members, frame_member(xyz, 1, cs))
    end
    fens, fes = merge_members(members; tolerance = radius / 1000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    tipl = selectnode(fens; box = Float64[0 0 radius radius 0 0], inflate = radius/1000)
    basel = selectnode(fens; box = Float64[radius radius 0 0 0 0], inflate = radius/1000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, basel, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 1)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    loadbdry = FESetP1(reshape(tipl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = fill(0.0, 6); q[direction] = force
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)
    # @show dchi.values[tipl, direction], deflex
    @test norm(dchi.values[tipl, direction] .- deflex) / deflex < 0.05

    # scaling = 1e0
    # dchi.values .*= scaling
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_space_box([[-radius -radius -radius]; [radius radius radius]]),
    #     plot_nodes(fens),
    #     plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    # pl = render(plots)

    true
end
test()
end # module
nothing

nothing


module mccbeamThin1
#  Clamped-clamped beam with distributed load. Soft axis bending

using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMRITBeamModule
using FinEtoolsFlexStructures.FEMMRITBeamModule: FEMMRITBeam
stiffness = FEMMRITBeamModule.stiffness
distribloads_global = FEMMRITBeamModule.distribloads_global
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
# # using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function test(b = 2)
    E=30002*1000.0;#30002ksi
    # b=2;#in
    h=18;#in
    L=240;# in
    q=1/b;# 1200 lbf/ft

    # Cross-sectional properties
    A=b*h;#in^2
    I2=b*h^3/12;#cm^4
    I3=b^3*h/12;#cm^4
    k_s = 5/6
    nu=0.0;

    ##
    # Exact deformation under the load

    deflex=norm(q)*L^4/(384*E*I3);

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [-1.0, 0.0, 0.0], k_s)

    # Select the number of elements per leg.
    n = 16;
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
    tipl = selectnode(fens; box = Float64[0 0 0 0 0 0], inflate = L/10000)
    l1 = selectnode(fens; box = Float64[0 0 -L/2 -L/2 0 0], inflate = L/10000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 0 +L/2 +L/2 0 0], inflate = L/10000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMRITBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    fi = ForceIntensity(FFlt[q, 0, 0]);
    F = distribloads_global(femm, geom0, u0, Rfield0, dchi, fi);

    # Solve the static problem
    solve!(dchi, K, F)
    return dchi.values[tipl, 1][1] / deflex * 100, b/L
    # @test norm(dchi.values[tipl, 1] .- deflex) / deflex < 1.0e-5


    # scaling = 1e4
    # dchi.values .*= scaling
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_space_box([[-L/2 -L/2 -L/2]; [L/2 L/2 L/2]]),
    # plot_nodes(fens),
    # plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    # pl = render(plots)
    # true
end
for b in [2  0.02 0.0002  ]
    @test abs(test(b)[1] - 98.44) / 100 <= 1e-2
end
end # module


module mccbeamThin2
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
# # using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function test(b = 2)
    E=30002*1000.0;#30002ksi
    # b=2;#in
    h=18;#in
    L=240;# in
    q=1/b;# 1200 lbf/ft

    # Cross-sectional properties
    A=b*h;#in^2
    I2=b*h^3/12;#cm^4
    I3=b^3*h/12;#cm^4
    k_s = 5/6
    nu=0.0;

    ##
    # Exact deformation under the load

    deflex=norm(q)*L^4/(384*E*I3);

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [-1.0, 0.0, 0.0], k_s)

    # Select the number of elements per leg.
    n = 16;
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
    tipl = selectnode(fens; box = Float64[0 0 0 0 0 0], inflate = L/10000)
    l1 = selectnode(fens; box = Float64[0 0 -L/2 -L/2 0 0], inflate = L/10000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 0 +L/2 +L/2 0 0], inflate = L/10000)
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
    return dchi.values[tipl, 1][1] / deflex * 100, b/L
    # @test norm(dchi.values[tipl, 1] .- deflex) / deflex < 1.0e-5


    # scaling = 1e4
    # dchi.values .*= scaling
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_space_box([[-L/2 -L/2 -L/2]; [L/2 L/2 L/2]]),
    # plot_nodes(fens),
    # plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    # pl = render(plots)
    # true
end
for b in [2  0.02 0.0002  ]
    # @show test(b)
    @test abs(test(b)[1] - 100) / 100 <= 1e-2
end
end # module

nothing



"""
CURVED BEAM WITH STATIC LOADS
PROBLEM DESCRIPTION
In this example a curved cantilever beam, modeled with shell elements, is
subjected to unit forces at the tip in the in-plane and out-of-plane directions, that
is, the Y and Z directions, respectively. The in-plane and out-of-plane loads are
applied in different load cases. The tip displacements in the direction of the load
are compared with independent hand calculated results.
The geometry, properties and loading are as suggested in MacNeal and Harder
1985. The cantilever beam is bent into a 90° arc. It has a 4.12 inch inner radius
and a 4.32 inch outer radius. Thus it is 0.2 inch wide and approximately 6.63 inch
long at its centerline. The beam is 0.1 inch thick in the Y direction. For modeling
in SAP2000, the curved beam is meshed into six area objects, each subtending a
15° arc.

Material Properties
E = 10,000,000 lb/in2
ν = 0.25
G = 4,000,000 lb/in2
Section Properties
Thickness = 0. 1 in
"""
module mcurvedbeam2a

using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMLinBeamModule
using FinEtoolsFlexStructures.FEMMLinBeamModule: FEMMLinBeam
stiffness = FEMMLinBeamModule.stiffness
distribloads_global = FEMMLinBeamModule.distribloads_global
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
# using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function test(direction = 2)
    E = 10000000.0 #  lb/in2
    nu = 0.25
    # Section Properties
    b = Thickness = 0.1 # in
    h = Depth = 0.2 # in
    radius = (4.32 + 4.12) / 2 # in
    k_s = 5/6 # shear correction factor
    force = 1.0 # lb
    # Select the number of elements per leg.
    nel = 4
    if direction == 2
        deflex = 0.0886
    else
    # direction = 3
        deflex = 0.5
    end

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s)

    ang=90
    members = []
    xyz=fill(0.0, 2, 3)
    for i in 1:nel
        xyz[1, :] .= (radius*cos((i-1)/nel*ang/360*2*pi), radius*sin((i-1)/nel*ang/360*2*pi), 0)
        xyz[2, :] .= (radius*cos((i)/nel*ang/360*2*pi), radius*sin((i)/nel*ang/360*2*pi), 0)
        push!(members, frame_member(xyz, 1, cs))
    end
    fens, fes = merge_members(members; tolerance = radius / 1000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    tipl = selectnode(fens; box = Float64[0 0 radius radius 0 0], inflate = radius/1000)
    basel = selectnode(fens; box = Float64[radius radius 0 0 0 0], inflate = radius/1000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, basel, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMLinBeam(IntegDomain(fes, GaussRule(1, 1)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    loadbdry = FESetP1(reshape(tipl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = fill(0.0, 6); q[direction] = force
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)
    # @show dchi.values[tipl, direction], deflex
    @test norm(dchi.values[tipl, direction] .- deflex) / deflex < 0.05

    # scaling = 1e0
    # dchi.values .*= scaling
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_space_box([[-radius -radius -radius]; [radius radius radius]]),
    #     plot_nodes(fens),
    #     plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    # pl = render(plots)

    true
end
test(2)
test(3)
end # module

nothing



"""
CURVED BEAM WITH STATIC LOADS
PROBLEM DESCRIPTION
In this example a curved cantilever beam, modeled with shell elements, is
subjected to unit forces at the tip in the in-plane and out-of-plane directions, that
is, the Y and Z directions, respectively. The in-plane and out-of-plane loads are
applied in different load cases. The tip displacements in the direction of the load
are compared with independent hand calculated results.
The geometry, properties and loading are as suggested in MacNeal and Harder
1985. The cantilever beam is bent into a 90° arc. It has a 4.12 inch inner radius
and a 4.32 inch outer radius. Thus it is 0.2 inch wide and approximately 6.63 inch
long at its centerline. The beam is 0.1 inch thick in the Y direction. For modeling
in SAP2000, the curved beam is meshed into six area objects, each subtending a
15° arc.

Material Properties
E = 10,000,000 lb/in2
ν = 0.25
G = 4,000,000 lb/in2
Section Properties
Thickness = 0. 1 in
"""
module mcurvedbeam1a

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
# using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function test(direction = 2)
    E = 10000000.0 #  lb/in2
    nu = 0.25
    # Section Properties
    b = Thickness = 0.1 # in
    h = Depth = 0.2 # in
    radius = (4.32 + 4.12) / 2 # in
    k_s = 5/6 # shear correction factor
    force = 1.0 # lb
    # Select the number of elements per leg.
    nel = 4
        if direction == 2
        deflex = 0.0886
    else
        # direction = 3
        deflex = 0.5
    end

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s)

    ang=90
    members = []
    xyz=fill(0.0, 2, 3)
    for i in 1:nel
        xyz[1, :] .= (radius*cos((i-1)/nel*ang/360*2*pi), radius*sin((i-1)/nel*ang/360*2*pi), 0)
        xyz[2, :] .= (radius*cos((i)/nel*ang/360*2*pi), radius*sin((i)/nel*ang/360*2*pi), 0)
        push!(members, frame_member(xyz, 1, cs))
    end
    fens, fes = merge_members(members; tolerance = radius / 1000);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    tipl = selectnode(fens; box = Float64[0 0 radius radius 0 0], inflate = radius/1000)
    basel = selectnode(fens; box = Float64[radius radius 0 0 0 0], inflate = radius/1000)
    for i in [1,2,3,4,5,6] # clamped
        setebc!(dchi, basel, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 1)), material)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    loadbdry = FESetP1(reshape(tipl, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = fill(0.0, 6); q[direction] = force
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)
    # @show dchi.values[tipl, direction], deflex
    @test norm(dchi.values[tipl, direction] .- deflex) / deflex < 0.05

    # scaling = 1e0
    # dchi.values .*= scaling
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_space_box([[-radius -radius -radius]; [radius radius radius]]),
    #     plot_nodes(fens),
    #     plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    # pl = render(plots)

    true
end
test(2)
test(3)
end # module


"""
gravity loading on a straight beam
Inconsistent loading with distribloads
"""
module mgravitybeamincons2

using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMLinBeamModule
using FinEtoolsFlexStructures.FEMMLinBeamModule: FEMMLinBeam
stiffness = FEMMLinBeamModule.stiffness
inspectintegpoints = FEMMLinBeamModule.inspectintegpoints
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
# using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function test(nel = 2)
    E = 12.0*phun("GPa")
    nu = 0.3
    # Section Properties
    b = Thickness = 0.2*phun("m")
    h = Depth = 0.4*phun("m")
    I = b * h^3 / 12
    ecc = 0.10*phun("m")
    L = 10.0*phun("m")
    k_s = 5/6 # shear correction factor
    w = 2000.0*phun("kg/m^3") * 10.0*phun("m/sec^2")
    uniform_eccentricity = [0.0, 0.0, 0.0, ecc]
    # uniform_eccentricity = [0.0, 0.0, 0.0*phun("m"), 0.0]

    # Cross-sectional properties
    # cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s) # Timoshenko
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0]) # Bernoulli

    xyz = [[L/2 0 0]; [-L/2 0 0]]
    fens, fes = frame_member(xyz, nel, cs)
    # @show count(fes)

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    leftl = selectnode(fens; box = initbox!([], xyz[1, :]), inflate = L/1000)
    rightl = selectnode(fens; box = initbox!([], xyz[2, :]), inflate = L/1000)
    for i in [1,2,3,4,5] # pin
        setebc!(dchi, leftl, true, i)
        setebc!(dchi, rightl, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMLinBeam(IntegDomain(fes, GaussRule(1, 1), b * h), material, uniform_eccentricity)

    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    q = fill(0.0, 6); q[2] = w
    fi = ForceIntensity(q)
    F = distribloads(femm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)
# @show dchi.values

    K_df = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:df]
    F_d =  K_df * gathersysvec(dchi, :f)
    H =  F_d[dchi.dofnums[leftl[1], 1]-nfreedofs(dchi)]
    M0 = H * ecc
    deflw = (5 * w * b * h * L^4) / (384 * E * I)
    deflH = (M0 * L^2) / (8 * E * I)

    @test abs(maximum(dchi.values[:, 2]) - 1.02307e-02) / 1.02307e-02 < 1.0e-3

    function inspector(idat, i, conn, ecoords, elvecf, loc)
        if conn[1] == leftl[1] || conn[2] == leftl[1]
            # @info "element at support $(xyz[1, :])"
            # @show elvecf
            @test abs(elvecf[1]) ≈ H
        end
    end

    inspectintegpoints(femm, geom0, dchi, NodalField([1.0]), 1:count(fes), inspector, nothing)

    # scaling = 2e2
    # dchi.values .*= scaling
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_space_box([[-L -L -L]; [L L L]]),
    #     plot_nodes(fens),
    #     plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    # pl = render(plots)

    true
end
test(4)
end # module

nothing



"""
gravity loading on a straight beam
Inconsistent loading with distribloads
Alternative definition of the cross section orientation
"""
module mgravitybeamincons3

using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMLinBeamModule
using FinEtoolsFlexStructures.FEMMLinBeamModule: FEMMLinBeam
stiffness = FEMMLinBeamModule.stiffness
inspectintegpoints = FEMMLinBeamModule.inspectintegpoints
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
# using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function test(nel = 2)
    E = 12.0*phun("GPa")
    nu = 0.3
    # Section Properties: The beam is being bent in the f1-f3 plane
    b = 0.4*phun("m")
    h = 0.2*phun("m")
    I = b^3 * h / 12
    ecc = 0.10*phun("m")
    L = 10.0*phun("m")
    k_s = 5/6 # shear correction factor
    w = 2000.0*phun("kg/m^3") * 10.0*phun("m/sec^2")
    uniform_eccentricity = [0.0, 0.0, ecc, 0.0]
    # uniform_eccentricity = [0.0, 0.0, 0.0*phun("m"), 0.0]

    # Cross-sectional properties
    # cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0], k_s) # Timoshenko
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 1.0, 0.0]) # Bernoulli

    xyz = [[L/2 0 0]; [-L/2 0 0]]
    fens, fes = frame_member(xyz, nel, cs)
    # @show count(fes)

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    leftl = selectnode(fens; box = initbox!([], xyz[1, :]), inflate = L/1000)
    rightl = selectnode(fens; box = initbox!([], xyz[2, :]), inflate = L/1000)
    for i in [1,2,3,4,5] # pin
        setebc!(dchi, leftl, true, i)
        setebc!(dchi, rightl, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMLinBeam(IntegDomain(fes, GaussRule(1, 1), b * h), material, uniform_eccentricity)

    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    q = fill(0.0, 6); q[2] = w
    fi = ForceIntensity(q)
    F = distribloads(femm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve!(dchi, K, F)
# @show dchi.values

    K_df = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:df]
    F_d =  K_df * gathersysvec(dchi, :f)
    H =  F_d[dchi.dofnums[leftl[1], 1]-nfreedofs(dchi)]
    M0 = H * ecc
    deflw = (5 * w * b * h * L^4) / (384 * E * I)
    deflH = (M0 * L^2) / (8 * E * I)
    @test abs(maximum(dchi.values[:, 2]) - 1.02307e-02) / 1.02307e-02 < 1.0e-3

    function inspector(idat, i, conn, ecoords, elvecf, loc)
        if conn[1] == leftl[1] || conn[2] == leftl[1]
            # @info "element at support $(xyz[1, :])"
            # @show elvecf
            @test abs(elvecf[1]) ≈ H
        end
    end

    inspectintegpoints(femm, geom0, dchi, NodalField([1.0]), 1:count(fes), inspector, nothing)

    # scaling = 2e2
    # dchi.values .*= scaling
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_space_box([[-L -L -L]; [L L L]]),
    #     plot_nodes(fens),
    #     plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values); dims = 1)
    # pl = render(plots)

    true
end
test(4)
end # module



nothing
