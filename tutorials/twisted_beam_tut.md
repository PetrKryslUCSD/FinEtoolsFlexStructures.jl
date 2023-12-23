Static analysis of clamped twisted beam

Source code: [`twisted_beam_tut.jl`](twisted_beam_tut.jl)

Last updated: 12/23/23

## Description

The initially twisted cantilever beam is one of the standard test
problems for verifying the finite-element accuracy [1]. The beam is
clamped at one end and loaded either with unit in-plane or
unit out-of-plane force at the other. The centroidal axis of the beam is
straight at the undeformed  configuration, while its cross-sections are
twisted about the centroidal axis from 0 at the clamped end to pi/2 at
the free end.

Reference deflection in the direction of the applied force.
| Cross section thickness       | Loading in the Z direction | Loading in the Y direction |
| -------   |     -------  |  ------- |
| t = 0.32 | 0.005425   |    0.001753  |
| t = 0.0032 |  0.005256  |  0.001294    |

References:

[1] MacNeal,  R. H., and R. L. Harder, “A Proposed Standard Set of Problems to
Test Finite Element Accuracy,” Finite Elements in Analysis Design, vol. 11, pp.
3–20, 1985.

[2] Simo,  J. C., D. D. Fox, and M. S. Rifai, “On a Stress Resultant
Geometrically Exact Shell Model. Part II: The Linear Theory; Computational
Aspects,” Computational Methods in Applied Mechanical Engineering, vol. 73,
pp. 53–92, 1989.

[3] Zupan D, Saje M (2004) On "A proposed standard set of problems to test
finite element accuracy": the twisted beam. Finite Elements in Analysis
and Design 40: 1445-1451.

## Goals

# - Introduce definition of the shell model.
# - Calculate the discrete model quantities and solve the static equilibrium problem.
# - Demonstrate visualization of the resultant section forces and moments.

## Definition of the basic inputs

The finite element code relies on the basic functionality implemented in this
package.

````julia
using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
````

The inputs are defined in consistent units.
The elastic properties are:

````julia
E = 0.29e8;
nu = 0.22;
````

The material is elastic isotropic.

````julia
mater = MatDeforElastIso(DeforModelRed3D, E, nu)
````

Here are the cross-section width and the cantilevered length of the beam.

````julia
W = 1.1;
L = 12.0;
````

## Reference Solutions

The reference solutions are defined by these tuples of values: thickness,
magnitude of the force, direction of the force, and the reference deflection
along the force.

````julia
params_thicker_dir_3 = (t = 0.32, force = 1.0, dir = 3, uex = 0.005424534868469);
params_thicker_dir_2 = (t = 0.32, force = 1.0, dir = 2, uex = 0.001753248285256);

params_thinner_dir_3 = (t = 0.0032, force = 1.0e-6, dir = 3, uex = 0.005256);
params_thinner_dir_2 = (t = 0.0032, force = 1.0e-6, dir = 2, uex = 0.001294);
````

Now select one particular simulation. Here we go with the thicker shell and
the direction of the load 3:

````julia
params = params_thicker_dir_3
````

## Mesh generation

The mesh is initially generated for a rectangular 2d domain. The element size
can be controlled with these two variables:

````julia
nL, nW = 8, 2
````

This represents a fairly coarse mesh. We may try something more refined:
nL, nW = 4*8, 4*2

````julia
fens, fes = T3block(L, W, nL, nW, :a);
````

The 2d mesh is expanded into a three dimensional domain, and the
locations of the nodes are modified to produce the pre-youtwisted shape.

````julia
fens.xyz = xyz3(fens)
fens = let
    for i in 1:count(fens)
        x = fens.xyz[i, 1]
        a = x / L * (pi / 2)
        y = fens.xyz[i, 2] - (W / 2)
        z = fens.xyz[i, 3]
        fens.xyz[i, :] = [x, y * cos(a) - z * sin(a), y * sin(a) + z * cos(a)]
    end
    fens
end
````

The implementation of the 3-node triangle shell element is in this module. We
will refer to the functions that we need from this module by referencing them
relative to the module name.

````julia
t3ffm = FEMMShellT3FFModule
````

The shell elements have a type `FESetShellT3`. The type of the mesh is the
plain isoparametric triangle, which is instructed to delegate to the shell
element.

````julia
accepttodelegate(fes, FESetShellT3())
````

Use the convenience function to make an instance of the finite element model
machine for the T3FF shell.

````julia
femm = t3ffm.make(IntegDomain(fes, TriRule(1), params.t), mater)
````

Construct the requisite fields, geometry and displacement.
Initialize configuration variables. Displacements are all zero,
the rotation matrices are all identities.

````julia
geom0 = NodalField(fens.xyz)
u0 = NodalField(zeros(size(fens.xyz, 1), 3))
Rfield0 = initial_Rfield(fens)
dchi = NodalField(zeros(size(fens.xyz, 1), 6))
````

Apply supports. The clamped end is selected based on the X coordinate.
We will select all nodes within a box, and the box is slightly inflated using
the geometrical tolerance based on the spacing of the nodes.

````julia
tolerance = min(W / nW, L / nL) / 100

l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
for i in 1:6
    setebc!(dchi, l1, true, i)
end
numberdofs!(dchi);
````

Associate the finite element model machine with geometry. The shell
formulation requires knowledge of the normals to the shell surface, and that
information is computed from the positions of the nodes.

````julia
t3ffm.associategeometry!(femm, geom0)
````

Assemble the system stiffness matrix.

````julia
K = t3ffm.stiffness(femm, geom0, u0, Rfield0, dchi);
````

The load is a concentrated force applied at the center of the beam. First
we select the node.

````julia
nl = selectnode(fens; box = Float64[L L 0 0 0 0], tolerance = tolerance)
````

Next we create a mesh of the loaded boundary (i.e. the selected node), so that
we can integrate along it.

````julia
loadbdry = FESetP1(reshape(nl, 1, 1))
````

Now we create a finite element model machine that knows how to integrate.
Since our boundary consists of a single point, we use a point rule.

````julia
lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
````

Now we create a force intensity to represent the loading.

````julia
v = float.([0, 0, 0, 0, 0, 0])
v[params.dir] = params.force
fi = ForceIntensity(v);
````

Finally we computed the load vector corresponding to the force intensity.

````julia
F = distribloads(lfemm, geom0, dchi, fi, 3);
````

Extract the free-free block of the matrix, and the free block for the
right-hand side vector.

````julia
K_ff = matrix_blocked(K, nfreedofs(dchi))[:ff]
F_f = vector_blocked(F, nfreedofs(dchi))[:f]
````

The system of linear algebraic equations of balance is solved for the displacements and rotations.

````julia
U = K_ff \ F_f
````

The vector of unknowns is now distributed into a field.

````julia
scattersysvec!(dchi, U[:])
````

The deflection in the correct direction at the loaded node is now extracted.

````julia
tipdefl = dchi.values[nl, params.dir][1]
@info "Normalized deflection: $(tipdefl / params.uex * 100)%"
````

Generate a graphical display of resultants. The resultants only make sense in
a coordinate system aligned with the shell surface. Here we construct such a
coordinate system by using the tangent directions to the surface to construct
the normal, and the vector along the global X axis will serve as the first
basis vector of the local cartesian coordinate system.

````julia
function updatecsys!(csmatout, XYZ, tangents, feid, qpid)
````

first the normal

````julia
    cross3!(view(csmatout, :, 3), view(tangents, :, 1), view(tangents, :, 2))
    csmatout[:, 3] ./= norm(view(csmatout, :, 3))
    csmatout[:, 1] .= (1.0, 0.0, 0.0)
    cross3!(view(csmatout, :, 2), view(csmatout, :, 3), view(csmatout, :, 1))
    return csmatout
end
ocsys = CSys(3, 3, updatecsys!)
scalars = []
for nc in 1:3
    fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
````

fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)

````julia
    push!(scalars, ("m$nc", fld.values))
end
vtkwrite("twisted_beam-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
scalars = []
for nc in 1:3
    fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
````

fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)

````julia
    push!(scalars, ("n$nc", fld.values))
end
vtkwrite("twisted_beam-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
scalars = []
for nc in 1:2
    fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
````

fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)

````julia
    push!(scalars, ("q$nc", fld.values))
end
vtkwrite("twisted_beam-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

true
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

