# Modal analysis of a free-floating steel circle: three-dimensional solid model

Source code: [`circle_3d_modal_conv_tut.jl`](circle_3d_modal_conv_tut.jl)

Last updated: 04/23/24

## Description

Vibration analysis of a free-floating steel ring. This is a
benchmark from the NAFEMS Selected Benchmarks for Natural Frequency Analysis,
publication: Test VM09: Circular Ring --  In-plane and Out-of-plane
Vibration.

The results can be compared with analytical expressions, but the main purpose
is to compute data for extrapolation to the limit to predict the true natural
frequencies.

## Reference frequencies

There will be 6 rigid body modes (zero natural frequencies).

The numerical results are due to the publication:
NAFEMS Finite Element Methods & Standards, Abbassian, F., Dawswell, D. J., and
Knowles, N. C. Selected Benchmarks for Natural Frequency Analysis, Test No.
6. Glasgow: NAFEMS, Nov., 1987.

The reference values were analytically determined (Blevins, FORMULAS FOR
DYNAMICS, ACOUSTICS AND VIBRATION, Table 4.16). Note that shear flexibility
was neglected when computing the reference values.

| Mode number | Shape                | Reference Value (Hz)  |  NAFEMS Target Value (Hz) |
| --------   |  -------------------  |  -------------------- |  ------------------------ |
| 7, 8       | (out of plane)        |        51.85          |           52.29           |
| 9, 10      |  (in plane)           |        53.38          |           53.97           |
| 11, 12     |  (out of plane)       |       148.8           |          149.7            |
| 13, 14     |  (in plane)           |       151.0           |          152.4            |
| 15, 16     |  (out of plane)       |       287.0           |          288.3            |
|  17, 18    |  (in plane)           |       289.5           |          288.3            |

## Goals

- Show convergence relative to reference values.
- Compute data for extrapolation to the limit to predict the true natural
  frequencies.

## Definition of the basic inputs

The finite element code realize on the basic functionality implemented in this
package.

````julia
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, solve_blocked!
using FinEtoolsDeforLinear
using SymRCM
using Arpack
using LinearAlgebra
````

The material parameters may be defined with the specification of the units.
The elastic properties are:

````julia
E = 200.0 * phun("GPa")
nu = 0.3;
````

The mass density is

````julia
rho = 8000 * phun("kg/m^3")
````

The geometry consists of the radius of the circle, and the diameter of the
circular cross-section.

````julia
radius = 1.0 * phun("m"); diameter = 0.1 * phun("m");
````

We shall calculate these eigenvalues, but we are mostly interested in the
first three  natural frequencies.

````julia
neigvs = 18;
````

The mass shift needs to be applied since the structure is free-floating.

````julia
oshift = (2*pi*15)^2
````

We will generate this many elements per radius of the cross-section, and along
the length of the circular ring.

````julia
nrs, nLs = 2 .^(1:3), 20*2 .^(1:3)
results = let
    results = []
    for (nperradius, nL) in zip(nrs, nLs)
        tolerance = diameter/nperradius/100
````

Generate the mesh in a straight cylinder.

````julia
        fens, fes = H8cylindern(diameter/2, 2*pi, nperradius, nL)
````

Select the nodes in the bases of the cylinder for future merging of the nodes.

````julia
        z0l = selectnode(fens, plane = [0, 0, 1, 0], inflate = tolerance)
        z2l = selectnode(fens, plane = [0, 0, 1, 2*pi], inflate = tolerance)
````

Twist the straight cylinder into a ring.

````julia
        for i in 1:count(fens)
            a = fens.xyz[i, 3]
            x = fens.xyz[i, 1]
            y = fens.xyz[i, 2] + radius
            fens.xyz[i, :] .= (x, y*cos(a), y*sin(a))
        end
````

Merge the nodes of the bases, which involves renumbering the connectivity.

````julia
        fens, fes = mergenodes(fens, fes, tolerance, vcat(z0l, z2l))
````

If desired, visualize the mesh.
File = "ring.vtk"
vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" $File`)

Generate the material and the FEM machine.

````julia
        MR = DeforModelRed3D
        material = MatDeforElastIso(MR, rho, E, nu, 0.0)
        femm = FEMMDeforLinearESNICEH8(MR, IntegDomain(fes, NodalTensorProductRule(3)), material)
````

Set up the nodal fields for the geometry in the displacements.

````julia
        geom = NodalField(fens.xyz)
        u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
````

Renumber the nodes to produce quicker linear equation solves.

````julia
        p = let
            C = connectionmatrix(femm, count(fens))
            p = symrcm(C)
        end
        numberdofs!(u, p)
````

Now set up the discrete model.

````julia
        associategeometry!(femm,  geom)
        K  = stiffness(femm, geom, u)
        M = mass(femm, geom, u)
````

The matrix partitioning now must be enforced to accommodate the prescribed
displacements and rotations.

````julia
        K_ff = matrix_blocked(K, nfreedofs(u))[:ff]
        M_ff = matrix_blocked(M, nfreedofs(u))[:ff]
````

Solve the free vibration problem.

````julia
        evals, evecs, nconv = eigs(Symmetric(K_ff+oshift*M_ff), Symmetric(M_ff); nev=neigvs, which=:SM, explicittransform = :none)
````

Correct for the mass shift.

````julia
        evals = evals .- oshift;
        sigdig(n) = round(n * 10000) / 10000
        fs = real(sqrt.(complex(evals)))/(2*pi)
        println("Eigenvalues: $(sigdig.(fs)) [Hz]")
        push!(results, (fs[6+1], fs[6+3], fs[6+5]))
    end
    results # return these values from the block
end

@show results
````

## Richardson extrapolation

Here we will use Richardson extrapolation from the three sets of data. This
will allow us to predict the convergence rate and the true solution for each
of the three frequencies (or rather the pairs of frequencies, 7 and 8, 9 and
10, and 11 and 12).

We will immediately set up the convergence plots. We will extrapolate and then
compute from that the normalized error to be plotted with respect to the
refinement factors (which in this case are 4, 2, and 1). We use the
refinement factor as a convenience: we will calculate the element size by
dividing the circumference of the ring with a number of elements generated
circumferentially.

Present the convergence data in graphical form.

````julia
using FinEtools.AlgoBaseModule: richextrapol
using PlotlyJS

x = 2*pi*radius./nLs ./ radius;
````

Modes 7 and 8

````julia
sols = [r[1] for r in results]
resextrap = richextrapol(sols, [4.0, 2.0, 1.0])
print("Predicted frequency 7 and 8: $(resextrap[1])\n")
y = abs.(sols .- resextrap[1])./resextrap[1]
tc78 = scatter(; x=x, y=y, mode="markers+lines", name="Mode 7, 8", line_color="rgb(215, 15, 15)",
    marker=attr(size=9, symbol="diamond-open"))
````

Modes 9 and 10

````julia
sols = [r[2] for r in results]
resextrap = richextrapol(sols, [4.0, 2.0, 1.0])
print("Predicted frequency 9 and 10: $(resextrap[1])\n")
y = abs.(sols .- resextrap[1])./resextrap[1]
tc910 = scatter(; x=x, y=y, mode="markers+lines", name="Mode 9, 10", line_color="rgb(15, 215, 15)",
    marker=attr(size=9, symbol="square-open"))
````

Modes 11 and 12

````julia
sols = [r[3] for r in results]
resextrap = richextrapol(sols, [4.0, 2.0, 1.0])
print("Predicted frequency 11 and 12: $(resextrap[1])\n")
y = abs.(sols .- resextrap[1])./resextrap[1]
tc1112 = scatter(; x=x, y=y, mode="markers+lines", name="Mode 11, 12", line_color="rgb(15, 15, 215)",
    marker=attr(size=9, symbol="circle-open"))
````

Set up the layout:

````julia
layout = Layout(;
    xaxis=attr(title="Normalized Element Size [ND]", type="log"),
    yaxis=attr(title="Normalized error [ND]", type="log"),
    title="3D: Convergence of modes 7, ..., 12")
````

Plot the graphs:

````julia
config  = PlotConfig(plotlyServerURL="https://chart-studio.plotly.com", showLink=true)
pl = plot([tc78, tc910, tc1112], layout; config = config)
display(pl)

nothing
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

