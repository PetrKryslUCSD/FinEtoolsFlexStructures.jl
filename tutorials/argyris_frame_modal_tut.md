# Modal analysis of Argyris frame: effect of prestress

Source code: [`argyris_frame_modal_tut.jl`](argyris_frame_modal_tut.jl)

Last updated: 12/23/23

## Description

Vibration analysis of a L-shaped frame under a loading.
The fundamental vibration frequency depends on the prestress force.

## Goals

- Construct an L-shaped frame by merging individual members.
- Compute the geometric stiffness.
- Evaluate the effect of prestress on the fundamental frequency of vibration.

````julia
using LinearAlgebra
````

The finite element code relies on the basic functionality implemented in this
package.

````julia
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, solve_blocked!
````

The linear deformation code will be needed to evaluate the loading.

````julia
using FinEtoolsDeforLinear
````

The functionality for the beam model comes from these modules.

````julia
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.RotUtilModule:  update_rotation_field!
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
using FinEtoolsFlexStructures.FEMMCorotBeamModule
````

The co-rotational beam functionality is in the `FEMMCorotBeamModule` module.
We need to refer to it by the `module.function` notation. `CB` is a handy
abbreviation.

````julia
CB = FEMMCorotBeamModule
````

Parameters: Young's modules and Poisson ratio. Note that we can supply inputs
in particular physical units.

````julia
E = 71240.0 * phun("MPa")
nu = 0.31;
rho = 5000 * phun("kg/m^3");
````

Cross-sectional dimensions and length of each leg in millimeters:

````julia
b = 0.6 * phun("mm"); h = 30.0 * phun("mm"); L = 240.0 * phun("mm");
````

Magnitude of the total applied force.

````julia
magn = 1e-5 * phun("N");
````

Cross-sectional properties are represented by the object as a functions of the
distance along the curve (here those are returning constants). Dimensions of
the cross section and the orientation of the section are defined:

````julia
cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 1.0, 0.0])
````

## Generate the discrete model

Select the number of elements per leg.

````julia
n = 8;
members = Tuple{FENodeSet, AbstractFESet}[]
push!(members, frame_member([0 0 L; L 0 L], n, cs))
push!(members, frame_member([L 0 L; L 0 0], n, cs))
fens, fes = merge_members(members; tolerance = L / 10000)
````

Construct the requisite fields (geometry and displacement).
Initialize configuration variables.

````julia
geom0 = NodalField(fens.xyz)
u0 = NodalField(zeros(size(fens.xyz,1), 3))
Rfield0 = initial_Rfield(fens)
dchi = NodalField(zeros(size(fens.xyz,1), 6))
````

Apply EBC's (supports): one leg is clamped.

````julia
l1 = selectnode(fens; box = [0 0 0 0 L L], tolerance = L / 10000)
for i in [1, 2, 3, 4, 5, 6]
    setebc!(dchi, l1, true, i)
end
applyebc!(dchi)
numberdofs!(dchi);
````

Material properties. This material object represents an isotropic material.

````julia
material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
````

## Solve the static problem to find the internal forces

````julia
@info "Solving the static problem"
````

Assemble the global discrete system. The stiffness and mass matrices are
computed and assembled.

````julia
femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
K = CB.stiffness(femm, geom0, u0, Rfield0, dchi);
M = CB.mass(femm, geom0, u0, Rfield0, dchi);
````

These matrices have dimensions that correspond to the total number of degrees
of freedom, not just the free degrees of freedom (the actual unknowns).

Construct force intensity, select the loaded boundary, and assemble the load
vector.

````julia
tipn = selectnode(fens; box=[L L 0 0  0 0], tolerance=L/n/1000)[1]
loadbdry = FESetP1(reshape([tipn], 1, 1))
lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
fi = ForceIntensity(Float64[-magn, 0, 0, 0, 0, 0]);
F = CB.distribloads(lfemm, geom0, dchi, fi, 3);
````

Solve for the displacement under the static load. This convenience function
solves the system by extracting the partitions based on the number of free
degrees of freedom.

````julia
solve_blocked!(dchi, K, F)
````

Update deflections and rotations so that the initial stress can be computed.
First the displacements:

````julia
u1 = deepcopy(u0)
u1.values .= dchi.values[:, 1:3]
````

Then the rotations:

````julia
Rfield1 = deepcopy(Rfield0)
update_rotation_field!(Rfield1, dchi)
````

The static deflection defined by `u1` and `Rfield1` is now used to compute the
internal forces which in turn lead to the geometric stiffness.

````julia
Kg = CB.geostiffness(femm, geom0, u1, Rfield1, dchi);
````

## Solve the eigenvalue problems for a range of loading factors

````julia
@info "Solving the free vibration problems"
````

Solution of the eigenvalue free-vibration problem is obtained with the Arnoldi method.

````julia
using Arpack
````

We will solve for this many natural frequencies. Then, since they are ordered
by magnitude, we will pick the fundamental by taking the first from the
list.

````julia
neigvs = 4
````

The matrix partitioning now must be enforced to take into account the number
of free degrees of freedom so that we can solve the eigenvalue
(vibration) problem.

````julia
K_ff = matrix_blocked(K, nfreedofs(dchi))[:ff]
M_ff = matrix_blocked(M, nfreedofs(dchi))[:ff]
Kg_ff = matrix_blocked(Kg, nfreedofs(dchi))[:ff]
````

First we will  sweep through the loading factors that are positive, meaning
the force points in the direction in which it was defined (towards the
clamped end of the frame).

````julia
lfp = linearspace(0.0, 68000.0, 400)
fsp = let
    fsp = Float64[]
    for load_factor in lfp
        evals, evecs, nconv = eigs(
            Symmetric(K_ff + load_factor .* Kg_ff),
            Symmetric(M_ff); nev = neigvs, which = :SM, explicittransform = :none)
        e = real(evals[1])
        f = e > 0.0 ? sqrt(e) / (2 * pi) : 0.0
        push!(fsp, f)
    end
    fsp
end
````

Next, we will sweep through a range of negative load factors: this simply
turns the force around so that it points away from the clamped end. This can
also buckle the frame, but the magnitude is higher.

````julia
lfm = linearspace(-109000.0, 0.0, 400)
fsm = let
    fsm = Float64[]
    for load_factor in lfm
        evals, evecs, nconv = eigs(
            Symmetric(K_ff + load_factor .* Kg_ff),
            Symmetric(M_ff); nev = neigvs, which = :SM, explicittransform = :none)
        e = real(evals[1])
        f = e > 0.0 ? sqrt(e) / (2 * pi) : 0.0
        push!(fsm, f)
    end
    fsm
end
````

## Plot of the fundamental frequency as it depends on the loading factor

````julia
using PlotlyJS
````

We concatenate the ranges for the load factors and the calculated fundamental
frequencies and present them in a single plot.

Plot the amplitude of the response curves. We output two curves.
The first for the driving-point FRF:

````julia
x, y = cat(collect(lfp), collect(lfm); dims=1), cat(fsp, fsm; dims=1)
trace1 = scatter(; x=x, y=y, mode="markers", name="Fundamental frequency", line_color="rgb(15, 15, 215)",
    marker=attr(size=4, symbol="diamond-open"))
````

Set up the layout:

````julia
layout = Layout(;
    xaxis=attr(title="Loading factor P", type="linear"),
    yaxis=attr(title="Frequency(P) [Hz]", type="linear"),
    title="Transfer function")
````

Plot the graphs:

````julia
config  = PlotConfig(plotlyServerURL="https://chart-studio.plotly.com", showLink=true)
pl = plot([trace1, ], layout; config = config)
display(pl)
````

Clearly, the curve giving the dependence of the fundamental frequency on the
loading factor consists of two branches. These two branches correspond to two
different buckling modes: one for the positive orientation of the force and
one for the negative orientation.

## Visualize some fundamental mode shapes

````julia
@info "Visualizing the vibration modes"
````

Here we visualize the fundamental vibration modes for different values of the
loading factor.

The package `VisualStructures` specializes in the dynamic visualization of
beam and shell structures.

````julia
using VisualStructures: plot_space_box, plot_solid, render, react!, default_layout_3d, save_to_json
scale = 0.005

vis(loading_factor, evec) = let
    tbox = plot_space_box(reshape(inflatebox!(boundingbox(fens.xyz), 0.5 * L), 2, 3))
    tenv0 = plot_solid(fens, fes; x=geom0.values, u=0.0 .* dchi.values[:, 1:3], R=Rfield0.values, facecolor="rgb(125, 155, 125)", opacity=0.3);
    plots = cat(tbox, tenv0; dims=1)
    layout = default_layout_3d(;width=600, height=600)
    layout[:scene][:aspectmode] = "data"
    pl = render(plots; layout=layout, title = "Loading factor $(loading_factor)")
    sleep(0.5)
    scattersysvec!(dchi, evec)
    scale = L/3 /  max(maximum(abs.(dchi.values[:, 1])), maximum(abs.(dchi.values[:, 2])), maximum(abs.(dchi.values[:, 3])))
    for xscale in scale .* sin.(collect(0:1:89) .* (2 * pi / 21))
        scattersysvec!(dchi, xscale .* evec)
        u1 = deepcopy(u0)
        u1.values .= dchi.values[:, 1:3]
        Rfield1 = deepcopy(Rfield0)
        update_rotation_field!(Rfield1, dchi)
        tenv1 = plot_solid(fens, fes; x=geom0.values, u=dchi.values[:, 1:3], R=Rfield1.values, facecolor="rgb(50, 55, 125)");
        plots = cat(tbox, tenv0, tenv1; dims=1)
        pl.plot.layout[:title] = "Loading factor $(loading_factor)"
        react!(pl, plots, pl.plot.layout)
        sleep(0.115)
    end
end
````

This is the vibration mode in the lead up to the buckling mode for the
positive orientation of the force.

````julia
load_factor = 60000
evals, evecs, nconv = eigs(
    Symmetric(K_ff + load_factor .* Kg_ff),
    Symmetric(M_ff); nev = neigvs, which = :SM, explicittransform = :none)
vis(load_factor, evecs[:, 1])
````

This is the same vibration mode for the negative orientation of the force, but
note that the associated fundamental frequency increased due to the effect of
the force upon the stiffening of the clamped leg of the frame that is now in
tension, and therefore stiffer.

````julia
load_factor = -50000
evals, evecs, nconv = eigs(
    Symmetric(K_ff + load_factor .* Kg_ff),
    Symmetric(M_ff); nev = neigvs, which = :SM, explicittransform = :none)
vis(load_factor, evecs[:, 1])
````

Increasing the load factor in the negative orientation further, the
fundamental frequency will switch: it will be a different mode shape, the one
that is close to the buckling mode shape for this orientation of the force.

````julia
load_factor = -100000
evals, evecs, nconv = eigs(
    Symmetric(K_ff + load_factor .* Kg_ff),
    Symmetric(M_ff); nev = neigvs, which = :SM, explicittransform = :none)
vis(load_factor, evecs[:, 1])

nothing
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

