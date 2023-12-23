# # Modal analysis of simply supported/clamped beam

# Source code: [`beam_modal_tut.jl`](beam_modal_tut.jl)

# Last updated: 12/23/23

# ## Description

# Vibration analysis of a beam simply supported in one plane, and clamped
# in another. The results are compared with analytical expressions.

# ## Goals

# - Introduce definition of cross-section.
# - Show generation of finite element mesh of beams.
# - Describe geometry, displacement, and rotation fields.
# - Describe application of support conditions.
# - Calculate the discrete model quantities and solve the free vibration problem.
# - Demonstrate visualization of the free vibrations.

# ## Definition of the basic inputs

# We will probably need some linear algebra functions.
using LinearAlgebra

# The finite element code relies on the basic functionality implemented in this
# package.
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked

# The material parameters may be defined with the specification of the units.
# The elastic properties are:
E = 30002.0 * phun("ksi") 
nu = 0.0;

# The mass density is expressed in customary units as
rho = 0.28 * phun("lbm/in^3")

# Here are the cross-sectional dimensions and the length of the beam between supports.
b = 1.8 * phun("in"); h = 1.8 * phun("in"); L = 100 * phun("in");

# ## Analytical frequencies 

# The analytical frequencies were taken from table 8-1 of Formulas for natural
# frequency and mode shape, Robert D. Blevins, Krieger publishing company,
# Malabar Florida, reprint edition 2001. Available also in FORMULAS FOR
# DYNAMICS, ACOUSTICS AND VIBRATION, by the same author, Wiley, 2016, Table 4.2.
# 
# The beam is aligned with the $Y$ global Cartesian coordinate. The beam
# behaves as a simply supported beam in the vertical plane (global Cartesian
# $YZ$), while in the horizontal plane (global Cartesian $XY$) it behaves
# as a clamped beam.
# 
# The cross-sectional properties are:
A = b * h;
I2 = b * h^3 / 12;
I3 = b^3 * h / 12;
# Then the analytical vibration frequencies for the first two modes are:
@show analyt_freq = [(1 * pi)^2, (4.73004074)^2] .* (sqrt(E * I2 / rho / A) / (2 * pi * L^2));

# The purpose of the numerical model is to calculate approximation to these two
# analytical natural frequencies.
neigvs = length(analyt_freq);
 

# ## Cross-section definition

# Cross-sectional properties are incorporated in the cross-section object. The
# three arguments supplied are functions. All are returning "constants", as
# appropriate for a uniform cross section beam. In particular the first two
# functions each return the dimension of the cross-section as a constant; the
# third function defines the orientation of the cross-section in the global
# Cartesian coordinates. `[1.0, 0.0, 0.0]` is the vector that together with the
# tangent to the midline curve of the beam spans the $x_1-x_2$ plane of the
# local coordinates for the beam.
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
cs = CrossSectionRectangle(s -> b, s -> h, s -> [1.0, 0.0, 0.0])

# We can compare the analytical values of the cross-section properties with
# those stored in the `cs` structure:
@show A, I2, I3
@show cs.parameters(0.0)

# ## Mesh generation

# Now we generate the mesh of the beam. The locations of its two endpoints are:
xyz = [[0 -L / 2 0]; [0 L / 2 0]]
# We will generate
n = 4
# beam elements along the member. The frame member is subdivided into a given
# number of finite elements, which are given the cross section properties
# defined in `cs`.
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member
fens, fes = frame_member(xyz, n, cs);
# The mesh definition consists of the nodes
@show fens
# and the finite elements
@show fes
# Note  that the cross-sectional properties are incorporated through `cs`.

# ## Material

# Material properties can be now used to create a material: isotropic elasticity
# model of the `FinEtoolsDeforLinear` package is instantiated.
using FinEtoolsDeforLinear
material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
# Since we are interested in dynamics, the mass density must be supplied.

# ## Fields

# Now we start constructing the discrete finite element model.
# We begin by constructing the requisite fields, geometry and displacement.
# These are the so-called "configuration variables", all initialized to zero.

# This is the geometry field. It describes the locations of the nodes.
geom0 = NodalField(fens.xyz)
# This is the displacement field, three unknown displacements per node.
u0 = NodalField(zeros(size(fens.xyz, 1), 3))
# This is the rotation (orientation) field, where the three unknown rotations
# per node are represented with a rotation matrix, in total nine numbers. The
# utility function `initial_Rfield` is used to construct this initial
# orientation field where each orientation matrix is the identity.
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield
Rfield0 = initial_Rfield(fens)
# Here we verify the number of nodes and the number of degrees of freedom in the
# rotation field per node.
@show nents(Rfield0)
@show ndofs(Rfield0)

# Finally, this is the displacement and rotation field for incremental changes,
# incremental displacements and incremental rotations. In total, 6 unknowns per
# node. All degrees of freedom are vector components (displacement vector and
# rotation vector).
dchi = NodalField(zeros(size(fens.xyz, 1), 6))


# ## Support conditions

# Now we apply the essential boundary conditions (EBCs) to enforce the action of
# the supports at the ends of the beam. 

# First we select the node at the location  `[0 -L/2 0]`. We do that by
# searching for all nodes that are located inside a tiny box centered at this
# location. The tolerance of a fraction of the length of an element (i. e. the
# distance between the nodes) is used to take the location and blow it up into
# a nonzero-volume box.
l1 = selectnode(fens; box=[0 0 -L / 2 -L / 2 0 0], tolerance=L / n / 1000)
# The nodes in this list should consist of a single node (which is why the
# tolerance is so small, in order to limit the selection to a single node near
# the location given).
@show length(l1) == 1
# The boundary condition at this point dictates zero displacements (degrees of
# freedom 1, 2, and 3) and zero rotations about $Y$ (5) and $Z$ (6). This
# leaves rotation about the $X$ free.
for i in [1,2,3,5,6]
    setebc!(dchi, l1, true, i)
end
# Similarly, the node next to the other end of the beam is selected.
l1 = selectnode(fens; box=[0 0 L / 2  L / 2 0 0], tolerance=L / n / 1000)
# And an identical boundary condition combination is enforced.
for i in [1,2,3,5,6]
    setebc!(dchi, l1, true, i)
end

# The essential boundary conditions will  reduce the number of free
# (unknown) degrees of freedom. All degrees of freedom are assigned numbers.
numberdofs!(dchi);
# The number of three degrees of freedom is
@show nfreedofs(dchi)
# Here we inspect the degrees of freedom in the incremental
# displacement/rotation field:
@show dchi.dofnums
# Note that the degrees of freedom are actually carried by the incremental
# field, not by the displacement or the rotation fields. There are therefore 6
# degrees of freedom per node in the incremental displacement field.



# ## Assemble the global discrete system

# The integration domain of the beam is one-dimensional. The quantities in the
# co-rotational beam are really all evaluated analytically, so the numerical
# integration does have no effect.
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material);

# For disambiguation we will refer to the `stiffness` and `mass` functions by
# qualifying them with the corotational-beam module, `FEMMCorotBeamModule`.
# We will use the abbreviation `CB`. 
using FinEtoolsFlexStructures.FEMMCorotBeamModule
CB = FEMMCorotBeamModule

# Thus we can construct the stiffness and mass matrix as follows: Note that the
# finite element machine (FEMM) is the first argument. This provides access to
# the integration domain. The next argument is the geometry field, followed by
# the displacement, rotations, and incremental displacement/rotation fields.
K = CB.stiffness(femm, geom0, u0, Rfield0, dchi);
M = CB.mass(femm, geom0, u0, Rfield0, dchi);
# We can compare the size of the stiffness matrix with the number of degrees of
# freedom that are unknown (`nfreedofs(dchi)`).
@show size(K)
# We can see that the system matrices have a row and column for each degree of
# freedom in the system, whether known or unknown.

# The matrix partitioning now must be enforced to reflect the prescribed
# displacements and rotations.

K_ff = matrix_blocked(K, nfreedofs(dchi))[:ff]
M_ff = matrix_blocked(M, nfreedofs(dchi))[:ff]

# ## Solve the free-vibration problem

# The Arnoldi algorithm implemented in the well-known `Arpack` package is used
# to solve the generalized eigenvalue problem with the sparse matrices. As is
# common in structural dynamics, we request the smallest eigenvalues in
# absolute value (`:SM`). In order to help the function decide which algorithm
# is suitable, we explicitly designate the matrices as symmetric.
using Arpack
evals, evecs, nconv = eigs(Symmetric(K_ff), Symmetric(M_ff); nev=neigvs, which=:SM, explicittransform = :none);

# We should check that the requested eigenvalues actually converged:
@show nconv == neigvs

# The eigenvalues (i. e. the squares of the angular frequencies) are returned in
# the vector `evals`. The mode shapes constitute the columns of the matrix
# `evecs`.
@show size(evecs)

# The natural frequencies are obtained from the squares of the angular
# frequencies. We note the use of `sqrt.` which broadcasts the square root over
# the array `evals`.
fs = sqrt.(evals) / (2 * pi);


# ## Comparison of computed and analytical results

# The approximate and analytical frequencies are now reported.
println("Approximate frequencies: $fs [Hz]")
println("Analytical frequencies: $analyt_freq [Hz]")

  
# Close agreement between the approximate and analytical frequencies can be
# observed: The error of the numerical solution is a fraction of a percent.
errs = abs.(analyt_freq .- fs) ./ analyt_freq
println("Relative errors of frequencies: $errs [ND]")


# ## Visualize vibration modes

# The animation will show one of the vibration modes overlaid on the undeformed
# geometry.

# The visualization utilities take advantage of the `PlotlyJS` library through
# the convenience package `VisualStructures`.

using VisualStructures: plot_space_box, plot_solid, render, react!, default_layout_3d, save_to_json

# The configuration during the animation needs to reflect rotations.
# The function `update_rotation_field!` will update the rotation field given a
# vibration mode.
using FinEtoolsFlexStructures.RotUtilModule: update_rotation_field!

# The magnitude of the vibration modes (displacements  and rotations) will be
# amplified with this scale factor:
scale = 1.2

# This is the mode that will be animated:
mode = 1

# In order to handle variables inside loops correctly, we create a local scope
# with the `let end` block.
let
    # The extents of the box will be preserved during animation in order to
    # eliminate changes in the viewing parameters.
    tbox = plot_space_box([[-0.2 * L -L / 2 -0.2 * L]; [+0.2 * L L / 2 +0.2 * L]])
    # This is the geometry of the structure without deformation (undeformed). It
    # is displayed as gray, partially transparent.
    tenv0 = plot_solid(fens, fes; x=geom0.values, u=0.0 .* dchi.values[:, 1:3], R=Rfield0.values, facecolor="rgb(125, 155, 125)", opacity=0.3);
    # Initially the plot consists of the box and the undeformed geometry.
    plots = cat(tbox, tenv0; dims=1)
    # Create the layout for the plot. Set the size of the window.
    layout = default_layout_3d(; options = Dict(:responsive=>true))
    # Set the aspect mode to get the correct proportions.
    layout[:scene][:aspectmode] = "data"
    # Render the undeformed structure
    pl = render(plots; layout=layout, title="Mode $(mode)")
    sleep(2.115)
    # This is the animation loop. 
    # - Distribute a fraction of the selected eigenvector into the incremental displacement/rotation field.
    # - Create the deformed configuration by defining displacement field `u1` and rotation field `Rfield1`. 
    # - Create the plot for the deformed configuration, and add it to the list of plots.
    # - Call the `react!` function to update the display. Sleep for a brief period of time to give the display a chance to become current.
    for xscale in scale .* sin.(collect(0:1:89) .* (2 * pi / 21))
        scattersysvec!(dchi, xscale .* evecs[:, mode])
        u1 = deepcopy(u0)
        u1.values .= dchi.values[:, 1:3]
        Rfield1 = deepcopy(Rfield0)
        update_rotation_field!(Rfield1, dchi)
        tenv1 = plot_solid(fens, fes; x=geom0.values, u=dchi.values[:, 1:3], R=Rfield1.values, facecolor="rgb(50, 55, 125)");
        plots = cat(tbox, tenv0, tenv1; dims=1)
        react!(pl, plots, pl.plot.layout)
        sleep(0.115)
    end
    # Save the plot to a Json file. It can be then re-displayed later.
    save_to_json(pl, "deformed_plot.json")
end

# Load the plot from a file.
using VisualStructures: plot_from_json
plot_from_json("deformed_plot.json")

# Nothing to return.
nothing
