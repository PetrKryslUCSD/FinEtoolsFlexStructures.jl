# # Modal analysis of a free-floating steel circle

# Source code: [`circle_modal_tut.jl`](circle_modal_tut.jl)

# ## Description

# Vibration analysis of a beam simply supported in one plane, and clamped in
# another. The results are compared with analytical expressions. This is a
# benchmark from the NAFEMS Selected Benchmarks for Natural Frequency Analysis,
# publication: Test VM09: Circular Ring --  In-plane and Out-of-plane
# Vibration.

# ## Reference frequencies 

# There will be 6 rigid body modes (zero natural frequencies).

# The numerical results are due to the publication: 
# NAFEMS Finite Element Methods & Standards, Abbassian, F., Dawswell, D. J., and
# Knowles, N. C. Selected Benchmarks for Natural Frequency Analysis, Test No.
# 6. Glasgow: NAFEMS, Nov., 1987. 

# The reference values were analytically determined (Blevins, FORMULAS FOR
# DYNAMICS, ACOUSTICS AND VIBRATION, Table 4.16). Note that shear flexibility
# was neglected when computing the reference values.

# | Mode       |         Reference Value (Hz)  |  NAFEMS Target Value (Hz) | 
# | -------   |     -------  |  ------- | 
# | 7, 8 | (out of plane)   |        51.85          |         52.29  | 
# | 9, 10 |  (in plane)       |       53.38         |          53.97  | 
# | 11, 12 |  (out of plane)   |     148.8          |         149.7  | 
# | 13, 14 |  (in plane)       |     151.0          |         152.4  | 
# | 15, 16 |  (out of plane)   |     287.0          |         288.3  | 
# |  17, 18 |  (in plane)      |      289.5         |          288.3  | 

# ## Goals

# - Show convergence relative to reference values. 
# - Demonstrate the optimization of eigenvalue accuracy by choosing mass type.

# ## Definition of the basic inputs

# The finite element code realize on the basic functionality implemented in this
# package.
using FinEtools

using LinearAlgebra
# The material parameters may be defined with the specification of the units.
# The elastic properties are:
E = 200.0 * phun("GPa") 
nu = 0.3;

# The mass density is
rho = 8000 * phun("kg/m^3")
# Here are the cross-sectional dimensions and the length of the beam between
# supports.
radius = 1.0 * phun("m"); diameter = 0.1 * phun("m"); 

# ## Cross-section

# Cross-sectional properties are incorporated in the cross-section property. The
# three arguments supplied are functions. All are returning "constants". In
# particular the first two functions each return the dimension of the
# cross-section as a constant(the beam has a uniform cross-section); the third
# function defines the orientation of the cross-section in the global Cartesian
# coordinates. `[1.0, 0.0, 0.0]` is the vector that together with the tangent
# to the midline curve of the beam spans the $x_1x_2$ plane of the local
# coordinates for the beam.
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionCircle
cs = CrossSectionCircle(s -> diameter/2, s -> [1.0, 0.0, 0.0])
@show cs.parameters(0.0)

# For instance the the first out of plane mode is listed in the reference cited
# above as
R = radius
I = cs.parameters(0.0)[4]
m = rho * cs.parameters(0.0)[1]
J = cs.parameters(0.0)[2]
G = E/2/(1+nu)
i = 2 # the first non-rigid body mode
@show i*(i^2-1)/(2*pi*R^2)*sqrt(E*I/m/(i^2+E*I/G/J))

# The first "ovaling" (in-plane) mode is:
i = 2 # the first ovaling mode
@show i*(i^2-1)/(2*pi*R^2*(i^2+1)^(1/2))*sqrt(E*I/m)



# The purpose of the numerical model is to calculate approximation to the reference frequencies.

neigvs = 18;

# We will generate
n = 20
# beam elements along the member.
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member
tolerance = radius/n/1000;
fens, fes = frame_member([0 0 0; 2*pi 0 0], n, cs)
for i in 1:count(fens)
    a = fens.xyz[i, 1]
    fens.xyz[i, :] .= (radius+radius*cos(a), radius*sin(a), 0)
end
fens, fes = mergenodes(fens, fes, tolerance, [1, n+1])

# ## Material

# Material properties can be now used to create a material: isotropic elasticity model of the `FinEtoolsDeforLinear` package is instantiated.
using FinEtoolsDeforLinear
material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

# ## Fields

# Now we start constructing the discrete finite element model.
# We begin by constructing the requisite fields, geometry and displacement.
# These are the so-called "configuration variables", all initialized to 0.
# This is that geometry field.
geom0 = NodalField(fens.xyz)
# This is the displacement field, three unknown displacements per node.
u0 = NodalField(zeros(size(fens.xyz, 1), 3))
# This is the rotation field, three unknown rotations per node are represented
# with a rotation matrix, in total nine numbers. The utility function
# `initial_Rfield`
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield
Rfield0 = initial_Rfield(fens)
# Here we verify the number of nodes and the number of degrees of freedom in the
# rotation field per node.
@show nents(Rfield0)
@show ndofs(Rfield0)

# Finally, this is the displacement and rotation field for incremental changes,
# incremental displacements and incremental rotations. In total, 6 unknowns per
# node.
dchi = NodalField(zeros(size(fens.xyz, 1), 6))

# There are no support conditions.
applyebc!(dchi)
# The  the number of free
# (unknown) degrees of freedom is equal to the total number of degrees of freedom in the system.
numberdofs!(dchi);


# ## Assemble the global discrete system

using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material);

# For disambiguation we will refer to the stiffness and mass functions by qualifying them with the corotational-beam module, `FEMMCorotBeamModule`.
using FinEtoolsFlexStructures.FEMMCorotBeamModule
CB = FEMMCorotBeamModule
# Thus we can construct the stiffness and mass matrix as follows:
# Note that the finite element machine is the first argument. This provides
# access to the integration domain. The next argument is the geometry field,
# followed by the displacement, rotations, and incremental
# displacement/rotation fields. 
K = CB.stiffness(femm, geom0, u0, Rfield0, dchi);
M = CB.mass(femm, geom0, u0, Rfield0, dchi);
# We can compare the size of the stiffness matrix with the number of degrees of
# freedom that are unknown (20).
@show size(K)

# ## Solve the free-vibration problem

oshift = (2*pi*15)^2

# The Arnoldi algorithm implemented in the well-known `Arpack` package is used
# to solve the generalized eigenvalue problem with the sparse matrices. As is
# common in structural dynamics, we request the smallest eigenvalues in
# absolute value (`:SM`). 
using Arpack
evals, evecs, nconv = eigs(Symmetric(K + oshift * M), Symmetric(M); nev=neigvs, which=:SM, explicittransform = :none);
# First  we should check that the requested eigenvalues actually converged:
@show nconv == neigvs

# Make sure the eigenvalues and eigenvectors are stripped of the imaginary part.
evals = real.(evals)
evecs = real.(evecs)

# The eigenvalues (i. e. the squares of the angular frequencies) are returned in
# the vector `evals`. The mode shapes constitute the columns of the matrix `evecs`.
@show size(evecs)
# The natural frequencies are obtained from the squares of the angular
# frequencies. We note the use of `sqrt.` which broadcast the square root over
# the array `evals`.
fs = sqrt.([max(0, e - oshift) for e in evals]) / (2 * pi);

# ## Comparison of computed and analytical results

# The approximate and analytical frequencies are now reported.
sigdig(n) = round(n * 10000) / 10000
println("Approximate frequencies: $(sigdig.(fs)) [Hz]")

# ## Set up the visualization of the vibration modes

# The animation will show one of the vibration modes overlaid on the undeformed geometry. The configuration during the animation needs to reflect rotations. The function `update_rotation_field!` will update the rotation field given a vibration mode.
using FinEtoolsFlexStructures.RotUtilModule: update_rotation_field!

# The visualization utilities take advantage of the PlotlyJS library.
using PlotlyJS
using VisualStructures: plot_space_box, plot_solid, render, react!, default_layout_3d, save_to_json

# The magnitude of the vibration modes (displacements  and rotations) will be amplified with this scale factor:
scale = 1.5

# In order to handle variables inside loops correctly, we create a local scope with the `let end` block.
vis(mode) = let
    # The extents of the box will be preserved during animation in order to eliminate changes in the viewing parameters.
    tbox = plot_space_box(reshape(inflatebox!(boundingbox(fens.xyz), 0.3 * radius), 2, 3))
    # This is the geometry of the structure without deformation (undeformed). It is displayed as gray, partially transparent.
    tenv0 = plot_solid(fens, fes; x=geom0.values, u=0.0 .* dchi.values[:, 1:3], R=Rfield0.values, facecolor="rgb(125, 155, 125)", opacity=0.3);
    # Initially the plot consists of the box and the undeformed geometry.
    plots = cat(tbox, tenv0; dims=1)
    # Create the layout for the plot. Set the size of the window.
    layout = default_layout_3d(;width=600, height=600)
    # Set the aspect mode to get the correct proportions.
    layout[:scene][:aspectmode] = "data"
    # Render the undeformed structure
    pl = render(plots; layout=layout, title="Mode $(mode)")
    sleep(2.115)
    # This is the animation loop. 
    # 1. Distribute a fraction of the selected eigenvector into the incremental displacement/rotation field.
    # 2. Create the deformed configuration by defining displacement field `u1` and rotation field `Rfield1`.
    # 3. Create the plot for the deformed configuration, and add it to the list of plots.
    # 4. Call the `react!` function to update the display. Sleep for a brief period of time to give the display a chance to become current.
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
end

# ## Visualize vibration mode

# Animate the harmonic motion of the mode given as argument:
# vis(7)

using FinEtoolsFlexStructures.FEMMCorotBeamModule: MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA,
MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA, 
MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA, 
MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA

results = let
    results = Dict()
    for mtype in [
        MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA, 
        MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA, 
        MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA, 
        MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA]
        M = CB.mass(femm, geom0, u0, Rfield0, dchi; mass_type = mtype);
        
        evals, evecs, nconv = eigs(Symmetric(K + oshift * M), Symmetric(M); nev=neigvs, which=:SM, explicittransform = :none);
        evals = real.(evals)
        evecs = real.(evecs)
        results[mtype] = evals, evecs
    end
    results
end

colors = [
"rgb(125, 15, 15)", 
"rgb(15, 155, 15)", 
"rgb(15, 15, 155)", 
"rgb(125, 115, 115)"
]

evals = results[MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA][1]
x = 1:length(evals); y = sqrt.([max(0, e - oshift) for e in evals]) / (2 * pi);
tc0 = scatter(;x=x, y=y, mode="markers", name = "cons, wo", line_color = "rgb(215, 15, 15)", marker = attr(size = 9, symbol = "diamond-open"))
evals = results[MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA][1]
x = 1:length(evals); y = sqrt.([max(0, e - oshift) for e in evals]) / (2 * pi);
tc1 = scatter(;x=x, y=y, mode="markers", name = "cons, w", line_color = "rgb(15, 215, 15)", marker = attr(size = 9, symbol = "triangle-down"))
evals = results[MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA][1]
x = 1:length(evals); y = sqrt.([max(0, e - oshift) for e in evals]) / (2 * pi);
tc2 = scatter(;x=x, y=y, mode="markers", name = "lumped, wo", line_color = "rgb(15, 15, 215)", marker = attr(size = 9, symbol = "x-open"))
evals = results[MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA][1]
x = 1:length(evals); y = sqrt.([max(0, e - oshift) for e in evals]) / (2 * pi);
tc3 = scatter(;x=x, y=y, mode="markers", name = "lumped, w", line_color = "rgb(165, 165, 15)", marker = attr(size = 9, symbol = "square"))

# 7, 8 (out of plane)         51.85                 52.29 
# 9, 10 (in plane)            53.38                 53.97 
# 11, 12 (out of plane)      148.8                 149.7 
# 13, 14 (in plane)          151.0                 152.4 
# 15, 16 (out of plane)      287.0                 288.3 
# 17, 18 (in plane)          289.5                 288.3 
rfs = vec(Float64[0 0 0 0 0 0 51.85 51.85 53.38 53.38 148.8 148.8 151.0 151.0 287.0 287.0 289.5 289.5])

rtc = scatter(;x=collect(1:length(rfs)), y=rfs, mode="lines", name = "ref", line_color = "rgb(15, 15, 15)")

# Set up the layout:
layout = Layout(;xaxis=attr(title="Mode", type = "linear"), yaxis=attr(title="Frequency [hertz]", type = "linear"), title = "Comparison of mass types")
# Plot the graphs:
config  = PlotConfig(plotlyServerURL="https://chart-studio.plotly.com", showLink=true)
pl = plot([rtc, tc0, tc1, tc2, tc3], layout; config = config)
display(pl)


# ## "Mixed" mass matrix

# As the graph suggests we can try mixing together mass matrices computed from
# different assumptions (lumped versus consistent, with or without rotation
# inertia). No particular justification can be provided for these numbers other
# than best low-frequency matching for linear bars suggests that 50% can be
# taken lumped and 50% consistent. Here it is obviously more complicated due to
# the presence of rotation inertia.
M = 0.5 .* CB.mass(femm, geom0, u0, Rfield0, dchi; mass_type = MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA) + 
    0.5 .* CB.mass(femm, geom0, u0, Rfield0, dchi; mass_type = MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA);
# With this mixed mass matrix we solve the free vibration problem again.
evals, evecs, nconv = eigs(Symmetric(K + oshift * M), Symmetric(M); nev=neigvs, which=:SM, explicittransform = :none);

# Plotting the newly obtained data on top of the previously presented data, we
# can observe sometimes substantial improvement of accuracy of the mixed-matrix
# formulation relative to the individual mass matrix types.
x = 1:length(evals); y = sqrt.([max(0, e - oshift) for e in evals]) / (2 * pi);
mtc = scatter(;x=x, y=y, mode="markers", name = "mixed", line_color = "rgb(215, 15, 215)", marker = attr(size = 9, symbol = "circle"))
config  = PlotConfig(plotlyServerURL="https://chart-studio.plotly.com", showLink=true)
pl = plot([rtc, tc0, tc1, tc2, tc3, mtc], layout; config = config)
display(pl)

nothing

