# # Prestressed simply-supported column modal analysis

# Source code: [`prestressed_column_modal_tut.jl`](prestressed_column_modal_tut.jl)

# ## Description

# Vibration analysis of a simply supported column loaded with axial force. 
# The fundamental vibration frequency depends on the prestress force.

# ## Goals

# - Demonstrate the change in the fundamental vibration frequency due to the
#   prestress through an axial force.
# - Show how to evaluate the geometric stiffness from a static solution.
# - Solve the eigenvalue free vibration problem with the inclusion of the
#   geometric stiffness matrix.

# ## Definition of the basic inputs

# The finite element code realize on the basic functionality implemented in this
# package.
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked

# The material parameters may be defined with the specification of the units.
# The elastic properties are:
E = 30002.0 * phun("ksi") 
nu = 0.0;

# The mass density is expressed in customary units as
g = 32.17*12 * phun("in/sec^2")
rho = 4.65 * phun("oz/in ^3") / g
# Here are the cross-sectional dimensions and the length of the beam between supports.
b = 1.8 * phun("in"); h = 1.8 * phun("in"); L = 300 * phun("in");

# ## Cross-section

# Cross-sectional properties are incorporated in the cross-section property. The
# three arguments supplied are functions. All are returning "constants". In
# particular the first two functions each return the dimension of the
# cross-section as a constant(the beam has a uniform cross-section); the third
# function defines the orientation of the cross-section in the global Cartesian
# coordinates. `[1.0, 0.0, 0.0]` is the vector that together with the tangent
# to the midline curve of the beam spans the $x_1x_2$ plane of the local
# coordinates for the beam.
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
cs = CrossSectionRectangle(s -> b, s -> h, s -> [1.0, 0.0, 0.0])

# Here we retrieve the cross-sectional properties at the arc length 0.0.
@show A, J, I1, I2, I3 = cs.parameters(0.0)

# ## Analytical frequencies 

# The analytical frequencies were taken from table 8-1 of Formulas for
# natural frequency and mode shape, Robert D. Blevins, Krieger publishing
# company, Malabar Florida, reprint edition 2001.
# 
# The beam has cylindrical supports at either end. 
# 
# Simply supported column without a pre-stressing force has a fundamental frequency of
@show analyt_freq = (1*pi)^2/(2*pi*L^2)*sqrt(E*I2/rho/A);


# The critical Euler buckling load (simple-support): (pi^2*E*I2/L^2).
@show PEul = (pi^2*E*I2/L^2);

# Ps=linspace(-1, 1, 100);

# The purpose of the numerical model is to calculate approximation to the fundamental
# analytical natural frequency.
neigvs = 1;

# Now we generate the mesh of the beam. The locations of its two endpoints are:
xyz = [[0 -L/2 0]; [0 L/2 0]]
# We will generate
n = 20
# beam elements along the member.
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member
fens, fes = frame_member(xyz, n, cs);


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

# Finally, this is the displacement and rotation field for incremental changes,
# incremental displacements and incremental rotations. In total, 6 unknowns per
# node.
dchi = NodalField(zeros(size(fens.xyz, 1), 6))

# ## Support conditions

# Now we apply the essential boundary conditions (EBCs) to enforce the action of
# the supports at the ends of the beam. 

# First we select the node at the location  `[0 -L/2 0]`. This is the immovable node.
immovable = selectnode(fens; box=[0 0 -L/2 -L/2 0 0], tolerance=L/n/1000)[1]

# The boundary condition at this point dictates zero displacements (degrees of
# freedom 1, 2, and 3) and zero rotations about the axis of the beam (degree of freedom 5). 
for i in [1,2,3,5]
    setebc!(dchi, [immovable], true, i)
end

# Similarly, the node next to the other end of the beam is selected. This node
# is to be exposed to the external loading through the axial force.
movable = selectnode(fens; box=[0 0 L/2  L/2 0 0], tolerance=L/n/1000)[1]

# This time the transverse displacements and the axial rotation are suppressed.
for i in [1,3,5]
    setebc!(dchi, [movable], true, i)
end

# These boundary conditions now need to be "applied". This simply means that the
# prescribed values of the degrees of freedom are copied into the active degrees
# of freedom.
applyebc!(dchi)
# The essential boundary conditions will also reduce the number of free
# (unknown) degrees of freedom.
numberdofs!(dchi);
# Here we inspect the degrees of freedom in the incremental
# displacement/rotation field:
@show dchi.dofnums
# Note that the degrees of freedom are actually carried by the incremental
# field, not by the displacement or the rotation fields. 


# ## Assemble the global discrete system

using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material);

# For disambiguation we will refer to the stiffness and mass functions by
# qualifying them with the corotational-beam module, `FEMMCorotBeamModule`.
using FinEtoolsFlexStructures.FEMMCorotBeamModule
CB = FEMMCorotBeamModule

# Thus we can construct the stiffness matrix as follows:
# Note that the finite element machine is the first argument. This provides
# access to the integration domain. The next argument is the geometry field,
# followed by the displacement, rotations, and incremental
# displacement/rotation fields. 
K = CB.stiffness(femm, geom0, u0, Rfield0, dchi);

# Extract the free-free block of the matrix.
K_ff = matrix_blocked(K, nfreedofs(dchi))[:ff]

# Now we construct the means of applying the concentrated force of prestress at
# the movable node. The node is the "boundary" of the domain of the column.
loadbdry = FESetP1(reshape([movable], 1, 1))
# The concentrated force can be considered a distributed loading at a single
# point. This distributed loading can be integrated with a quadrature rule
# suitable for a single point domains.
lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))

# We assume that the magnitude of the force is unity. Note that the force is
# applied in the positive direction along the beam, meaning as tensile.
# Therefore, buckling will be caused by applying  a negative buckling factor.
# Numerically, the buckling factor corresponding to failure is then of the
# magnitude of the Euler force. 
fi = ForceIntensity(Float64[0, 1.0, 0, 0, 0, 0]);

# The distributed loading is now integrated over the "volume" of the integration
# domain.
F = CB.distribloads(lfemm, geom0, dchi, fi, 3);
F_f = vector_blocked(F, nfreedofs(dchi))[:f]

# Solve for the displacement under @show the static load
scattersysvec!(dchi, K_ff\F_f);

# Update deflections so that the initial stress can be computed. First the
# displacements:
u1 = deepcopy(u0)
u1.values .= dchi.values[:, 1:3]
# Then the rotations:
Rfield1 = deepcopy(Rfield0)
using FinEtoolsFlexStructures.RotUtilModule:  update_rotation_field!
update_rotation_field!(Rfield1, dchi)


# The static deflection is now used to compute the internal forces
# which in turn lead to the geometric stiffness matrix.
Kg = CB.geostiffness(femm, geom0, u1, Rfield1, dchi);

# Now we can evaluate the mass matrix,
M = CB.mass(femm, geom0, u0, Rfield0, dchi);
# and we have the complete discrete model for the solution of the free vibration problem.


# Extract the free-free block of the matrices.
M_ff = matrix_blocked(M, nfreedofs(dchi))[:ff]

Kg_ff = matrix_blocked(Kg, nfreedofs(dchi))[:ff]


# ## Solve the free-vibration problem

# The Arnoldi algorithm implemented in the well-known `Arpack` package is used
# to solve the generalized eigenvalue problem with the sparse matrices. As is
# common in structural dynamics, we request the smallest eigenvalues in
# absolute value (`:SM`). 
using Arpack


# The prestress-modified natural frequencies are computed in the loop for 
# a number of prestress force values.
Ps = collect(linearspace(-1.0, 1.0, 50)).*PEul
freqs = let freqs =[];
    for P in Ps
        # Note that we take the complete stiffness matrix: elastic plus prestress (initial stress).
        evals, evecs, nconv = eigs(K_ff + P.*Kg_ff, M_ff; nev=neigvs, which=:SM, explicittransform = :none);
        # The fundamental frequency:
        f = sqrt(evals[1]) / (2 * pi);
        push!(freqs, f);
    end
    freqs
end

# Show the normalized force and the fundamental frequencies.
sigdig(n) = round(n * 1000) / 1000
@show sigdig.(Ps./PEul), sigdig.(freqs)

# ## Present a plot

using Gnuplot


@gp  "set terminal windows 0 "  :-

@gp  :- Ps./PEul freqs./analyt_freq " lw 2 lc rgb 'red' with p title 'Fundamental frequency' "  :-

@gp  :- "set xlabel 'P/P_{Euler}'" :-
@gp  :- "set ylabel 'Frequency(P)/Frequency(0) [Hz]'" :-
@gp  :- "set title 'Prestressed column'"

nothing
