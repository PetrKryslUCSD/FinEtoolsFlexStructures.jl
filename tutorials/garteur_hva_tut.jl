# # GARTEUR SM-AG19 Testbed: Harmonic Vibration Analysis 

# Source code: [`garteur_hva_tut.jl`](garteur_hva_tut.jl)

# ## Description

# This virtual test application is based on the test article  
# used by the GARTEUR Structures & Materials Action Group 19  
# which organized a Round Robin exercise where 12 European laboratories  
# tested a single structure between 1995 and 1997. The benchmark structure   
# was a laboratory structure built to simulate the dynamic behaviour  
# of an aeroplane. The structure was initially built for a benchmark  
# study on experimental modal analysis conducted by the  
# Structures and Materials Action Group (SM-AG19) of the Group  
# for Aeronautical Research and Technology in EURope (GARTEUR).  
# The test-bed was designed and manufactured by ONERA, France.

# ![](IMAC97photo.png)

# ### References

# [1] Ground Vibration Test Techniques, compiled by A Gravelle, GARTEUR
# Structures & Materials Action Group 19 Technical report TP-115, 1999.
# [2] Etienne Balmes, Jan R. Wright, GARTEUR GROUP ON GROUND VIBRATION TESTING |
# RESULTS FROM THE TEST OF A SINGLE STRUCTURE BY 12 LABORATORIES IN EUROPE,
# Proceedings of DETC'97, 1997 ASME Design Engineering Technical Conferences,
# September 14-17, 1997, Sacramento, California.
# [3] 3M(TM) Viscoelastic Damping Polymer 112 Series,  Technical Data, May 2017.

# ## Goals

# - Show how to construct model from multiple connected beams.
# - Demonstrate the use of massless connectors.
# - Demonstrate the use of point masses.
# - Demonstrate the use of grounded springs.  
# - Illustrate verification of the solution of the free vibration problem. 

# ## Geometry of the testbed airplane.

# It was a rather simple structure which was reasonably dynamically
# representative of a simple airplane structure. It was composed of several beams
# simulating a fuselage with wings and a tail. Wing tip drums allowed to adjust
# bending and torsion frequencies similarly to airplane ones, with some very
# close modal frequencies. 

# ![](garteur-geom.png)

# The script included below defines the geometry of the structure, the
# cross-sectional properties, the connectivity, and the location of the nodes.
include("garteur_geometry_tut.jl")

# The geometry is visualized in the tutorial 
# [garteur_geometry_vis_tut](garteur_geometry_vis_tut.jl).


# ## Material

# Material properties can be now used to create a material: isotropic elasticity
# model of the `FinEtoolsDeforLinear` package is instantiated.
using FinEtoolsDeforLinear
# The material of the structure is aluminum.
# The elastic modulus:
E = 70000.0 * phun("MPa") 
nu = 0.31;
# The mass density:
rho = 2700 * phun("kg/m^3")
alu = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

# Material for the massless connectors has the mass density set to zero;
# otherwise it has the same properties as the aluminum material  of the
# structure.
massless = MatDeforElastIso(DeforModelRed3D, 0.0, alu.E, alu.nu, 0.0)

# This simple function returns material based on the label of the beam elements.
getmaterial(labl) = begin
    if labl >= 7 
        return massless
    end
    return alu
end



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
# node. Note that the the incremental displacements are in general complex.
dchi = NodalField(0.0im .* zeros(size(fens.xyz, 1), 6))

# There are no support conditions.
applyebc!(dchi)
# The  the number of free(unknown) degrees of freedom is equal to the total
# number of degrees of freedom in the system.
numberdofs!(dchi);

# ## Identify support points and locations of sensors

# Suspension points are at these nodes:
suspln = selectnode(fens; box = initbox!(Float64[], vec([0.0*L 0.0*L 0.805*L])), inflate = tolerance)
susprn = selectnode(fens; box = initbox!(Float64[], vec([0.0*L -0.0*L 0.805*L])), inflate = tolerance)
suspbn = selectnode(fens; box = initbox!(Float64[], vec([-2.0*L 0.0*L 0.0*L])), inflate = tolerance)

# Find out at which nodes the sensors are:
sensors = Dict()
sensors = let 
    # The sensors at the tip of the left and right wing drum
    sensor112n = selectnode(fens; box = initbox!(Float64[], vec([+1.8*L 9.8*L 0.96*L])), inflate = tolerance)
    sensors[112] = sensor112n
    sensor12n = selectnode(fens; box = initbox!(Float64[], vec([+1.8*L -9.8*L .96*L])), inflate = tolerance)
    sensors[12] = sensor12n
    sensor111n = selectnode(fens; box = initbox!(Float64[], vec([-1.8*L 9.8*L 0.96*L])), inflate = tolerance)
    sensors[111] = sensor111n
    sensor11n = selectnode(fens; box = initbox!(Float64[], vec([-1.8*L -9.8*L .96*L])), inflate = tolerance)
    sensors[11] = sensor11n
        # The joint between the horizontal and vertical tail parts
    sensor202n = selectnode(fens; box = initbox!(Float64[], vec([-8*L 0 3.8*L])), inflate = tolerance)
    sensors[202] = sensor202n
    sensors
end

# ## Assemble the global discrete system


# For disambiguation we will refer to the stiffness and mass functions by
# qualifying them with the corotational-beam module, `FEMMCorotBeamModule`.
using FinEtoolsFlexStructures.FEMMCorotBeamModule
CB = FEMMCorotBeamModule

# Note that we have an array of finite element sets. We compute the matrices for
# each set separately and accumulate them into the final overall matrix. Thus
# we can construct the stiffness and mass matrix as follows.
using  SparseArrays

# Loop over all the finite element sets and add up their contributions, the
# stiffness in the mass matrix for each.
Kf, Kd, M = let
    Kf = spzeros(nfreedofs(dchi), nfreedofs(dchi))
    Kd = spzeros(nfreedofs(dchi), nfreedofs(dchi))
    M = spzeros(nfreedofs(dchi), nfreedofs(dchi))
    for fes in fesa
        labl  = fes.label[1]
        femm = CB.FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), getmaterial(labl));
        if labl == 7 # connectors representing the damping layer
            Kd += CB.stiffness(femm, geom0, u0, Rfield0, dchi);
        else
            Kf += CB.stiffness(femm, geom0, u0, Rfield0, dchi);
        end
        M += CB.mass(femm, geom0, u0, Rfield0, dchi);
    end
    Kf, Kd, M
end

# ## Additional concentrated masses.


using LinearAlgebra

using FinEtoolsFlexStructures.FEMMPointMassModule
PM = FEMMPointMassModule

# There is a sensor on the tail.
femmcm1 =  PM.FEMMPointMass(IntegDomain(FESetP1(reshape(sensors[202], 1, 1)), PointRule()), Matrix{Float64}(2*L*L/5*L/5*2*rho*LinearAlgebra.I(3)));

# These are the forward/interior locations on the wing drums where the
# compensation masses are attached.
mass1n = selectnode(fens; box = initbox!(Float64[], vec([1.8*L 9.2*L .96*L])), inflate = tolerance)
mass2n = selectnode(fens; box = initbox!(Float64[], vec([1.8*L -9.2*L .96*L])), inflate = tolerance)
femmcm2 =  PM.FEMMPointMass(IntegDomain(FESetP1(reshape([mass1n; mass2n;], 2, 1)), PointRule()), Matrix{Float64}(0.2*phun("kg")*LinearAlgebra.I(3)));

Mp = PM.mass(femmcm1, geom0, u0, Rfield0, dchi) + PM.mass(femmcm2, geom0, u0, Rfield0, dchi);

# ## Bungee supports 

# This is the assumed stifffness of the bungee cords (each one separately).
bungeecoefficient = 4000*phun("N/m");

using LinearAlgebra

using FinEtoolsFlexStructures.FEMMPointGroundedSpringModule
BS = FEMMPointGroundedSpringModule

# There are three suspension points at the top of the fuselage. We assume that these bungee supports exert only reaction in the vertical direction.
femmbs =  BS.FEMMPointGroundedSpring(IntegDomain(FESetP1(reshape([suspln; susprn; suspbn;], 3, 1)), PointRule()), 
    Matrix{Float64}([bungeecoefficient*[0;0;1]*[0;0;1]' 0*LinearAlgebra.I(3); 0*LinearAlgebra.I(3) 0*LinearAlgebra.I(3)]));

Kb = BS.stiffness(femmbs, geom0, u0, Rfield0, dchi)



# ## Damping

# The overall damping ratios measured in the physical experiments were
# approximately 1%. Damping levels thought appropriate for obtaining realistic
# damping were obtained through the use of a viscoelastic layer glued to the
# wing beam with an aluminum constraining plate on top. 


# ### Viscoelastic damping layer

# The properties of the material of the viscoelastic layer were described in the
# technical specification [3]. Since the properties are frequency dependent, we
# take as a representative value the numbers obtained for 20 Hz.

# The viscoelastic used was the 3M acrylic viscoelastic polymer ISD 112 in the
# form of a tape of 76 mm width with thickness of 50 microm.  This material is
# particularly well suited for the testbed operating range of 5-50 Hz and at 20
# degrees C where the loss factor is near its peak of 0.4 [3]. These are
# representative quantities taken at 20 Hz.

eta = 0.5
omega_5_50 = 2*pi*20


# The connectors between the wing beam and the constraining plate are taken as
# representative of the stiffness of the constraining layer. The viscoelastic
# damping is then taken as proportional to this stiffness.
Cd = (eta / omega_5_50) .* Kd

# ### Aluminum structure

# We assume the loss factor of the fuselage, wing, and tail to be 0.01. This can
# be used to derive the damping model for these parts of the aircraft structure
# in the form of Rayleigh damping.
zeta= 0.01/2; # damping ratio
zeta_1= zeta; # damping ratio for mode 1
zeta_2= zeta; # damping ratio for mode 2
# Stiffness and mass proportional damping parameters.
omega_1=2*pi*6; # Guess
omega_2=2*pi*30; # Guess
rdmass = (2*1.0/(omega_2^2-omega_1^2)).*[omega_2^2 -omega_1^2]*[zeta_1*omega_1; zeta_2*omega_2];
rdstiffness = (2*1.0/(omega_2^2-omega_1^2)).*[-1 1]*[zeta_1*omega_1; zeta_2*omega_2];

# The airframe damping matrix is a mixture of the mass matrix and the stiffness
# matrix.
Cf = rdmass .* M + rdstiffness .* Kf;


# These are the system matrices. The stiffness consists of the contribution of
# the airframe, the connectors between the wing components, and the bungees of
# the suspension.
Kt = Kf + Kd + Kb
# The mass matrix is the contribution of the structure and the attached point
# masses (compensation masses).
Mt = M + Mp
# The damping matrix is composed of the damping of the overall structure, and of
# the damping layer (the connectors between the wing components).
Ct = Cf + Cd

# ## Loading

# Here we assume that the stinger  was attached at the location of the sensor
# 12. The force is vertical (in the Z direction) and the magnitude of the force
# is arbitrary.
forceat = 12
fmagn= 1.0;
loadbdry = FESetP1(reshape(sensors[forceat], 1, 1))
lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
# The force is applied in the vertical direction, and we assume it is positive
# upwards.
fi = ForceIntensity(Float64[0, 0.0, -fmagn, 0, 0, 0]);

# The force loading is now integrated over the "volume" of the integration
# domain.
F = CB.distribloads(lfemm, geom0, dchi, fi, 3);

# ## Solve the Harmonic vibration problem

# The frequency sweep will start at `fromf` and continue through the frequency
# `tof`.  The frequencies will be logarithmically distributed throughout this
# range.
fromf = 3.0
tof = 70.0
nf = 150
frequencies = logspace(log10(fromf), log10(tof), nf);

# Only a single number per sensor will be collected, the Z direction
# displacement.
receptance12 = fill(0.0im, nf)
receptance112 = fill(0.0im, nf)

# Now Loop over 
for   fi in eachindex(frequencies)
    f =  frequencies[fi];
    om = 2*pi*f;
    # Solve the system of complex  equations of balance:
    U = (-om^2*Mt + 1im*om*Ct + Kt) \ F
    # Distribute the vector of the solution:
    scattersysvec!(dchi, U)
    # Now sample the solution at the locations of the two sensors, 12 and 112:
    v = fill(0.0im, 6)
    gathervalues_asvec!(dchi, v, sensors[12])
    p_d = v[1:3];
    receptance12[fi] = p_d[3]/fmagn;
    gathervalues_asvec!(dchi, v, sensors[112])
    p_d = v[1:3];
    receptance112[fi] = p_d[3]/fmagn;
end

# The receptances were calculated above. The mobility and accelerance may be
# obtained in a postprocessing step.
oms = (2*pi) .* frequencies;
mobility12 = receptance12 .* (-1im*oms);
mobility112 = receptance112 .* (-1im*oms);
accelerance12 = receptance12 .* (-oms.^2);
accelerance112 = receptance112 .* (-oms.^2);

results = Dict()
results[12] = Dict("receptance"=>receptance12, "mobility"=>mobility12, "accelerance"=>accelerance12)
results[112] = Dict("receptance"=>receptance112, "mobility"=>mobility112, "accelerance"=>accelerance112)

# ## Present the results graphically

using Gnuplot

# Plot the amplitude of the response curves. We output two curves.
# The first for the driving-point FRF:
quantity = "accelerance"; units = "m/s^2/N"
outputat = 12
y = abs.(results[outputat][quantity]) / phun(units)
@gp  "set terminal windows 0 "  :-
@gp  :- frequencies y " lw 2 lc rgb 'blue' with lines title 'output at $(outputat)' "  :-


# The second for the cross transfer:
outputat = 112
y = abs.(results[outputat][quantity]) / phun(units)
@gp  :- frequencies y " lw 2 lc rgb 'red' with lines title 'output at $(outputat)' "  :-
@gp  :- "set logscale y" :-
@gp  :- "set xlabel 'Frequency [Hz]'" :-
@gp  :- "set ylabel 'abs(H) [$(units)]'" :-
@gp  :- "set title 'Force at $(forceat), $(quantity)'"



# Plot the phase shift of the response curves. Again we output two curves, 
# the first for the driving-point FRF:
outputat = 12
y = atan.(imag(results[outputat][quantity]), real(results[outputat][quantity]))/pi*180 
@gp  "set terminal windows 1 "  :-
@gp  :- frequencies y " lw 2 lc rgb 'blue' with lines title 'output at $(outputat)' "  :-
# The second for the cross transfer:
outputat = 112
y = atan.(imag(results[outputat][quantity]), real(results[outputat][quantity]))/pi*180 
@gp  :- frequencies y " lw 2 lc rgb 'red' with lines title 'output at $(outputat)' "  :-
@gp  :- "set xlabel 'Frequency [Hz]'" :-
@gp  :- "set ylabel 'Phase shift [deg]'" :-
@gp  :- "set title 'Force at $(forceat), $(quantity)'"


nothing
