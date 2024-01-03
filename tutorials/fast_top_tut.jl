# #  Fast Lagrangian top.

# Source code: [`fast_top_tut.jl`](fast_top_tut.jl)

# ## Description

# Fast-spinning Lagrange  top simulated by a beam model. The reference
# solution is available in Krysl, P., Endres, L.: Explicit Newmark/Verlet 
# algorithm for time integration of the rotational dynamics of rigid bodies.
# Int. J. Numer. Meth. Eng. 62, 2154–2177 (2005).

# A beam is pinned at the bottom, and spins in a gravity field, executing
# precession and nutation. Its top node inscribes a curve in the 3D space which
# is displayed as the equations of motion are integrated.

# ![](fast_top.png)


# ## Goals

# - Illustrate integration of the nonlinear equations 
#   of motion with the Newmark algorithm.
# - Show the power of online visualization.

using LinearAlgebra
using PlotlyJS

# ## Definition of the basic inputs

# The finite element code realize on the basic functionality implemented in this
# package.

using FinEtools

# The material parameters may be defined with the specification of the units.
# The elastic properties and the mass density are:
E = 71240.0 * phun("MPa")
nu = 0.31
rho = 2.7e3 * phun("kg/m^3")

# The top is a block of square cross-section with dimensions
Width = 60 * phun("mm"); Height = 60 * phun("mm"); 
Length = 4*Width

Mass = Width*Height*Length*rho;
IL =  (1/12)*Mass*(Width^2+Height^2);
Ib =  (1/12)*Mass*(Length^2+Height^2)+Mass*(Length/2)^2;
Ih =  (1/12)*Mass*(Width^2+Length^2)+Mass*(Length/2)^2;
Omega0 = 313*pi;
R0 = rotmat3([0.05 0 0]);
g = 9.81 * phun("m/sec^2"); 
q = [0,0,-g*Width*Height*rho];
utol = 1e-3 * phun("m");
maxit = 12
dt = min(2*pi/norm(Omega0)/10, 0.005);
tend = 0.8 * phun("sec");
ng = 1/2; nb = 1/4*(1/2+ng)^2;
# Choose the mass formulation:
mass_type=1;


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
cs = CrossSectionRectangle(s -> Width, s -> Width, s -> [1.0, 0.0, 0.0])

# Select the number of elements per leg.
spin_vector = R0*[0, 0, Omega0];
X=[0 0 0;    reshape(R0*[0,0,Length], 1, 3)];
n=2;
tolerance=Length/n/100;

using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member
members = []
push!(members, frame_member(X, n, cs))

using FinEtoolsFlexStructures.MeshFrameMemberModule: merge_members
fens, fes = merge_members(members; tolerance = tolerance);

# Material properties
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
dchi = NodalField(zeros(size(fens.xyz,1), 6))

# ## Support conditions

# The "bottom" of the top is pinned.
box = fill(0.0, 6)
initbox!(box, X[1,:])
supportn = selectnode(fens; box = box, inflate = tolerance)
initbox!(box, X[2,:])
tipn = selectnode(fens; box = box, inflate = tolerance)
for i in [1, 2, 3]
    setebc!(dchi, supportn, true, i)
end
applyebc!(dchi)
numberdofs!(dchi);

# Initial conditions
v0 = deepcopy(dchi) # need the numbering of the degrees of freedom
for i in 1:size(v0.values, 1)
    v0.values[i, 4:6] = spin_vector
end


# For disambiguation we will refer to the stiffness and mass functions by qualifying them with the corotational-beam module, `FEMMCorotBeamModule`.
using FinEtoolsFlexStructures.FEMMCorotBeamModule
CB = FEMMCorotBeamModule
femm = CB.FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
fi = ForceIntensity(q);

# Here we prepare matrix and vector assemblers that only work with a free
# degrees of freedom.
massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using DelimitedFiles

function integrate(tend, CB, geom0, u0, Rfield0, dchi, v0, report)
    # Make sure we don't clobber any of these variable fields
    geom0, u0, Rfield0, dchi, v0 = deepcopy((geom0, u0, Rfield0, dchi, v0))
    # Additional fields
    stepdchi = deepcopy(dchi)
    u1 = deepcopy(u0)
    v1 = deepcopy(dchi)
    a1 = deepcopy(dchi) # zero out the Acceleration
    a0 = deepcopy(dchi) # zero out the Acceleration
    Rfield1 = deepcopy(Rfield0)
    v0v = gathersysvec(v0)
    a0v = gathersysvec(a0)
    vpv = gathersysvec(v0);
    dchipv = gathersysvec(dchi);
    stepdchiv = gathersysvec(dchi);
    rhs = gathersysvec(dchi);
    TMPv = deepcopy(rhs)
    utol = 1e-13*nfreedofs(dchi);
    
    t = 0.0;
    step = 0;
    while (t <= tend)
        t = t + dt;
        (mod(step, 50)==0) && println("Time $(t)"); # pause
        # Initialization
        applyebc!(dchi) # Apply boundary conditions
        u1.values[:] = u0.values[:]; # guess
        Rfield1.values[:] = Rfield0.values[:]; # guess
        stepdchi.values[:] .= 0.0;# Total increment in current step
        a1.values[:] = -(1/nb/dt)*v0.values[:] -(1/2-nb)/nb*a0.values[:];
        v1.values[:] = v0.values[:] + dt*((1-ng)*a0.values[:] + ng*a1.values[:]);
        gathersysvec!(v0, v0v)
        gathersysvec!(a0, a0v)
        dchipv = dt*v0v + (dt^2/2*(1-2*nb))*a0v
        vpv = v0v +(dt*(1-ng))*a0v;
        
        iter = 1;
        while true
            F = CB.distribloads_global(femm, vassem, geom0, u1, Rfield1, dchi, fi) +
             CB.restoringforce(femm, vassem, geom0, u1, Rfield1, dchi);
            rhs .= F
            K = CB.stiffness(femm, massem, geom0, u1, Rfield1, dchi);
            M = CB.mass(femm, massem, geom0, u1, Rfield1, dchi);
            G = CB.gyroscopic(femm, massem, geom0, u1, Rfield1, v1, dchi);
            gathersysvec!(stepdchi, stepdchiv)
            @. TMPv = ((-1/(nb*dt^2))*stepdchiv+(1/(nb*dt^2))*dchipv)
            rhs .+= M*TMPv
            @. TMPv = ((-ng/nb/dt)*stepdchiv+(ng/nb/dt)*dchipv - vpv)
            rhs .+= G*TMPv;
            dchi = scattersysvec!(dchi, (K+(ng/nb/dt)*G+(1/(nb*dt^2))*M)\rhs); # Disp. incr
            u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
            stepdchi.values[:] += dchi.values[:]
            v1.values[:] += (ng/nb/dt)*dchi.values[:];
            a1.values[:] += (1/nb/dt^2)*dchi.values[:];
            update_rotation_field!(Rfield1, dchi)
            if maximum(abs.(dchi.values[:])) < utol# convergence check
                break; 
            end
            if (iter > maxit)# bailout for failed convergence
                error("Possible failed convergence");
            end
            iter += 1;
        end
        u0.values[:] = u1.values[:];       # update the displacement
        Rfield0.values[:] = Rfield1.values[:]; # update the rotations
        v0.values[:] = v1.values[:];       # update the velocities
        a0.values[:] = a1.values[:];       # update the accelerations
        
        report(step, u1, Rfield1)
        
        step=step+1;
    end
end

# The visualization utilities take advantage of the PlotlyJS library.
using PlotlyJS
using VisualStructures: plot_space_box, plot_midline, plot_solid, render, react!, default_layout_3d, save_to_json

# Display the graph of the motion of the tip of the top.

function updategraph(step, u1, Rfield1)
    if (mod(step,20)==0)
        push!(tipx, X[2,1]+u1.values[tipn[1], 1])
        push!(tipy, X[2,2]+u1.values[tipn[1], 2])
        plots = cat(tbox, tref, scatter(;x=tipx./Length, y=tipy./Length, mode="markers", name = "Sol", line_color = "rgb(155, 15, 15)"); dims = 1)
        react!(pl, plots, pl.plot.layout)
        sleep(0.01)
    end
end

# integrate(tend, CB, geom0, u0, Rfield0, dchi, v0, updategraph)

tipx = Float64[]
tipy = Float64[] 
tipz = Float64[] 
layout = Layout(; scene=attr(
    xaxis = attr(title="X"),
    yaxis = attr(title="Y"),
    zaxis = attr(title="Z"),
    camera = attr(
        up=attr(x=-0.8, y=-0.6, z=0.07),
        center=attr(x=0.0, y=0.0, z=0.0),
        eye=attr(x=0.1, y=0.12, z=2.17),
        projection = attr(type = "orthographic")
        )), showlegend = false)
tbox = plot_space_box([[-1.1*Width -1.1*Width 0]; [1.1*Width 1.1*Width 1.1*Length]])
tshape0s = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
tshape0m = plot_midline(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], color = "rgb(125, 105, 175)", lwidth = 4)
plots = cat(tbox,  tshape0s, tshape0m; dims = 1)
pl2 = render(plots; layout = layout)
sleep(0.5)

function updateplot(step, u1, Rfield1)
    push!(tipx, X[2,1]+u1.values[tipn[1], 1])
    push!(tipy, X[2,2]+u1.values[tipn[1], 2])
    push!(tipz, X[2,3]+u1.values[tipn[1], 3])
    if (mod(step,13)==0)
        curv = scatter3d(;x=tipx, y=tipy, z=tipz, mode="lines", name = "Sol", line_color = "rgb(155, 15, 15)")
        tshape1s = plot_solid(fens, fes; x = geom0.values, u = u1.values, R = Rfield1.values, facecolor = "rgb(125, 15, 15)");
        tshape1m = plot_midline(fens, fes; x = geom0.values, u = u1.values, color = "rgb(125, 15, 15)", lwidth = 4)
        plots = cat(tbox,  tshape0s, tshape0m,  curv, tshape1s, tshape1m; dims = 1)
        react!(pl2, plots, pl2.plot.layout)
        sleep(0.12)
    end
end


integrate(4.5*tend, CB, geom0, u0, Rfield0, dchi, v0, updateplot)
