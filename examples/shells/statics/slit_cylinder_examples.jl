"""
Slit cylinder under torsional loads.

Reference
Morris, A. J., NAFEMS: shell finite element evaluation tests
Reference P12
Section 3.3 In Extensional Bending
Slit Cylinder Under Applied Torque T = 4*pi*D*(1-nu)*C/R
x = axial coordinate
displacements: axial u = -C phi, circumferential v = C x, radial w = 0,
twisting moment M_x,phi =  D*(1-nu)*C/R^2
all other stress resultants (section forces and moments) are zero.

Note: our coordinates are different. The axial coordinate is y.
"""
module slit_cylinder_examples

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

E = 3e6;
nu = 0.3;
R = 300.0;
L = 3*R;
thickness = 3.0;

# T = 4*pi*D*(1-nu)*C/R  -> C = T*R/(4*pi*D*(1-nu))
T = 1000.0
D = E*thickness^3/12/(1-nu^2)
C = T*R/(4*pi*D*(1-nu))
# The maximum displacement axially is:
vmax = C * pi
# Twisting moment magnitude
M_x_phi =  D*(1-nu)*C/R^2

cylindrical!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) = begin
    r = vec(XYZ); r[2] = 0.0;
    csmatout[:, 3] .= vec(r)/norm(vec(r))
    csmatout[:, 2] .= (0.0, 1.0, 0.0) #  this is along the axis
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

function _execute(n = 2, visualize = true)
    formul = FEMMShellT3FFModule

    n = Int(round(n/2)*2)
    @info "Mesh: $n elements per side"

    # Mesh
    tolerance = min(R, L) / n  / 100
    
    tolerance = R/n/1000
    fens, fes = T3block(360.0, L, n, n);
    bfes = meshboundary(fes)
    al0 = selectelem(fens, bfes, box = [0.0 0.0 -Inf Inf], inflate = tolerance)
    al360 = selectelem(fens, bfes, box = [360.0 360.0 -Inf Inf], inflate = tolerance)
    ll0 = selectelem(fens, bfes, box = [-Inf Inf 0.0 0.0], inflate = tolerance)
    llL2 = selectelem(fens, bfes, box = [-Inf Inf L L], inflate = tolerance)
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a/180*pi), y-L/2, R*cos(a/180*pi))
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    femm.drilling_stiffness_scale = 0.1
    associategeometry! = formul.associategeometry!
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Fix one of the nodes to remove rigid body displacements
    @show     nl = selectnode(fens, nearestto = Float64[0.0, L/2, -R])
    for i in 1:6
        setebc!(dchi, nl, true, i)
    end

    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    loadbdry = subset(bfes, ll0)
    lfemm = FEMMBase(IntegDomain(loadbdry, GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, 0, 0, 0, T/(2*pi*R), 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 1);
    loadbdry = subset(bfes, llL2)
    lfemm = FEMMBase(IntegDomain(loadbdry, GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[0, 0, 0, 0, -T/(2*pi*R), 0]);
    F += distribloads(lfemm, geom0, dchi, fi, 1);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])

    @show vmax, maximum(dchi.values[connectednodes(subset(bfes, llL2)), 2])

    fld = fieldfromintegpoints(femm, geom0, dchi, :moment, 3, outputcsys = ocsys)
    @show -M_x_phi, minimum(fld.values[:]), maximum(fld.values[:])


        # Generate a graphical display of resultants
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
        push!(scalars, ("m$nc", fld.values))
        fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
        push!(scalars, ("em$nc", fld.values))
    end
    vtkwrite("slit_cylinder-$(n)-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
        push!(scalars, ("n$nc", fld.values))
        fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
        push!(scalars, ("en$nc", fld.values))
    end
    vtkwrite("slit_cylinder-$(n)-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
        push!(scalars, ("q$nc", fld.values))
        fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
        push!(scalars, ("eq$nc", fld.values))
    end
    vtkwrite("slit_cylinder-o-$(n)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

    vtkwrite("slit_cylinder-$(n)-uur.vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])

    return true
end

function test_convergence()
    @info "Slit cylinder"
    for n in [8, 16, 32] # 3:64 #
        _execute(n, !false)
    end
    return true
end

function allrun()
    println("#####################################################")
    println("# test_convergence ")
    test_convergence()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing

