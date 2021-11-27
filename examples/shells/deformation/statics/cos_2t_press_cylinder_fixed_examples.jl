# Pressurized cylinder clamped around the ends.
# 
# Example introduced in
# @article{Lee2004,
#    author = {Lee, P. S. and Bathe, K. J.},
#    title = {Development of MITC isotropic triangular shell finite elements},
#    journal = {Computers & Structures},
#    volume = {82},
#    number = {11-12},
#    pages = {945-962},
#    ISSN = {0045-7949},
#    DOI = {10.1016/j.compstruc.2004.02.004},
#    year = {2004},
#    type = {Journal Article}
# }

module cos_2t_press_cylinder_fixed_examples

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

# Parameters:
E = 2.0e5
nu = 1/3;
pressure = 1.0;
R = 1.0;
L = 2.0;

# The cylinder axis is parallel to Y

function cylindrical!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) 
    n = cross(tangents[:, 1], tangents[:, 2]) 
    n = n/norm(n)
    # r = vec(XYZ); r[2] = 0.0
    csmatout[:, 3] .= n
    csmatout[:, 2] .= (0.0, 1.0, 0.0)
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

function computetrac!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    r = vec(XYZ); r[2] = 0.0
    r .= vec(r)/norm(vec(r))
    theta = atan(r[3], r[1])
    forceout[1:3] = r*pressure*cos(2*theta)
    forceout[4:6] .= 0.0
    # @show dot(n, forceout[1:3])
    return forceout
end

function _execute_dsg_model(formul, n = 8, thickness = R/100, visualize = true)
    tolerance = R/n/1000
    fens, fes = T3block(90/360*2*pi,L/2,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*cos(a))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # plane of symmetry perpendicular to Z
    l1 = selectnode(fens; box = Float64[-Inf Inf -Inf Inf 0 0], inflate = tolerance)
    for i in [3,4,5]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    # clamped edge perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Midpoint of the fixed edge
    # nl = selectnode(fens; box = Float64[R R L/2 L/2 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    
    fi = ForceIntensity(FFlt, 6, computetrac!);
    F = distribloads(lfemm, geom0, dchi, fi, 2);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    strainenergy = 1/2 * U' * K * U
    @info "Strain Energy: $(round(strainenergy, digits = 9))"

    # Generate a graphical display of resultants
    ocsys = CSys(3, 3, cylindrical!)
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("m$nc", fld.values))
    end
    vtkwrite("cos_2t_press_cylinder_fixed-$(n)-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("n$nc", fld.values))
    end
    vtkwrite("cos_2t_press_cylinder_fixed-$(n)-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("q$nc", fld.values))
    end
    vtkwrite("cos_2t_press_cylinder_fixed-$(n)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

    # Visualization
    if visualize
        scattersysvec!(dchi, (L/8)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
            #plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end

    return strainenergy
end

function test_convergence(formul, thickness = R/100)
    @info "Pressurized Cylindrical shell, fixed ends, formulation=$(formul)"
    for n in [8, 16, 32, 64, 128]
        _execute_dsg_model(formul, n, thickness, false)
    end
    return true
end

end # module

using .cos_2t_press_cylinder_fixed_examples
using FinEtoolsFlexStructures.FEMMShellT3FFModule
cos_2t_press_cylinder_fixed_examples.test_convergence(FEMMShellT3FFModule)
