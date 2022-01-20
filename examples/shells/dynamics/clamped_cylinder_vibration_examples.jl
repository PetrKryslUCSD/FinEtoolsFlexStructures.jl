"""
Clamped cylinder, free vibration.

Reference
http://www.adina.com/newsgH53.shtml
"""
module clamped_cylinder_vibration_examples

using Arpack
using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

E = 200e9;
nu = 0.3;
rho = 7800.0
R = 0.1
L = 0.4

cylindrical!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) = begin
    r = vec(XYZ); r[2] = 0.0;
    csmatout[:, 3] .= vec(r)/norm(vec(r))
    csmatout[:, 2] .= (0.0, 1.0, 0.0) #  this is along the axis
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

function _execute(n = 2, thickness = 0.01, visualize = true)
    formul = FEMMShellT3FFModule

    n = Int(round(n/2)*2)
    @info "Mesh: $n elements per side"

    # Mesh
    tolerance = min(R, L) / n  / 100
    
    tolerance = R/n/1000
    fens, fes = T3block(360.0, L, n, n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a/180*pi), y-L/2, R*cos(a/180*pi))
    end
    fens, fes = mergenodes(fens, fes, tolerance)
    bfes = meshboundary(fes)
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    
    associategeometry! = formul.associategeometry!
    stiffness = formul.stiffness
    mass = formul.mass

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Clamp the cross sections
    nl = connectednodes(bfes)
    for i in 1:6
        setebc!(dchi, nl, true, i)
    end

    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    M = mass(femm, geom0, dchi);

    # Solve
    OmegaShift = 0.0*2*pi
    neigvs = 20
    evals, evecs, nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    evals[:] = evals .- OmegaShift;
    fs = real(sqrt.(complex(evals)))/(2*pi)
    @show fs

    for ev in 1:length(fs)
        U = evecs[:, ev]
        scattersysvec!(dchi, 1.0/maximum(abs.(U)).*U)
        vtkwrite("clamped_cylinder_vibration-$thickness-mode-$(ev).vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
    end

    return true
end

function test_convergence()
    @info "Clamped cylinder vibration"
    for n in [64] # 3:64 #
        _execute(n, 0.01, !false)
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

