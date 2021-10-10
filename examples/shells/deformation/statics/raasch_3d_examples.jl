"""
The geometry consists of a “hook” in the form of a curved strip rigidly clamped
at one end and loaded with a unit in-plane shear along the width at the other
end. It has two circular segments that are connected at the tangent point. The
smaller segment has a mean radius of 0.3556 m (14 inches) and spans 60° from
the clamped end to the tangent point. The larger segment spans 150° from the
tangent point to the free end and has a mean radius of 1.1684 m (46 inches).
The hook is 0.0508 m (2 inches) thick and 0.508 m (20 inches) wide, modeled as
linear elastic with an elastic modulus of 22.77 MPa (3300 psi) and a Poisson's
ratio of 0.35. In most tests the shear force is applied through the use of a
distributing coupling constraint. The coupling constraint provides coupling
between a reference node on which the load is prescribed and the nodes located
on the free end. The distributed nodal loads on the free end are equivalent to
a uniformly distributed load of 8.7563 N/m (0.05 lb/in). In two of the tests an
equivalent shear force is applied as a distributed shear traction instead.


"""
module raasch_3d_examples

using LinearAlgebra, Statistics
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

using Infiltrator

function _execute(input = "raasch_s4_1x9.inp", drilling_stiffness_scale = 1.0, visualize = true, nL = 9, nW = 1, nT = 4)
    E = 3300.0;
    nu = 0.35;
    thickness  =  2.0;
    tolerance = min(46/nL/100, 20/nW/100, thickness/nT/10)
    # analytical solution for the vertical deflection under the load
    analyt_sol = 5.09;
    R = 46.0;

    fens, fes = H8block(210.0, 20.0, thickness, nL, nW, nT)
    fens.xyz[:, 3] .+= - thickness/2
    r = (r0, t) -> r0 + t 
    for k in 1:count(fens)
        a, z, t = fens.xyz[k, :]
        if a <= 60.0
            fens.xyz[k, :] .= (r(14, t)*sin(a*pi/180), (14-r(14, t)*cos(a*pi/180)), z)
        else
            a = a - 60.0 + 30.0
            xc = 14*sin(60*pi/180) + 46*cos(30*pi/180)
            yc = 14*(1-cos(60*pi/180)) - 46*sin(30*pi/180)
            fens.xyz[k, :] .= (xc-r(46, -t)*cos(a*pi/180), yc+r(46, -t)*sin(a*pi/180), z)
        end
    end
    input = "raasch_3d_$(nL)x$(nW)x$(nT)"

    @show count(fens), count(fes)

    vtkwrite("$input.vtu", fens, fes)

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, E, nu)
    
    @info "Mesh: $input"

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))

    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)
    femm = associategeometry!(femm, geom)

    # Apply EBC's
    # Clamped end
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,]
        setebc!(u0, l1, true, i)
    end
    
    applyebc!(u0)
    numberdofs!(u0);

    # Assemble the system matrix
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u0);

    # Load
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, box = [0 Inf -16 -16 0 20], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(2, 2)))
    fi = ForceIntensity(FFlt[0, 0, 0.05/thickness, 0, 0, 0]);
    F = distribloads(lfemm, geom, u0, fi, 2);
    
    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(u0, U[:])
    nl = selectnode(fens; box = Float64[0 Inf -16 -16 0 20], inflate = tolerance)
    targetu =  mean(u0.values[nl, 3])
    @info "Solution: $(round(targetu, digits=8)),  $(round(targetu/analyt_sol, digits = 4)*100)%"

    return targetu/analyt_sol
end

function test_convergence()
    
    @info "Raasch hook, 3D mesh"

    for n in 1:7
        _execute("", 1.0, false, 9*2^(n-1), 2^(n-1)+1, 2+(n-1))
    end
    return true
end

end # module

using .raasch_3d_examples
raasch_3d_examples.test_convergence()
