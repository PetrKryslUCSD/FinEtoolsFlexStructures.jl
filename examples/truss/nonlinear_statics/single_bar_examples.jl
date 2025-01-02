"""

"""
module single_bar_examples

using FinEtools
using FinEtools.AssemblyModule: SysmatAssemblerFFBlock, SysvecAssemblerFBlock
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtoolsFlexStructures.FEMMCorotTrussModule
using FinEtoolsFlexStructures.FEMMCorotTrussModule: FEMMCorotTruss
stiffness = FEMMCorotTrussModule.stiffness
mass = FEMMCorotTrussModule.mass
geostiffness = FEMMCorotTrussModule.geostiffness
restoringforce = FEMMCorotTrussModule.restoringforce
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json

function solve(visualize=false)
    # Parameters:
    E = 30_000.0 * phun("MPa")
    nu = 0.3
    W = 2500.0 * phun("mm")
    H = 250.0 * phun("mm")
    x = Float64[
        0.0 0.0 0.0 # 1
        W H 0.0  # 2
    ] 
    area = 5e7 / E
    P = 1.0 * phun("kilo*N")
    maxit = 10

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> sqrt(area), s -> sqrt(area), s -> [0.0, 0.0, -1.0])

    # Select the number of elements per leg.
    n = 1
    members = []
    for j in (
        [1, 2],
    )
        push!(members, frame_member(x[j, :], n, cs))
    end
    fens, fes = merge_members(members; tolerance=1 / 1000)

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, E, nu)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 3))

    # Apply EBC's
    pinned = selectnode(fens; box = boundingbox(x[1, :]), tolerance = H/100)
    for i in [1, 2, 3,]
        setebc!(dchi, pinned, true, i)
    end
    loaded = selectnode(fens; box = boundingbox(x[2, :]), tolerance = H/100)
    for i in [1, 3,]
        setebc!(dchi, loaded, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)
    
    # Assemble the global discrete system
    # massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    # vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    femm = FEMMCorotTruss(IntegDomain(fes, GaussRule(1, 2)), material)
    
    K = stiffness(femm, geom0, u0, Rfield0, dchi)
    loadbdry = FESetP1(reshape(loaded, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = zeros(3)
    q[2] = P
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3)
    
    # Solve the static problem
    fr = freedofs(dchi)
    
    # Auxiliary Variables 
    u1 = deepcopy(u0)
    rhs = gathersysvec(dchi, DOF_KIND_FREE);
    utol = 1e-13*nfreedofs(dchi);

    F = F[fr]

    load_parameters = 0.0:25:100;
    tip_displacement = fill(0.0, length(load_parameters), 1);
    step = 0
    for load_parameter in load_parameters
        applyebc!(dchi) # Apply boundary conditions
        u1.values[:] = u0.values[:]; # guess
        
        println("Load: $load_parameter")
        iter = 1;
        while true
            Fr = restoringforce(femm, geom0, u1, Rfield0, dchi);       # Internal forces
            rhs .= load_parameter .* F + Fr[fr];
            K = stiffness(femm, geom0, u1, Rfield0, dchi) + geostiffness(femm, geom0, u1, Rfield0, dchi);
            dchi = scattersysvec!(dchi, (K[fr, fr])\rhs); # Disp. incr
            u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
            print("$iter: ||du||=$(maximum(abs.(dchi.values[:])))\n")
            if maximum(abs.(dchi.values[:])) < utol # convergence check
                break;
            end
            if (iter > maxit)# bailout for failed convergence
                error("Possible failed convergence");
            end
            iter += 1;
        end
        u0.values[:] = u1.values[:];       # update the displacement
        
        step = step + 1
        
        tip_displacement[step] = u1.values[loaded[1], 2]
    end
    

    if visualize
        scaling = 1e1
        dchi.values .*= scaling
        radius = 20 * phun("in")
        plots = cat(plot_space_box([[-radius -radius -radius]; [radius radius radius]]),
                    plot_nodes(fens),
                    plot_solid(fens, fes;
                               x = geom0.values, u = u1.values[:, 1:3], );
                    dims = 1)
        pl = render(plots)
    end
    
    true
end # force_1

function allrun()
    println("#####################################################")
    println("# solve ")
    solve(true)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
