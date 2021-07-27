"""

"""
module single_element_examples

using LinearAlgebra
using Arpack
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3, local_frame!
using FinEtoolsFlexStructures.FEMMShellDSG3Module: FEMMShellDSG3, stiffness
    # using FinEtoolsFlexStructures.FEMMShellT3Module: FEMMShellT3, stiffness
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json, plot_triads

function single_dsg3()
    E = 200e3*phun("MPa")
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    thickness = 0.1*phun("m");
    L = 1.0*phun("m");

    # Mesh
    tolerance = L/1000
    fens, fes = T3block(L,L,1,1);
    fes = subset(fes, [1])
    fens.xyz = xyz3(fens)

    connected = findunconnnodes(fens, fes);
    fens, new_numbering = compactnodes(fens, connected);
    fes = renumberconn!(fes, new_numbering);
    @show fens.xyz
    @show fes

    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = FEMMShellDSG3(IntegDomain(fes, TriRule(1), thickness), mater)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # No EBC's
    l1 = collect(1:count(fens))
    for i in [1, 2, 3, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    @show dchi.dofnums

    # Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    K = Matrix(K)
    dec = eigen(K)
    @show dec.values
    v = dec.vectors

    # mfemm = FEMMDeforLinear(DeforModelRed3D, IntegDomain(fes, TriRule(1), thickness), mater)
    # M = mass(mfemm, geom0, dchi);

    # # Solve
    # OmegaShift = 0.1*2*pi
    # neigvs = 10
    # d, v, nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    # d[:] = d .- OmegaShift;
    # fs = real(sqrt.(complex(d)))/(2*pi)
    # @show fs

    for j in 1:6
        @show round.(v[:, j], digits=4)
    end
        
    # # Visualization
    for j in 1:length(dec.values)
        U = v[:, j]
        scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[-L -L -L]; [+L +L +L]]),
            plot_nodes(fens),
            plot_triads(fens; triad_length = 0.4, x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values)
            ;
            dims = 1)
        pl = render(plots, title = "$j: Eigenvalue $(dec.values[j])")
    end
end


function allrun()
    println("#####################################################")
    println("# single_dsg3 ")
    single_dsg3()
    return true
end # function allrun

end # module

using .single_element_examples
single_element_examples.allrun()