"""

"""
module single_element_examples

using LinearAlgebra
using Arpack
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3, local_frame!
using FinEtoolsFlexStructures.FEMMShellDSG3Module
using FinEtoolsFlexStructures.FEMMShellDSG3IModule
using FinEtoolsFlexStructures.FEMMShellDSG3IFModule
using FinEtoolsFlexStructures.FEMMShellCSDSG3Module
using FinEtoolsFlexStructures.FEMMShellT3Module
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json, plot_triads

function standard_single_dsg3()
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
    
    formul = FEMMShellCSDSG3Module
    formul = FEMMShellDSG3IFModule

    # @show Bb[:, [4, 5, 10, 11, 16, 17]]
    # @show Bs[:, [4, 5, 10, 11, 16, 17]]
    
    # Bb[:, [4, 5, 10, 11, 16, 17]] = [
    # 0.0 -1.0 0.0 1.0 0.0 0.0; 
    # 1.0 0.0 0.0 0.0 -1.0 0.0; 
    # 1.0 -1.0 -1.0 0.0 0.0 1.0]                                         
    # Bs[:, [4, 5, 10, 11, 16, 17]] = [
    # -1/6 1/3 1/6 1/2 0.0 1/6; 
    # -1/3 1/6 -1/6 0.0 -1/2 -1/6]  

    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    associate = formul.associategeometry!
        stiffness = formul.stiffness
        mass = formul.mass

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
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    K = Matrix(K)
    dec = eigen(K)
    @show dec.values
    v = dec.vectors

    # mfemm = FEMMDeforLinear(DeforModelRed3D, IntegDomain(fes, TriRule(1), thickness), mater)
    # M = mass(mfemm, geom0, dchi);

    # # Solve1
    # OmegaShift = 0.1*2*pi
    # neigvs = 10
    # d, v, nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    # d[:] = d .- OmegaShift;
    # fs = real(sqrt.(complex(d)))/(2*pi)
    # @show fs

    for j in 1:6
        @show dec.values[j]
        @show round.(v[:, j], digits=4)
    end
        
    # # Visualization
    for j in 1:length(dec.values)
        U = v[:, j]
        scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
        Rfield = deepcopy(Rfield0)
        update_rotation_field!(Rfield, dchi)
        plots = cat(plot_space_box([[-L -L -L]; [+L +L +L]]),
            plot_nodes(fens),
            plot_triads(fens; triad_length = 0.2, x = geom0.values, u = dchi.values[:, 1:3], R = Rfield.values),
            plot_triads(fens; triad_length = 0.4, x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values)
            ;
            dims = 1)
        pl = render(plots, title = "$j: Eigenvalue $(dec.values[j])")
    end
end


function distorted_single_dsg3()
    E = 200e3*phun("MPa")
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    thickness = 0.1*phun("m");
    L = 1.0*phun("m");

    # Mesh
    tolerance = L/1000
    fens, fes = T3block(L,L,1,1);
    fes = subset(fes, [1])
    fens.xyz = [
    0.0 -10.0 0.0; 
    13.0 -3.0 0.0;
    -1.0 7.0 0.0] 
    fens.xyz = xyz3(fens)

    connected = findunconnnodes(fens, fes);
    fens, new_numbering = compactnodes(fens, connected);
    fes = renumberconn!(fes, new_numbering);
    @show fens.xyz
    @show fes

    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    
    formul = FEMMShellCSDSG3Module
    formul = FEMMShellDSG3IFModule

    # @show Bb[:, [4, 5, 10, 11, 16, 17]]
    # @show Bs[:, [4, 5, 10, 11, 16, 17]]
    
    # Bb[:, [4, 5, 10, 11, 16, 17]] = [
    # 0.0 -1.0 0.0 1.0 0.0 0.0; 
    # 1.0 0.0 0.0 0.0 -1.0 0.0; 
    # 1.0 -1.0 -1.0 0.0 0.0 1.0]                                         
    # Bs[:, [4, 5, 10, 11, 16, 17]] = [
    # -1/6 1/3 1/6 1/2 0.0 1/6; 
    # -1/3 1/6 -1/6 0.0 -1/2 -1/6]  

    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    associate = formul.associategeometry!
        stiffness = formul.stiffness
        mass = formul.mass

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
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    K = Matrix(K)
    dec = eigen(K)
    @show dec.values
    v = dec.vectors

    # mfemm = FEMMDeforLinear(DeforModelRed3D, IntegDomain(fes, TriRule(1), thickness), mater)
    # M = mass(mfemm, geom0, dchi);

    # # Solve1
    # OmegaShift = 0.1*2*pi
    # neigvs = 10
    # d, v, nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    # d[:] = d .- OmegaShift;
    # fs = real(sqrt.(complex(d)))/(2*pi)
    # @show fs

    for j in 1:6
        @show dec.values[j]
        @show round.(v[:, j], digits=4)
    end

    box = reshape(boundingbox(fens.xyz), 1, 6)
    box = vcat(box[1, [1, 3, 5]]' .+ [-5, -5, -5], box[1, [2, 4, 6]]' .- [-5, -5, -5])    
    # # Visualization
    for j in 1:length(dec.values)
        U = v[:, j]
        scattersysvec!(dchi, (L/4)/maximum(abs.(U)).*U)
        Rfield = deepcopy(Rfield0)
        update_rotation_field!(Rfield, dchi)
        plots = cat(plot_space_box(box),
            plot_nodes(fens),
            plot_triads(fens; triad_length = 2.2, x = geom0.values, u = dchi.values[:, 1:3], R = Rfield.values),
            plot_triads(fens; triad_length = 3.4, x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values),
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
single_element_examples.distorted_single_dsg3()