module mmodal1
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMCorotBeamModule
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
mass = FEMMCorotBeamModule.mass
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
function test()
        # Parameters:
        E=71240.0;#MPa
        nu=0.31;# Poisson ratio
        rho=5e-9;
        b=3.0; h=30.0; L=240.0; # cross-sectional dimensions and length of each leg in millimeters
        # Choose the mass formulation:
        mass_type=2;
        scale = 0.4

        # Reference frequencies
        reffs = [11.2732, 30.5269]
        neigvs = 2;

        # Cross-sectional properties
        cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 1.0, 0.0])

        # Select the number of elements per leg.
        xyz =
        n = 4;
        members = []
        push!(members, frame_member([0 0 L; L 0 L], n, cs))
        push!(members, frame_member([L 0 L; L 0 0], n, cs))
        fens, fes = merge_members(members; tolerance = L / 10000);

        # Material properties
        material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

        # Construct the requisite fields, geometry and displacement
        # Initialize configuration variables
        geom0 = NodalField(fens.xyz)
        u0 = NodalField(zeros(size(fens.xyz,1), 3))
        Rfield0 = initial_Rfield(fens)
        dchi = NodalField(zeros(size(fens.xyz,1), 6))

        # Apply EBC's
        l1 = selectnode(fens; box = [0 0 0 0 L L], tolerance = L/10000)
        for i in [1,2,3,4,5,6]
            setebc!(dchi, l1, true, i)
        end
        applyebc!(dchi)
        numberdofs!(dchi);

        # Assemble the global discrete system
        femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
        K = stiffness(femm, geom0, u0, Rfield0, dchi);
        M = mass(femm, geom0, u0, Rfield0, dchi);
        K_ff = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:ff]
        M_ff = matrix_blocked(M, nfreedofs(dchi), nfreedofs(dchi))[:ff]

        # Solve the eigenvalue problem
        d,v,nev,nconv = eigs(K_ff, M_ff; nev=2*neigvs, which=:SM)
        fs = real(sqrt.(complex(d)))/(2*pi)
        # println("Natural frequencies: $fs [Hz]")
        # println("Reference: $reffs [Hz]")

        @test norm(reffs - fs[1:length(reffs)]) ./ norm(reffs) <= 3.0e-5

        # Visualize vibration modes
        # mode = 2
        # tbox = plot_space_box([[-L/2 -L/2 0]; [L/2 L/2 1.1*L]])
        # tenv0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
        # plots = cat(tbox, tenv0; dims = 1)
        # pl = render(plots; title = "Mode $(mode)")
        # Rfield1 = deepcopy(Rfield0)
        # for xscale in scale.*sin.(collect(0:1:72).*(2*pi/21))
        #     scattersysvec!(dchi, xscale.*v[:, mode])
        #     Rfield1 = deepcopy(Rfield0)
        #     update_rotation_field!(Rfield1, dchi)
        #     tenv1 = plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield1.values, facecolor = "rgb(25, 255, 25)");
        #     plots = cat(tbox, tenv0, tenv1; dims = 1)
        #     react!(pl, plots, pl.plot.layout)
        #     sleep(0.115)
        # end

        return true
    end # argyr_l_frame_modal_anim

end
using .mmodal1
mmodal1.test()

module mmassform11
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMCorotBeamModule
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
mass = FEMMCorotBeamModule.mass
geostiffness = FEMMCorotBeamModule.geostiffness
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using FinEtoolsFlexStructures.FEMMCorotBeamModule: local_frame_and_def!, local_mass_original!, local_stiffness!, natural_forces!, local_geometric_stiffness!, local_forces!
using FinEtoolsFlexStructures.FEMMCorotBeamModule: local_mass_CONSISTENT_WITH_ROTATION_INERTIA!, MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA

using Test

function test()
    # Parameters:
    E = 1000000.0;
    nu = 0.3;
    rho = 2700
    L1 = L0 = L =   30.0; # Length of the beam
    b =  0.5; # width
    h =  4.0; # height
    ecoords0 = fill(0.0, 2, 3); 
    ecoords1 = fill(0.0, 2, 3)
    edisp1 = fill(0.0, 2, 3); 
    evel1 = fill(0.0, 2, 6); 
    evel1f = fill(0.0, 2, 6)
    dofnums = zeros(FInt, 1, 12); 
    F0 = fill(0.0, 3, 3); 
    Ft = fill(0.0, 3, 3); 
    FtI = fill(0.0, 3, 3); 
    FtJ = fill(0.0, 3, 3)
    Te = fill(0.0, 12, 12)
    tempelmat1 = fill(0.0, 12, 12); 
    tempelmat2 = fill(0.0, 12, 12); 
    tempelmat3 = fill(0.0, 12, 12)
    elmat = fill(0.0, 12, 12);    
    elmatTe = fill(0.0, 12, 12);    
    elmato = fill(0.0, 12, 12)
    elvec = fill(0.0, 12);    
    elvecf = fill(0.0, 12)
    aN = fill(0.0, 6, 12)
    dN = fill(0.0, 6)
    DN = fill(0.0, 6, 6)
    PN = fill(0.0, 6)
    LF = fill(0.0, 12)
    R1I = fill(0.0, 3, 3);    
    R1J = fill(0.0, 3, 3);    
    OS = fill(0.0, 3, 3)

    # Reference frequencies
    reffs = [48.5475, 124.1839]

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [-1.0, 0.0, 0.0])

    # Select the number of elements per leg.
    n = 32;
    members = []
    push!(members, frame_member([0 0 0; 0 0 L], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000);

    dfes = delegateof(fes)
    A, I1, I2, I3, J, x1x2_vector = dfes.A, dfes.I1, dfes.I2, dfes.I3, dfes.J, dfes.x1x2_vector
    G = E / 2 / (1 + nu)
    i = 1

    A, I1, I2, I3, rho, L = 1.3, 35.1, 32.3, 16.5, 10000.0, 3.333
    MM = fill(0.0, 12, 12)
    local_mass_CONSISTENT_WITH_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    MMp = fill(0.0, 12, 12)
    local_mass_original!(MMp, A, I1, I2, I3, rho, L, MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA)
    @test norm(MM - MMp) / norm(MM) <= 1.0e-6
    true
end
end
using .mmassform11
mmassform11.test()

module mmassform10
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMCorotBeamModule
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
mass = FEMMCorotBeamModule.mass
geostiffness = FEMMCorotBeamModule.geostiffness
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using FinEtoolsFlexStructures.FEMMCorotBeamModule: local_frame_and_def!, local_mass_original!, local_stiffness!, natural_forces!, local_geometric_stiffness!, local_forces!
using FinEtoolsFlexStructures.FEMMCorotBeamModule: local_mass_CONSISTENT_WITH_ROTATION_INERTIA!, MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA, local_mass_CONSISTENT_NO_ROTATION_INERTIA!, MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA

using Test
function test()
    # Parameters:
    E = 1000000.0;
    nu = 0.3;
    rho = 2700
    L1 = L0 = L =   30.0; # Length of the beam
    b =  0.5; # width
    h =  4.0; # height
    ecoords0 = fill(0.0, 2, 3); 
    ecoords1 = fill(0.0, 2, 3)
    edisp1 = fill(0.0, 2, 3); 
    evel1 = fill(0.0, 2, 6); 
    evel1f = fill(0.0, 2, 6)
    dofnums = zeros(FInt, 1, 12); 
    F0 = fill(0.0, 3, 3); 
    Ft = fill(0.0, 3, 3); 
    FtI = fill(0.0, 3, 3); 
    FtJ = fill(0.0, 3, 3)
    Te = fill(0.0, 12, 12)
    tempelmat1 = fill(0.0, 12, 12); 
    tempelmat2 = fill(0.0, 12, 12); 
    tempelmat3 = fill(0.0, 12, 12)
    elmat = fill(0.0, 12, 12);    
    elmatTe = fill(0.0, 12, 12);    
    elmato = fill(0.0, 12, 12)
    elvec = fill(0.0, 12);    
    elvecf = fill(0.0, 12)
    aN = fill(0.0, 6, 12)
    dN = fill(0.0, 6)
    DN = fill(0.0, 6, 6)
    PN = fill(0.0, 6)
    LF = fill(0.0, 12)
    R1I = fill(0.0, 3, 3);    
    R1J = fill(0.0, 3, 3);    
    OS = fill(0.0, 3, 3)

    # Reference frequencies
    reffs = [48.5475, 124.1839]

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [-1.0, 0.0, 0.0])

    # Select the number of elements per leg.
    n = 32;
    members = []
    push!(members, frame_member([0 0 0; 0 0 L], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000);

    dfes = delegateof(fes)
    A, I1, I2, I3, J, x1x2_vector = dfes.A, dfes.I1, dfes.I2, dfes.I3, dfes.J, dfes.x1x2_vector
    G = E / 2 / (1 + nu)
    i = 1

    A, I1, I2, I3, rho, L = 1.3, 35.1, 32.3, 16.5, 10000.0, 3.333
    MM = fill(0.0, 12, 12)
    local_mass_CONSISTENT_NO_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    MMp = fill(0.0, 12, 12)
    local_mass_original!(MMp, A, I1, I2, I3, rho, L, MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA)
    @test norm(MM - MMp) / norm(MM) <= 1.0e-6
    true
end
end
using .mmassform10
mmassform10.test()

module mmassform12
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMCorotBeamModule
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
mass = FEMMCorotBeamModule.mass
geostiffness = FEMMCorotBeamModule.geostiffness
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using FinEtoolsFlexStructures.FEMMCorotBeamModule: local_frame_and_def!, local_mass_original!, local_stiffness!, natural_forces!, local_geometric_stiffness!, local_forces!
using FinEtoolsFlexStructures.FEMMCorotBeamModule: local_mass_CONSISTENT_WITH_ROTATION_INERTIA!, MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA, local_mass_LUMPED_DIAGONAL_NO_ROTATION_INERTIA!, MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA

using Test
function test()
    # Parameters:
    E = 1000000.0;
    nu = 0.3;
    rho = 2700
    L1 = L0 = L =   30.0; # Length of the beam
    b =  0.5; # width
    h =  4.0; # height
    ecoords0 = fill(0.0, 2, 3); 
    ecoords1 = fill(0.0, 2, 3)
    edisp1 = fill(0.0, 2, 3); 
    evel1 = fill(0.0, 2, 6); 
    evel1f = fill(0.0, 2, 6)
    dofnums = zeros(FInt, 1, 12); 
    F0 = fill(0.0, 3, 3); 
    Ft = fill(0.0, 3, 3); 
    FtI = fill(0.0, 3, 3); 
    FtJ = fill(0.0, 3, 3)
    Te = fill(0.0, 12, 12)
    tempelmat1 = fill(0.0, 12, 12); 
    tempelmat2 = fill(0.0, 12, 12); 
    tempelmat3 = fill(0.0, 12, 12)
    elmat = fill(0.0, 12, 12);    
    elmatTe = fill(0.0, 12, 12);    
    elmato = fill(0.0, 12, 12)
    elvec = fill(0.0, 12);    
    elvecf = fill(0.0, 12)
    aN = fill(0.0, 6, 12)
    dN = fill(0.0, 6)
    DN = fill(0.0, 6, 6)
    PN = fill(0.0, 6)
    LF = fill(0.0, 12)
    R1I = fill(0.0, 3, 3);    
    R1J = fill(0.0, 3, 3);    
    OS = fill(0.0, 3, 3)

    # Reference frequencies
    reffs = [48.5475, 124.1839]

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [-1.0, 0.0, 0.0])

    # Select the number of elements per leg.
    n = 32;
    members = []
    push!(members, frame_member([0 0 0; 0 0 L], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000);

    dfes = delegateof(fes)
    A, I1, I2, I3, J, x1x2_vector = dfes.A, dfes.I1, dfes.I2, dfes.I3, dfes.J, dfes.x1x2_vector
    G = E / 2 / (1 + nu)
    i = 1

    A, I1, I2, I3, rho, L = 1.3, 35.1, 32.3, 16.5, 10000.0, 3.333
    MM = fill(0.0, 12, 12)
    local_mass_LUMPED_DIAGONAL_NO_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    MMp = fill(0.0, 12, 12)
    local_mass_original!(MMp, A, I1, I2, I3, rho, L, MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA)
    @test norm(MM - MMp) / norm(MM) <= 1.0e-6
    true
end
end
using .mmassform12
mmassform12.test()

module mmassform13
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMCorotBeamModule
using FinEtoolsFlexStructures.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
mass = FEMMCorotBeamModule.mass
geostiffness = FEMMCorotBeamModule.geostiffness
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using FinEtoolsFlexStructures.FEMMCorotBeamModule: local_frame_and_def!, local_mass_original!, local_stiffness!, natural_forces!, local_geometric_stiffness!, local_forces!
using FinEtoolsFlexStructures.FEMMCorotBeamModule: local_mass_CONSISTENT_WITH_ROTATION_INERTIA!, MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA, local_mass_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA!, local_mass_original!, MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA

using Test
function test()
    # Parameters:
    E = 1000000.0;
    nu = 0.3;
    rho = 2700
    L1 = L0 = L =   30.0; # Length of the beam
    b =  0.5; # width
    h =  4.0; # height
    ecoords0 = fill(0.0, 2, 3); 
    ecoords1 = fill(0.0, 2, 3)
    edisp1 = fill(0.0, 2, 3); 
    evel1 = fill(0.0, 2, 6); 
    evel1f = fill(0.0, 2, 6)
    dofnums = zeros(FInt, 1, 12); 
    F0 = fill(0.0, 3, 3); 
    Ft = fill(0.0, 3, 3); 
    FtI = fill(0.0, 3, 3); 
    FtJ = fill(0.0, 3, 3)
    Te = fill(0.0, 12, 12)
    tempelmat1 = fill(0.0, 12, 12); 
    tempelmat2 = fill(0.0, 12, 12); 
    tempelmat3 = fill(0.0, 12, 12)
    elmat = fill(0.0, 12, 12);    
    elmatTe = fill(0.0, 12, 12);    
    elmato = fill(0.0, 12, 12)
    elvec = fill(0.0, 12);    
    elvecf = fill(0.0, 12)
    aN = fill(0.0, 6, 12)
    dN = fill(0.0, 6)
    DN = fill(0.0, 6, 6)
    PN = fill(0.0, 6)
    LF = fill(0.0, 12)
    R1I = fill(0.0, 3, 3);    
    R1J = fill(0.0, 3, 3);    
    OS = fill(0.0, 3, 3)

    # Reference frequencies
    reffs = [48.5475, 124.1839]

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [-1.0, 0.0, 0.0])

    # Select the number of elements per leg.
    n = 32;
    members = []
    push!(members, frame_member([0 0 0; 0 0 L], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000);

    dfes = delegateof(fes)
    A, I1, I2, I3, J, x1x2_vector = dfes.A, dfes.I1, dfes.I2, dfes.I3, dfes.J, dfes.x1x2_vector
    G = E / 2 / (1 + nu)
    i = 1

    A, I1, I2, I3, rho, L = 1.3, 35.1, 32.3, 16.5, 10000.0, 3.333
    MM = fill(0.0, 12, 12)
    local_mass_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    MMp = fill(0.0, 12, 12)
    local_mass_original!(MMp, A, I1, I2, I3, rho, L, MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA)
    @test norm(MM - MMp) / norm(MM) <= 1.0e-6
    true
end
end
using .mmassform13
mmassform13.test()

module mmodallin1
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMLinBeamModule
using FinEtoolsFlexStructures.FEMMLinBeamModule: FEMMLinBeam
stiffness = FEMMLinBeamModule.stiffness
mass = FEMMLinBeamModule.mass
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test

function test(mass_type = FEMMLinBeamModule.MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA)
        # Parameters:
        E=71240.0;#MPa
        nu=0.31;# Poisson ratio
        rho=5e-9;
        b=3.0; h=30.0; L=240.0; # cross-sectional dimensions and length of each leg in millimeters
        # Choose the mass formulation:
        scale = 0.4

        # Reference frequencies
        reffs = [11.2732, 30.5269]
        neigvs = 2;

        # Cross-sectional properties
        cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 1.0, 0.0])

        # Select the number of elements per leg.
        xyz =
        n = 4;
        members = []
        push!(members, frame_member([0 0 L; L 0 L], n, cs))
        push!(members, frame_member([L 0 L; L 0 0], n, cs))
        fens, fes = merge_members(members; tolerance = L / 10000);

        # Material properties
        material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

        # Construct the requisite fields, geometry and displacement
        # Initialize configuration variables
        geom0 = NodalField(fens.xyz)
        u0 = NodalField(zeros(size(fens.xyz,1), 3))
        Rfield0 = initial_Rfield(fens)
        dchi = NodalField(zeros(size(fens.xyz,1), 6))

        # Apply EBC's
        l1 = selectnode(fens; box = [0 0 0 0 L L], tolerance = L/10000)
        for i in [1,2,3,4,5,6]
            setebc!(dchi, l1, true, i)
        end
        applyebc!(dchi)
        numberdofs!(dchi);

        # Assemble the global discrete system
        femm = FEMMLinBeam(IntegDomain(fes, GaussRule(1, 2)), material)
        K = stiffness(femm, geom0, u0, Rfield0, dchi);
        M = mass(femm, geom0, u0, Rfield0, dchi; mass_type = mass_type);
        K_ff = matrix_blocked(K, nfreedofs(dchi), nfreedofs(dchi))[:ff]
        M_ff = matrix_blocked(M, nfreedofs(dchi), nfreedofs(dchi))[:ff]

        # Solve the eigenvalue problem
        d,v,nev,nconv = eigs(K_ff, M_ff; nev=2*neigvs, which=:SM)
        fs = real(sqrt.(complex(d)))/(2*pi)
        # println("Natural frequencies: $fs [Hz]")
        # println("Reference: $reffs [Hz]")

        @test norm(reffs - fs[1:length(reffs)]) ./ norm(reffs) <= 5.0e-2

        # Visualize vibration modes
        # mode = 2
        # tbox = plot_space_box([[-L/2 -L/2 0]; [L/2 L/2 1.1*L]])
        # tenv0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
        # plots = cat(tbox, tenv0; dims = 1)
        # pl = render(plots; title = "Mode $(mode)")
        # Rfield1 = deepcopy(Rfield0)
        # for xscale in scale.*sin.(collect(0:1:72).*(2*pi/21))
        #     scattersysvec!(dchi, xscale.*v[:, mode])
        #     Rfield1 = deepcopy(Rfield0)
        #     update_rotation_field!(Rfield1, dchi)
        #     tenv1 = plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield1.values, facecolor = "rgb(25, 255, 25)");
        #     plots = cat(tbox, tenv0, tenv1; dims = 1)
        #     react!(pl, plots, pl.plot.layout)
        #     sleep(0.115)
        # end

        return true
    end # argyr_l_frame_modal_anim
    for mass_type in [FEMMLinBeamModule.MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA,
        FEMMLinBeamModule.MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA,
        FEMMLinBeamModule.MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA,
        FEMMLinBeamModule.MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA]
        test(mass_type)
    end
end


