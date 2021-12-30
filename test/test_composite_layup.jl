module mcsys_0
using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    XYZ = [0.0 0.0 0.0]
    J0 =  [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    labl = 1
    CM = CompositeLayupModule
    cs = CM.cartesian_csys((1, 2, 3))  
    updatecsmat!(cs, reshape(XYZ, 1, 3), J0, labl);
    @test norm(cs.csmat - 
        [1.0 0.0 0.0; 
        0.0 1.0 0.0; 
        0.0 0.0 1.0]) < 1.0e-15
    cs = CM.cartesian_csys((2, -1, 3))  
    updatecsmat!(cs, reshape(XYZ, 1, 3), J0, labl);
    @test norm(cs.csmat - 
            [0.0 -1.0 0.0; 
            1.0 0.0 0.0; 
            0.0 0.0 1.0]) < 1.0e-15
    cs = CM.cartesian_csys((2, 1, -3))  
    updatecsmat!(cs, reshape(XYZ, 1, 3), J0, labl);
    @test norm(cs.csmat - 
            [0.0 1.0 0.0; 
            1.0 0.0 0.0; 
            0.0 0.0 -1.0]) < 1.0e-15
    true
end
end
using .mcsys_0
mcsys_0.test()

module mlayup1
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    E = 200e3*phun("MPa")
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    thickness = 0.005*phun("m");
    angle = 45.0
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    cm = CompositeLayupModule
    ply1 = cm.Ply("ply1", mater, thickness, angle)
    @test ply1.angle == 45
    true
end
end
using .mlayup1
mlayup1.test()

module mlayup2
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using Test

function test()
    E = 200e3*phun("MPa")
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    thickness = 0.005*phun("m");
    angle = 45.0
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    CM = CompositeLayupModule
    ply1 = CM.Ply("ply1", mater, thickness, 45.0)
    ply2 = CM.Ply("ply2", mater, thickness, -45.0)
    cl = CM.CompositeLayup("sample", [ply1, ply2], CM.cartesian_csys((1, 2, 3))  )
    @test cl.plies[2].angle == -45
    true
end
end
using .mlayup2
mlayup2.test()


module mlayup3
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    mtol = 1.0e-15
    for dangle in [5 34 68 73 535 0 90 45 -135. -87]
        angle = dangle/180*pi
        Tm = fill(0.0, 3, 3)
        CM.plane_stress_T_matrix!(Tm, angle)
        Tinvm = fill(0.0, 3, 3)
        CM.plane_stress_Tinv_matrix!(Tinvm, angle)
        @test norm(Tm*Tinvm - I) < mtol
    end
    true
end
end
using .mlayup3
mlayup3.test()


module mlayup4
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    mtol = 1.0e-15
    R = [1 0 0; 0 1 0; 0 0 2]
    Tm = fill(0.0, 3, 3)
    Tme = fill(0.0, 3, 3)
    for dangle in [5 34 68 73 535 0 90 45 -135. -87]
        angle = dangle/180*pi
        CM.plane_stress_T_matrix!(Tm, angle)
        CM.plane_stress_T_matrix_eng!(Tme, angle)
        # @show R * Tm / R, Tme
        @test norm(R * Tm / R - Tme) < mtol
        # Tinvm = fill(0.0, 3, 3)
        # CM.plane_stress_Tinv_matrix!(DeforModelRedModule.DeforModelRed2DStress, Tinvm, angle)
        # @test norm(Tm*Tinvm - I) < mtol
    end
    true
end
end
using .mlayup4
mlayup4.test()


module mlayup5
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    mtol = 1.0e-15
    R = [1 0 0; 0 1 0; 0 0 2]
    Tm = fill(0.0, 3, 3)
    Tme = fill(0.0, 3, 3)
    for dangle in [5 34 68 73 535 0 90 45 -135. -87]
        angle = dangle/180*pi
        CM.plane_stress_Tinv_matrix!(Tm, angle)
        CM.plane_stress_Tinv_matrix_eng!(Tme, angle)
        @test norm(R * Tm / R - Tme) < mtol
    end
    true
end
end
using .mlayup5
mlayup5.test()

module mlayup6
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    mtol = 1.0e-15
    R = [1 0 0; 0 1 0; 0 0 2]
    Tm = fill(0.0, 3, 3)
    Tinvm = fill(0.0, 3, 3)
    for dangle in [-8.]
        angle = dangle/180*pi
        CM.plane_stress_T_matrix_eng!(Tm, angle)
        CM.plane_stress_Tinv_matrix_eng!(Tinvm, angle)
        @test norm(Tm*Tinvm - I) < mtol
    end
    true
end
end
using .mlayup6
mlayup6.test()

module mlayup7
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    mtol = 1.0e-15
    R = [1 0 0; 0 1 0; 0 0 2]
    Tm = fill(0.0, 3, 3)
    Tinvm = fill(0.0, 3, 3)
    for dangle in [5 34 68 73 535 0 90 45 -135. -87]
        angle = dangle/180*pi
        CM.plane_stress_T_matrix_eng!(Tm, -angle)
        CM.plane_stress_Tinv_matrix_eng!(Tinvm, angle)
        @test norm(Tm - Tinvm) < mtol
    end
    true
end
end
using .mlayup7
mlayup7.test()

module mlayup8
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    mtol = 1.0e-15
    R = [1 0 0; 0 1 0; 0 0 2]
    Tm = fill(0.0, 3, 3)
    Tinvm = fill(0.0, 3, 3)
    for dangle in [-55] # see Barbero, Introduction to composite materials, Example 5.3
        angle = dangle/180*pi
        CM.plane_stress_T_matrix!(Tm, angle)
        @test norm(Tm - [0.32898992833716556 0.6710100716628343 -0.9396926207859083; 0.6710100716628343 0.32898992833716556 0.9396926207859083; 0.46984631039295416 -0.46984631039295416 -0.34202014332566877]   ) < mtol
    end
    true
end
end
using .mlayup8
mlayup8.test()

module mlayup9
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    mtol = 1.0e-15
    R = [1 0 0; 0 1 0; 0 0 2]
    Tme = fill(0.0, 3, 3)
    Tm = fill(0.0, 3, 3)
    Tinvm = fill(0.0, 3, 3)
    for dangle in [5 34 68 73 535 0 90 45 -135. -87]
        angle = dangle/180*pi
        CM.plane_stress_T_matrix_eng!(Tme, angle)
        CM.plane_stress_T_matrix!(Tm, angle)
        # Verify the identity above equation 5.54 in Barbero
        @test norm(inv(Tm') - Tme) < mtol
    end
    true
end
end
using .mlayup9
mlayup9.test()


module mlayup10
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    E = 200e3*phun("MPa")
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    thickness = 0.005*phun("m");
    
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    CM = CompositeLayupModule
    A = fill(0.0, 3, 3)
    B = fill(0.0, 3, 3)
    C = fill(0.0, 3, 3)
    Atrue = E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2] * thickness

    for angle in [0.0 47.0 90.0 129.0 180.0]
        ply1 = CM.Ply("ply1", mater, thickness, angle)
        cl = CM.CompositeLayup("sample", [ply1, ], CM.cartesian_csys((1, 2, 3)))
        @test cl.plies[1].angle == angle
        A, B, C = CM.laminate_stiffnesses!(cl, A, B, C)
        @test norm(A - Atrue) < 1.0e-15 * norm(Atrue)
    end
    true
end
end
using .mlayup10
mlayup10.test()


module mlayup11
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    E = 200e3*phun("MPa")
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    thickness = 0.005*phun("m");
    angle = 45.0
    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    CM = CompositeLayupModule
    # From Barbero's Introduction ... book Example 5.6
    ply1 = CM.Ply("ply1", mater, 1.0, -55.0)
    cl = CM.CompositeLayup("sample", [ply1, ], CM.cartesian_csys((1, 2, 3)))
    cl.plies[1]._Dps .= [20.874 3.260 0; 3.260 11.898 0; 0 0 3.789]
    A = fill(0.0, 3, 3)
    B = fill(0.0, 3, 3)
    C = fill(0.0, 3, 3)
    A, B, C = CM.laminate_stiffnesses!(cl, A, B, C)
    Atrue = [12.401509954148374 5.709503642606021 -1.217123826408347; 5.709503642606021 15.471482760639574 -3.000216655678807; -1.2171238264083473 -3.000216655678808 6.2385036426060205]
    @test norm(A - Atrue) < 1.0e-15 * norm(Atrue)
    true
end
end
using .mlayup11
mlayup11.test()

module mlayup12
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    # From Barbero's Introduction ... book Example 5.7
    E1 = 67192*phun("MPa")
    E2 = 12139*phun("MPa")
    G12 = 7638*phun("MPa")
    nu12 = 0.365;
    rho= 8000*phun("KG/M^3");
    thickness = 1.0*phun("m");
    mater = MatDeforElastOrtho(DeforModelRed3D, rho, E1, E2, E2, nu12, nu12, nu12, G12, G12, G12, 0.0, 0.0, 0.0)
    CM = CompositeLayupModule
    A = fill(0.0, 3, 3)
    B = fill(0.0, 3, 3)
    C = fill(0.0, 3, 3)
    
    for angle in [45]
        ply1 = CM.Ply("ply1", mater, thickness, angle)
        # @show ply1._Dps./phun("MPa")
        @test norm(ply1._Dps./phun("MPa") - [68849.10243290635 4540.006665496835 0.0; 
                                                4540.006665496834 12438.374426018723 0.0; 
                                                0.0 0.0 7637.999999999999]) < 1.0e-15 * E2
        cl = CM.CompositeLayup("sample", [ply1, ], CM.cartesian_csys((1, 2, 3)))
        @test cl.plies[1].angle == angle
        A, B, C = CM.laminate_stiffnesses!(cl, A, B, C)
        # @show A./phun("MPa")
        @test norm(A ./ phun("MPa") - [30229.872547479692 14953.872547479687 14102.68200172191;
                                       14953.872547479687 30229.872547479685 14102.682001721907; 
                                       14102.682001721909 14102.682001721909 18051.865881982852]) < 1.0e-15 * E2
    end
    true
end
end
using .mlayup12
mlayup12.test()

module mlayup13
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    # From Barbero's Introduction ... book Example 5.7
    E1 = 67192*phun("MPa")
    E2 = 12139*phun("MPa")
    G12 = 7638*phun("MPa")
    nu12 = 0.365;
    rho= 8000*phun("KG/M^3");
    thickness = 1.0*phun("m");
    CM = CompositeLayupModule
    mater = CM.lamina_material(E1, E2, nu12, G12, G12, G12)
    A = fill(0.0, 3, 3)
    B = fill(0.0, 3, 3)
    C = fill(0.0, 3, 3)
    
    for angle in [45]
        ply1 = CM.Ply("ply1", mater, thickness, angle)
        # @show ply1._Dps./phun("MPa")
        @test norm(ply1._Dps./phun("MPa") - [68849.10243290635 4540.006665496835 0.0; 
                                                4540.006665496834 12438.374426018723 0.0; 
                                                0.0 0.0 7637.999999999999]) < 1.0e-15 * E2
        cl = CM.CompositeLayup("sample", [ply1, ])
        @test cl.plies[1].angle == angle
        A, B, C = CM.laminate_stiffnesses!(cl, A, B, C)
        # @show A./phun("MPa")
        @test norm(A ./ phun("MPa") - [30229.872547479692 14953.872547479687 14102.68200172191;
                                       14953.872547479687 30229.872547479685 14102.682001721907; 
                                       14102.682001721909 14102.682001721909 18051.865881982852]) < 1.0e-15 * E2
    end
    true
end
end
using .mlayup13
mlayup13.test()

module mlayup14
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    # From Barbero's Introduction ... book Example 5.7
    E1 = 67192*phun("MPa")
    E2 = 12139*phun("MPa")
    G12 = 7638*phun("MPa")
    nu12 = 0.365;
    rho= 8000*phun("KG/M^3");
    thickness = 0.5*phun("m");
    CM = CompositeLayupModule
    mater = CM.lamina_material(E1, E2, nu12, G12, G12, G12)
    A = fill(0.0, 3, 3)
    B = fill(0.0, 3, 3)
    C = fill(0.0, 3, 3)
    # "balanced" fabric
    ply1 = CM.Ply("ply1", mater, thickness, 45)
    ply2 = CM.Ply("ply2", mater, thickness, -45)
    @test norm(ply1._Dps./phun("MPa") - 
        [68849.10243290635 4540.006665496835 0.0; 
        4540.006665496834 12438.374426018723 0.0; 
        0.0 0.0 7637.999999999999]) < 1.0e-15 * E2
    @test norm(ply2._Dps./phun("MPa") - 
        [68849.10243290635 4540.006665496835 0.0; 
        4540.006665496834 12438.374426018723 0.0; 
        0.0 0.0 7637.999999999999]) < 1.0e-15 * E2
    cl = CM.CompositeLayup("sample", [ply1, ply2])
    A, B, C = CM.laminate_stiffnesses!(cl, A, B, C)
        # @show A./phun("MPa")
    @test norm(A ./ phun("MPa") - 
        [30229.872547479692 14953.872547479687 0.0; 
        14953.872547479687 30229.872547479685 0.0; 
        0.0 0.0 18051.865881982852] ) < 1.0e-15 * E2
    
    true
end
end
using .mlayup14
mlayup14.test()

module mlayup15
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    # From Barbero's Introduction ... book Example 5.7
    E1 = 67192*phun("MPa")
    E2 = 12139*phun("MPa")
    G12 = 7638*phun("MPa")
    nu12 = 0.365;
    rho= 8000*phun("KG/M^3");
    thickness = 0.001*phun("m");
    CM = CompositeLayupModule
    mater = CM.lamina_material(E1, E2, nu12, G12, G12, G12)
    A = fill(0.0, 3, 3)
    B = fill(0.0, 3, 3)
    C = fill(0.0, 3, 3)
    # "cross-ply" fabric: no coupling of extension and bending
    ply1 = CM.Ply("ply1", mater, thickness, 0)
    ply2 = CM.Ply("ply2", mater, thickness, 90)
    ply3 = CM.Ply("ply3", mater, thickness, 0)
    ply4 = CM.Ply("ply4", mater, thickness, 90)
    ply5 = CM.Ply("ply5", mater, thickness, 0)
    cl = CM.CompositeLayup("sample", [ply1, ply2, ply3, ply4, ply5])
    A, B, C = CM.laminate_stiffnesses!(cl, A, B, C)
        # @show A./phun("MPa")
    @test norm(B) < 1.0e-15 * E2
    
    true
end
end
using .mlayup15
mlayup15.test()

module mlayup16
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    # From Barbero's Introduction ... book Example 5.1
    E1 = 19981*phun("MPa")
    E2 = 11389*phun("MPa")
    G12 = 3789*phun("MPa")
    G13 = G12
    nu12 = 0.274;
    nu23 = 0.385
    G23 = E2 / 2 / (1 + nu23)
    rho= 8000*phun("KG/M^3");
    thickness = 0.001*phun("m");
    CM = CompositeLayupModule
    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    A = fill(0.0, 3, 3)
    B = fill(0.0, 3, 3)
    C = fill(0.0, 3, 3)
    H = fill(0.0, 2, 2)
    # "cross-ply" fabric: no coupling of extension and bending
    ply1 = CM.Ply("ply1", mater, thickness, -55)
    cl = CM.CompositeLayup("sample", [ply1, ])
    A, B, C = CM.laminate_stiffnesses!(cl, A, B, C)
        # @show A./phun("GPa")/thickness
    @test norm(A./phun("GPa")/thickness - [
        12.401649598137476 5.70964888928013 -1.2171315230336193; 
        5.70964888928013 15.47166067503043 -3.000261532585858; 
        -1.2171315230336193 -3.0002615325858573 6.238554718598339]) < 1.0e-15 * E2
    H = CM.laminate_transverse_stiffness!(cl, H)
    @show H./phun("GPa")/thickness
    @test norm(H./phun("GPa")/thickness - [
        3.3378632276560145 -0.1262916916205752; 
        -0.12629169162057502 3.245930394485983]) < 1.0e-15 * E2
    true
end
end
using .mlayup16
mlayup16.test()


module mlayup17
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using Test

function test()
    CM = CompositeLayupModule
    # From Barbero's Finite Element Analysis using Abaqus ... book Example 3.1
    ax = ay = 2000*phun("mm")
    # ASFD/9310
    E1 = 133860*phun("MPa")
    E2 = 7706*phun("MPa")
    G12 = 4306*phun("MPa")
    G13 = G12
    nu12 = 0.301;
    nu23 = 0.396
    G23 = 2760*phun("MPa")
    npairs = 2
    thickness = 10*phun("mm");
    mcsys = CM.cartesian_csys((1, 2, 3))
    CM = CompositeLayupModule
    mater = CM.lamina_material(E1, E2, nu12, G12, G13, G23)
    plies = CM.Ply[]
    for p in npairs
        push!(plies, CM.Ply("ply_0_$p", mater, thickness, 0))
        push!(plies, CM.Ply("ply_90_$p", mater, thickness, 90))
    end
    cl = CM.CompositeLayup("example_3.1", plies, mcsys)
    
    true
end
end
using .mlayup17
mlayup17.test()