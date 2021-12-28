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
    ply1 = cm.CompositeLayupPly("ply1", mater, thickness, angle)
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
    ply1 = CM.CompositeLayupPly("ply1", mater, thickness, 45.0)
    ply2 = CM.CompositeLayupPly("ply2", mater, thickness, -45.0)
    cl = CM.CompositeLayup("sample", [ply1, ply2])
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
        CM.plane_stress_T_matrix!(DeforModelRedModule.DeforModelRed2DStress, Tm, angle)
        Tinvm = fill(0.0, 3, 3)
        CM.plane_stress_Tinv_matrix!(DeforModelRedModule.DeforModelRed2DStress, Tinvm, angle)
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
        CM.plane_stress_T_matrix!(DeforModelRedModule.DeforModelRed2DStress, Tm, angle)
        CM.plane_stress_T_matrix_eng!(DeforModelRedModule.DeforModelRed2DStress, Tme, angle)
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
        CM.plane_stress_Tinv_matrix!(DeforModelRedModule.DeforModelRed2DStress, Tm, angle)
        CM.plane_stress_Tinv_matrix_eng!(DeforModelRedModule.DeforModelRed2DStress, Tme, angle)
        # @show R * Tm / R, Tme
        @test norm(R * Tm / R - Tme) < mtol
        # Tinvm = fill(0.0, 3, 3)
        # CM.plane_stress_Tinv_matrix!(DeforModelRedModule.DeforModelRed2DStress, Tinvm, angle)
        # @test norm(Tm*Tinvm - I) < mtol
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
    # for dangle in [5 34 68 73 535 0 90 45 -135. -87]
        angle = dangle/180*pi
        CM.plane_stress_T_matrix_eng!(DeforModelRedModule.DeforModelRed2DStress, Tm, angle)
        CM.plane_stress_Tinv_matrix_eng!(DeforModelRedModule.DeforModelRed2DStress, Tinvm, angle)
        @test norm(Tm*Tinvm - I) < mtol
    end
    true
end
end
using .mlayup6
mlayup6.test()
