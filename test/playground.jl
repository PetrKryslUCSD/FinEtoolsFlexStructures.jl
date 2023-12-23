
module mlayup10
using LinearAlgebra: norm, Transpose, mul!, I
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CompositeLayupModule
using Test
function test()
    CM = CompositeLayupModule
    E = 200e3 * phun("MPa")
    nu = 0.3
    rho = 8000 * phun("KG/M^3")
    thickness = 0.005 * phun("m")

    mater = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
    CM = CompositeLayupModule
    A = fill(0.0, 3, 3)
    B = fill(0.0, 3, 3)
    C = fill(0.0, 3, 3)
    Atrue = E / (1 - nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2] * thickness

    for angle in [0.0 47.0 90.0 129.0 180.0]
        ply1 = CM.Ply("ply1", mater, thickness, angle)
        cl = CM.CompositeLayup("sample", [ply1], CM.cartesian_csys((1, 2, 3)))
        @test cl.plies[1].angle == angle
        A, B, C = CM.laminate_stiffnesses!(cl, A, B, C)
        @test norm(A - Atrue) < 1.0e-15 * norm(Atrue)
    end
    true
end
end
using .mlayup10
mlayup10.test()
