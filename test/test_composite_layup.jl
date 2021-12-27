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
    @show ply1
    true
end
end
using .mlayup1
mlayup1.test()

