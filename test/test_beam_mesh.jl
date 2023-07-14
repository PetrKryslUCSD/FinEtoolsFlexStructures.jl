module m2
using FinEtools
using FinEtoolsFlexStructures.FESetL2BeamModule: FESetL2Beam
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionCircle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member
using Test
function test()
    L = 42
    xyz = [0 0 0;
    0 L/4 L*1/4;
    L/4 L/4 L*2/4;
    L/4 0 L*3/4;
    0 0 L]
    nL = 20

    cs = CrossSectionCircle(s -> 5.9910, s -> [0.0, 0.0, 1.0])
    fens, fes = frame_member(xyz, nL, cs)
    @test count(fes) == 20
    true
end
end
using .m2
m2.test()
