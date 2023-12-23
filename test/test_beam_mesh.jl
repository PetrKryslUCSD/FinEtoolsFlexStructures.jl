module mframe1
using FinEtools
using FinEtoolsFlexStructures.FESetL2BeamModule: FESetL2Beam
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionCircle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member
using Test
function test()
    L = 42
    xyz = [
        0 0 0
        0 L/4 L*1/4
        L/4 L/4 L*2/4
        L/4 0 L*3/4
        0 0 L
    ]
    nL = 20

    cs = CrossSectionCircle(s -> 5.9910, s -> [0.0, 0.0, 1.0])
    fens, fes = frame_member(xyz, nL, cs)
    @test count(fes) == 20
    true
end
test()
end

module mframe2
using FinEtools
using FinEtoolsFlexStructures.FESetL2BeamModule: FESetL2Beam
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionCircle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member
using Test
function test()
    L = 42
    xyz = [
        0 0 0
        0 L/4 L*1/4
        L/4 L/4 L*2/4
        L/4 0 L*3/4
        0 0 L
    ]
    nL = 20

    cs = CrossSectionCircle(s -> 5.9910, s -> [0.0, 0.0, 1.0])
    fens, fes = frame_member(xyz, nL, cs)
    @test count(fes) == 20

    cs = CrossSectionCircle(s -> 5.9910, s -> [0.0, 0.0, 1.0])
    fens, fes = frame_member(xyz, nL, cs)
    list = 1:2
    fes = subset(fes, list)
    accepttodelegate(fes, subset(fes.delegateof, list))
    @test count(fes) == 2
    # @show fes
    true
end
test()
end

module mframe3
using FinEtools
using FinEtoolsFlexStructures.FESetL2BeamModule: FESetL2Beam
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionCircle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member
using Test
function test()
    L = 42
    xyz = [
        0 0 0
        0 L/4 L*1/4
        L/4 L/4 L*2/4
        L/4 0 L*3/4
        0 0 L
    ]
    nL = 20

    cs = CrossSectionCircle(s -> 5.9910, s -> [0.0, 0.0, 1.0])
    fes = FESetL2Beam(nL, cs)
    # @show fes
    true
end
test()
end
