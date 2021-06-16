module mcrosssection1
using FinEtools
using FinEtoolsFlexStructures.CrossSectionModule
using Test
function test()
    cs = CrossSectionModule.CrossSectionCircle(s -> 5.9910, s -> [0.0, 0.0, 1.0])
    par = cs.parameters(0.0)
    for (c, r) in zip((par.A, par.J, par.I1, par.I2, par.I3), (112.75829799164978, 2023.5649824692157, 2023.5649824692157, 1011.7824912346078, 1011.7824912346078))
        @test c ≈ r
    end 
    true
end
end
using .mcrosssection1
mcrosssection1.test()

module mcrosssection2
using FinEtools
using FinEtoolsFlexStructures.CrossSectionModule
using Test
function test()
    cs = CrossSectionModule.CrossSectionRectangle(s -> 42.0, s -> 4.2, s -> [0.0, 0.0, 1.0])
    par = cs.parameters(0.0)
    for (c, r) in zip((par.A, par.J, par.I1, par.I2, par.I3), (176.4, 970.849152, 26190.108000000004, 259.30800000000005, 25930.800000000003))
        @test c ≈ r
    end 
    true
end
end
using .mcrosssection2
mcrosssection2.test()

module mcrosssection3
using FinEtools
using FinEtoolsFlexStructures.CrossSectionModule
using Test
function test()
    cs = CrossSectionModule.CrossSectionRectangle(s -> 1.3 * 4.2, s -> 4.2, s -> [0.0, 0.0, 1.0])
    par = cs.parameters(0.0)
    for (c, r) in zip((par.A, par.J, par.I1, par.I2, par.I3), (22.932000000000006, 71.66370760423138, 90.68000760000004, 33.71004000000001, 56.969967600000025))
         @test c ≈ r
    end 
 true
end
end
using .mcrosssection3
mcrosssection3.test()

module mcrosssection4
using FinEtools
using FinEtoolsFlexStructures.CrossSectionModule
using Test
function test()
    R = 0.5
    cs = CrossSectionModule.CrossSectionCircle(s -> (1/2+2*s)*R, s -> [0.0, 0.0, 1.0])
    par = cs.parameters(0.0)
    for (c, r) in zip((par.A, par.J, par.I1, par.I2, par.I3), (0.19634954084936207, 0.006135923151542565, 0.006135923151542565, 0.0030679615757712823, 0.0030679615757712823))
        @test c ≈ r
    end 
    true
end
end
using .mcrosssection4
mcrosssection4.test()