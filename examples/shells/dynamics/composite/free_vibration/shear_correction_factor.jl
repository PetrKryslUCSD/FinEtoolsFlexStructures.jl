layup_thickness  = 1.0
f(zs, ze) = 5 / 4 * (ze - zs - 4 / 3 * (ze^3 - zs^3) / layup_thickness^2)

@show f(-1/2, 1/2)
@show f(-1/2, 0) + f(0, 1/2)