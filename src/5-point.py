# -*- coding: utf-8 -*-
"""
Spyder Editor


struct Quad5Rule   <: AbstractIntegRule
    dim::Int
    npts::Int
    param_coords::Array{Float64,2}
    weights::Array{Float64,2}
end

function Quad5Rule()
    param_coords = [0.0 0.0; -1.0 -1.0; 1.0 -1.0; 1.0 1.0; -1.0 1.0]
    # a = 7 / 108 # 0.065
    a = 1/36
    weights = 4 .* [1 - 4*a, a, a, a, a]'
    Quad5Rule(2, 5, param_coords, weights)
end

"""

from sympy import *

var('r, s')
N = Matrix([(r-1)*(s-1)/4, 
            (r+1)*(s-1)/-4, 
            (r+1)*(s+1)/4,
            (r-1)*(s+1)/-4]).reshape(4,1)

# 5-point
# var('a')
a = 7/108
pc = Matrix([[0, 0], [-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
weights = 4 * Matrix([1 - 4*a, a, a, a, a])

# Gauss
pc = 1/sqrt(3) * Matrix([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
weights = Matrix([1, 1, 1, 1])


integrand = N[1] * N[0]
I = integrate(integrand, (r, -1, 1), (s, -1, 1))

Ih = 0.0
for q in range(pc.shape[0]):
    Ih += integrand.subs(r, pc[q, 0]).subs(s, pc[q, 1]) * weights[q]
    
# print(I, ' = ', Ih)
print(float(I - Ih))