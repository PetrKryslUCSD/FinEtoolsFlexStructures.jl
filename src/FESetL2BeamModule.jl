module FESetL2BeamModule

using FinEtools
import FinEtools: cat, subset
using ..CrossSectionModule: AbstractCrossSectionType
using LinearAlgebra: norm, Transpose, mul!

mutable struct FESetL2Beam{CT} <: AbstractFESet1Manifold{2}
    crosssection::CT
    A::FFltVec
    I1::FFltVec
    I2::FFltVec
    I3::FFltVec
    J::FFltVec
    A2s::FFltVec
    A3s::FFltVec
    x1x2_vector::Vector{FFltVec}
    dimensions::Vector{FFltVec}
end

function FESetL2Beam(N::IT, crosssection::CT) where {IT<:Integer, CT}
    par = crosssection.parameters(0.0)
    N = size(conn, 1)
    _A = fill(par.A, N)
    _I1 = fill(par.I1, N)
    _I2 = fill(par.I2, N)
    _I3 = fill(par.I3, N)
    _J = fill(par.J, N)
    _A2s = fill(par.A2s, N) 
    _A3s = fill(par.A3s, N)
    _x1x2_vector = [par.x1x2_vector for i in 1:N]
    _dimensions = [par.dimensions for i in 1:N]
    return  FESetL2Beam(crosssection, _A, _I1, _I2, _I3, _J, _A2s, _A3s, _x1x2_vector, _dimensions)
end

"""
    cat(self::T,  other::T) where {T<:FESetL2Beam}

Concatenate two sets of beam elements.
"""
function cat(self::T,  other::T) where {T<:FESetL2Beam}
    self.crosssection === other.crosssection  || error("Cannot concatenate sets with distinct cross-sections")
    result = deepcopy(self)
    result.A = vcat(self.A, other.A);
    result.I1 = vcat(self.I1, other.I1);
    result.I2 = vcat(self.I2, other.I2);
    result.I3 = vcat(self.I3, other.I3);
    result.J = vcat(self.J, other.J);
    result.A2s = vcat(self.A2s, other.A2s);
    result.A3s = vcat(self.A3s, other.A3s);
    result.x1x2_vector = vcat(self.x1x2_vector, other.x1x2_vector);
    result.dimensions = vcat(self.dimensions, other.dimensions);
    return result
end

"""
    subset(self::T, L::FIntVec) where {T<:FESetL2Beam}

Subset of a beam-element set
"""
function subset(self::T, L::FIntVec) where {T<:FESetL2Beam}
    result = deepcopy(self)
    result.A = self.A[L]
    result.I1 = self.I1[L]
    result.I2 = self.I2[L]
    result.I3 = self.I3[L]
    result.J = self.J[L] 
    result.A2s = self.A2s[L] 
    result.A3s = self.A3s[L] 
    result.x1x2_vector = self.x1x2_vector[L]
    result.dimensions = self.dimensions[L]
    return  result
end

function initial_local_frame!(F0, x0, x1x2_vector)
    # This is the element frame in the configuration t=0
    F0[:,1] = (x0[2,:]-x0[1,:]); 
    L0 = norm(@view F0[:,1]);
    F0[:,1] /= L0;
    #     F0(:,3)=skewmat(F0(:,1))*x1x2_vector;
    F0[:,3] .= (-F0[3,1]*x1x2_vector[2]+F0[2,1]*x1x2_vector[3],
                 F0[3,1]*x1x2_vector[1]-F0[1,1]*x1x2_vector[3],
                -F0[2,1]*x1x2_vector[1]+F0[1,1]*x1x2_vector[2]);
    F0[:,3] /= norm(@view F0[:,3]);
    #     F0(:,2)=skewmat(F0(:,3))*F0(:,1);
    F0[:,2] .= (-F0[3,3]*F0[2,1]+F0[2,3]*F0[3,1],
                 F0[3,3]*F0[1,1]-F0[1,3]*F0[3,1],
                -F0[2,3]*F0[1,1]+F0[1,3]*F0[2,1]);
    return L0, F0
end

end # module
