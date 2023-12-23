module FEMMRITBeamModule

using LinearAlgebra: norm, Transpose, mul!
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.MatrixUtilityModule: add_n1n2t!
using FinEtools.IntegDomainModule: IntegDomain
import FinEtoolsDeforLinear.MatDeforElastIsoModule: MatDeforElastIso
using ..FESetL2BeamModule: FESetL2Beam, initial_local_frame!

"""
    FEMMRITBeam{S<:FESetL2, F<:Function} <: AbstractFEMM

Type for linear reduced-integration beam finite element modeling machine.

Only linear kinematics is implemented at the moment, and only linear basis
functions are available (i.e. it is a two-node element). The beam stiffness is
shear-flexible(Timoshenko). The local beam stiffness is expressed analytically:
the one-point numerical integration is hardwired.
"""
mutable struct FEMMRITBeam{S<:FESetL2,F<:Function} <: AbstractFEMM
    integdomain::IntegDomain{S,F} # integration domain data
    material::MatDeforElastIso # material object
    # The attributes below are buffers used in various operations.
    _ecoords0::FFltMat
    _dofnums::FIntMat
    _F0::FFltMat
    _Te::FFltMat
    _elmat::FFltMat
    _elvec::FFltVec
    _elvecf::FFltVec
    _b_phi1::FFltVec
    _b_phi2::FFltVec
    _b_phi3::FFltVec
    _b_w1::FFltVec
    _b_w2::FFltVec
    _b_w3::FFltVec
    _i_phi2::FFltVec
    _i_phi3::FFltVec
    _LF::FFltVec
    _RI::FFltMat
    _RJ::FFltMat
    _OS::FFltMat
end

"""
    FEMMRITBeam(integdomain::IntegDomain{S, F}, material::MatDeforElastIso) where {S<:FESetL2, F<:Function}

Constructor.

Supply integration domain (mesh) and the material.
"""
function FEMMRITBeam(
    integdomain::IntegDomain{S,F},
    material::MatDeforElastIso,
) where {S<:FESetL2,F<:Function}
    typeof(delegateof(integdomain.fes)) <: FESetL2Beam ||
        error("Expected to delegate to FESetL2Beam")
    _ecoords0 = fill(0.0, 2, 3)
    _dofnums = zeros(FInt, 1, 12)
    _F0 = fill(0.0, 3, 3)
    _Te = fill(0.0, 12, 12)
    _elmat = fill(0.0, 12, 12)
    _elvec = fill(0.0, 12)
    _elvecf = fill(0.0, 12)
    _b_phi1 = fill(0.0, 12)
    _b_phi2 = fill(0.0, 12)
    _b_phi3 = fill(0.0, 12)
    _b_w1 = fill(0.0, 12)
    _b_w2 = fill(0.0, 12)
    _b_w3 = fill(0.0, 12)
    _i_phi2 = fill(0.0, 12)
    _i_phi3 = fill(0.0, 12)
    _LF = fill(0.0, 12)
    _RI = fill(0.0, 3, 3)
    _RJ = fill(0.0, 3, 3)
    _OS = fill(0.0, 3, 3)
    A, I1, I2, I3, J, A2s, A3s, x1x2_vector, dimensions = properties(integdomain.fes)
    any(x -> isinf(x), A2s) && error("The shear areas must be finite numbers")
    any(x -> isinf(x), A3s) && error("The shear areas must be finite numbers")

    return FEMMRITBeam(
        integdomain,
        material,
        _ecoords0,
        _dofnums,
        _F0,
        _Te,
        _elmat,
        _elvec,
        _elvecf,
        _b_phi1,
        _b_phi2,
        _b_phi3,
        _b_w1,
        _b_w2,
        _b_w3,
        _i_phi2,
        _i_phi3,
        _LF,
        _RI,
        _RJ,
        _OS,
    )
end

function properties(fes)
    d = delegateof(fes)
    # A, I1, I2, I3, J, A2s, A3s, x1x2_vector, dimensions = properties(fes)
    d.A, d.I1, d.I2, d.I3, d.J, d.A2s, d.A3s, d.x1x2_vector, d.dimensions
end

function b_phi1!(b, h, F0)
    d = (F0[1, 1], F0[2, 1], F0[3, 1])
    b .= 0.0
    b[4:6] .= -1 / h .* d
    b[10:12] .= 1 / h .* d
    return b
end

function b_phi2!(b, h, F0)
    d = (F0[1, 2], F0[2, 2], F0[3, 2])
    b .= 0.0
    b[4:6] .= -1 / h .* d
    b[10:12] .= 1 / h .* d
    return b
end

function b_phi3!(b, h, F0)
    d = (F0[1, 3], F0[2, 3], F0[3, 3])
    b .= 0.0
    b[4:6] .= -1 / h .* d
    b[10:12] .= 1 / h .* d
    return b
end

function b_w1!(b, h, F0)
    d = (F0[1, 1], F0[2, 1], F0[3, 1])
    b .= 0.0
    b[1:3] .= -1 / h .* d
    b[7:9] .= 1 / h .* d
    return b
end

function b_w2!(b, h, F0)
    d = (F0[1, 2], F0[2, 2], F0[3, 2])
    b .= 0.0
    b[1:3] .= -1 / h .* d
    b[7:9] .= 1 / h .* d
    return b
end

function b_w3!(b, h, F0)
    d = (F0[1, 3], F0[2, 3], F0[3, 3])
    b .= 0.0
    b[1:3] .= -1 / h .* d
    b[7:9] .= 1 / h .* d
    return b
end

function i_phi2!(i, h, F0)
    d = (F0[1, 2], F0[2, 2], F0[3, 2])
    i .= 0.0
    i[4:6] .= 1 / 2 .* d
    i[10:12] .= 1 / 2 .* d
    return i
end

function i_phi3!(i, h, F0)
    d = (F0[1, 3], F0[2, 3], F0[3, 3])
    i .= 0.0
    i[4:6] .= 1 / 2 .* d
    i[10:12] .= 1 / 2 .* d
    return i
end

"""
    stiffness(self::FEMMRITBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the material stiffness matrix.
"""
function stiffness(
    self::FEMMRITBeam,
    assembler::ASS,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {ASS<:AbstractSysmatAssembler,T<:Number,TI<:Number}
    fes = self.integdomain.fes
    beamfes = delegateof(fes)
    ecoords0, dofnums = self._ecoords0, self._dofnums
    F0, Te = self._F0, self._Te
    elmat = self._elmat
    _b_phi1 = self._b_phi1
    _b_phi2 = self._b_phi2
    _b_phi3 = self._b_phi3
    _b_w1 = self._b_w1
    _b_w2 = self._b_w2
    _b_w3 = self._b_w3
    _i_phi2 = self._i_phi2
    _i_phi3 = self._i_phi3
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)::Float64
    A, I1, I2, I3, J, A2s, A3s, x1x2_vector, dimensions = properties(fes)
    startassembly!(
        assembler,
        prod(size(elmat)) * count(fes),
        nalldofs(dchi),
        nalldofs(dchi),
    )
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i])
        h, F0 = initial_local_frame!(F0, ecoords0, x1x2_vector[i])
        b_phi1!(_b_phi1, h, F0)
        b_phi2!(_b_phi2, h, F0)
        b_phi3!(_b_phi3, h, F0)
        b_w1!(_b_w1, h, F0)
        b_w2!(_b_w2, h, F0)
        b_w3!(_b_w3, h, F0)
        i_phi2!(_i_phi2, h, F0)
        i_phi3!(_i_phi3, h, F0)
        fill!(elmat, 0.0) # Initialize element matrix
        # Axial Action
        add_n1n2t!(elmat, _b_w1, _b_w1, h * E * A[i])
        # Torsion Action
        add_n1n2t!(elmat, _b_phi1, _b_phi1, h * G * J[i])
        # Bending in the f1-f2 plane
        add_n1n2t!(elmat, _b_w2, _i_phi3, h * G * A2s[i])
        add_n1n2t!(elmat, _b_w2, _b_w2, h * G * A2s[i])
        add_n1n2t!(elmat, _b_phi3, _b_phi3, h * E * I3[i])
        add_n1n2t!(elmat, _i_phi3, _i_phi3, h * G * A2s[i])
        add_n1n2t!(elmat, _i_phi3, _b_w2, h * G * A2s[i])
        # Bending in the f1-f3 plane
        add_n1n2t!(elmat, _b_w3, _i_phi2, h * G * A3s[i])
        add_n1n2t!(elmat, _b_w3, _b_w3, h * G * A3s[i])
        add_n1n2t!(elmat, _b_phi2, _b_phi2, h * E * I2[i])
        add_n1n2t!(elmat, _i_phi2, _i_phi2, h * G * A3s[i])
        add_n1n2t!(elmat, _i_phi2, _b_w3, h * G * A3s[i])
        gatherdofnums!(dchi, dofnums, fes.conn[i]) # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)
    end # Loop over elements
    return makematrix!(assembler)
end

function stiffness(
    self::FEMMRITBeam,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {T<:Number,TI<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi)
end

function _transfmat!(Te, Ft)
    Te[1:3, 1:3] = Te[4:6, 4:6] = Te[7:9, 7:9] = Te[10:12, 10:12] = Ft
    return Te
end

"""
    distribloads_global(self::FEMMRITBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}, fi) where {T<:Number}
    
Compute the load vector due to distributed loads.

Compute the global load vector corresponding to applied distributed
load. Here it means force per unit length of the beam,
in the configuration u1, Rfield1. These are only forces, not moments.

!!! note

    The force intensity must be uniform across the entire element. The
    force intensity is given in the global coordinates.
"""
function distribloads_global(
    self::FEMMRITBeam,
    assembler::ASS,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
    fi,
) where {ASS<:AbstractSysvecAssembler,T<:Number,TI<:Number}
    fes = self.integdomain.fes
    ecoords0, dofnums = self._ecoords0, self._dofnums
    F0, Te = self._F0, self._Te
    elmat = self._elmat
    elvec, elvecf = self._elvec, self._elvecf
    Lforce = fill(0.0, 3)
    ignore = fill(0.0, 0, 0)
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)
    A, I1, I2, I3, J, A2s, A3s, x1x2_vector, dimensions = properties(fes)
    startassembly!(assembler, nalldofs(dchi))
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i])
        h, F0 = initial_local_frame!(F0, ecoords0, x1x2_vector[i])
        force = updateforce!(fi, ignore, ignore, i, 0) # retrieve the applied load
        Lforce = F0' * force # local distributed load components
        elvecf[1] = Lforce[1] * h / 2
        elvecf[2] = Lforce[2] * h / 2
        elvecf[3] = Lforce[3] * h / 2
        elvecf[4] = 0
        elvecf[5] = 0
        elvecf[6] = 0
        elvecf[7] = Lforce[1] * h / 2
        elvecf[8] = Lforce[2] * h / 2
        elvecf[9] = Lforce[3] * h / 2
        elvecf[10] = 0
        elvecf[11] = 0
        elvecf[12] = 0
        _transfmat!(Te, F0)
        mul!(elvec, Te, elvecf)
        gatherdofnums!(dchi, dofnums, fes.conn[i]) # degrees of freedom
        assemble!(assembler, elvec, dofnums)
    end # Loop over elements
    return makevector!(assembler)
end

function distribloads_global(
    self::FEMMRITBeam,
    geom0::NodalField{FFlt},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
    fi,
) where {T<:Number,TI<:Number}
    assembler = SysvecAssembler()
    return distribloads_global(self, assembler, geom0, u1, Rfield1, dchi, fi)
end

end # module
