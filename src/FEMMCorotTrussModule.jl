"""
Module for construction of the linear algebra matrix and vector quantities for
the non-linear co-rotational beam.
"""
module FEMMCorotTrussModule

using LinearAlgebra: norm, Transpose, mul!, dot
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.IntegDomainModule: IntegDomain
import FinEtoolsDeforLinear.MatDeforElastIsoModule: MatDeforElastIso
using ..FESetL2BeamModule: FESetL2Beam, initial_local_frame!


"""
    mutable struct FEMMCorotTruss{ID<:IntegDomain{S} where {S<:FESetL2}, M<:MatDeforElastIso} <: AbstractFEMM

Type for co-rotational beam finite element modeling machine (FEMM).
"""
mutable struct FEMMCorotTruss{ID<:IntegDomain{S} where {S<:FESetL2},M<:MatDeforElastIso} <:
               AbstractFEMM
    integdomain::ID # integration domain data
    material::M # material object
    # The attributes below are buffers used in various operations.
    _ecoords0::FFltMat
    _ecoords1::FFltMat
    _edisp1::FFltMat
    _evel1::FFltMat
    _dofnums::FIntMat
    _F0::FFltMat
    _Ft::FFltMat
    _Te::FFltMat
    _tempelmat1::FFltMat
    _tempelmat2::FFltMat
    _tempelmat3::FFltMat
    _elmat::FFltMat
    _elmatTe::FFltMat
    _elmato::FFltMat
    _elvec::FFltVec
    _OS::FFltMat
end

"""
    FEMMCorotTruss(integdomain::IntegDomain{S, F}, material::MatDeforElastIso) where {S<:FESetL2, F<:Function}

Constructor of finite element machine for nonlinear truss.

Supply the integration domain and the material.
"""
function FEMMCorotTruss(
    integdomain::ID,
    material::M,
) where {ID<:IntegDomain{S} where {S<:FESetL2},M<:MatDeforElastIso}
    typeof(delegateof(integdomain.fes)) <: FESetL2Beam ||
        error("Expected to delegate to FESetL2Beam")
    _ecoords0 = fill(0.0, 2, 3)
    _ecoords1 = fill(0.0, 2, 3)
    _edisp1 = fill(0.0, 2, 3)
    _evel1 = fill(0.0, 2, 3)
    _dofnums = zeros(FInt, 1, 6)
    _F0 = fill(0.0, 3, 3)
    _Ft = fill(0.0, 3, 3)
    _Te = fill(0.0, 6, 6)
    _tempelmat1 = fill(0.0, 6, 6)
    _tempelmat2 = fill(0.0, 6, 6)
    _tempelmat3 = fill(0.0, 6, 6)
    _elmat = fill(0.0, 6, 6)
    _elmatTe = fill(0.0, 6, 6)
    _elmato = fill(0.0, 6, 6)
    _elvec = fill(0.0, 6)
    _OS = fill(0.0, 3, 3)
    return FEMMCorotTruss(
        integdomain,
        material,
        _ecoords0,
        _ecoords1,
        _edisp1,
        _evel1,
        _dofnums,
        _F0,
        _Ft,
        _Te,
        _tempelmat1,
        _tempelmat2,
        _tempelmat3,
        _elmat,
        _elmatTe,
        _elmato,
        _elvec,
        _OS,
    )
end

function properties(fes)
    d = delegateof(fes)
    return d.A
end

function _transfmat!(Te, Ft)
    fill!(Te, 0.0)
    Te[1:3, 1:3] = Te[4:6, 4:6] = Ft
    return Te
end

"""
    local_frame!(Ft, xt)

Compute the initial or current local frame from the positions of the nodes.
"""
function local_frame!(Ft, xt)
    Ft = fill(0.0, 3, 3)
    # Compute the element frame in configuration t
    Ft[:, 1] = (xt[2, :] - xt[1, :])
    Lt = norm(Ft[:, 1])
    Ft[:, 1] /= Lt
    x1x2_vectort = eltype(Ft)[0.0, 0.0, 1.0]'
    if dot(Ft[:, 1], x1x2_vectort) > 0.999
        x1x2_vectort .= (0.0, 1.0, 0.0)
    end
    Ft[:, 3] .= ( # in-place cross product
        -Ft[3, 1] * x1x2_vectort[2] + Ft[2, 1] * x1x2_vectort[3],
        Ft[3, 1] * x1x2_vectort[1] - Ft[1, 1] * x1x2_vectort[3],
        -Ft[2, 1] * x1x2_vectort[1] + Ft[1, 1] * x1x2_vectort[2],
    )
    Ft[:, 3] /= norm(@view Ft[:, 3])
    Ft[:, 2] .= ( # in-place cross product
        -Ft[3, 3] * Ft[2, 1] + Ft[2, 3] * Ft[3, 1],
        Ft[3, 3] * Ft[1, 1] - Ft[1, 3] * Ft[3, 1],
        -Ft[2, 3] * Ft[1, 1] + Ft[1, 3] * Ft[2, 1],
    )
    return Lt, Ft
end

"""
    local_geometric_stiffness!(SM, PN, L)

Compute the local geometric stiffness matrix.

# Arguments
`PN`= axial force in the bar,
`L`= current length of the element,

# Outputs
`SM` = local geometric stiffness matrix, 6 x 6
"""
function local_geometric_stiffness!(SM, PN, L)
    N = PN[1]
    fill!(SM, 0.0)
    for i in [2, 3, 5, 6]
        SM[i, i] = N / L
    end
    for i in [5, 6]
        for j in [2, 3]
            SM[i, j] = SM[j, i] = -N / L
        end
    end
    return SM
end

"""
    local_mass!(MM, A, rho, L)

Mass matrix of the beam.

# Arguments
`A`= cross-sectional area,
`rho`=mass density,
`L`= initial length of the element,

# Outputs
`MM` = local mass matrix, 6 x 6
In the element frame the mass matrix is constant.
"""
function local_mass!(MM, A, rho, L)
    fill!(MM, 0.0)
    c1 = (rho * A * L)::Float64
    MM[1, 1] = c1 / 2
    MM[2, 2] = c1 / 2
    MM[3, 3] = c1 / 2
    MM[4, 4] = c1 / 2
    MM[5, 5] = c1 / 2
    MM[6, 6] = c1 / 2
    return MM
end

"""
    local_stiffness!(SM, E, A, L)

Compute the local elastic stiffness matrix.

# Arguments
`E` = Young's modulus,
`A`= cross-sectional area,
`L`= current length of the element,

# Outputs
`SM` = local stiffness matrix, 6 x 6
"""
function local_stiffness!(SM, E, A, L)
    SM .= 0.0
    SM[1, 1] = E * A / L
    SM[4, 4] = E * A / L
    SM[1, 4] = -E * A / L
    SM[4, 1] = -E * A / L
    return SM
end

natural_forces!(L0, Lt, E, A) = (Lt - L0) / L0 * E * A

"""
    local_forces!(FL, PN)

Compute forces through which the element acts on the nodes in the
local coordinate system.

# Arguments
`PN` = axial force.

# Outputs
FL = vector of forces acting on the nodes in the local coordinate system
"""
function local_forces!(FL, PN)
    FL .= 0.0
    FL[1] = -PN
    FL[4] = +PN
    return FL
end

"""
    mass(
        self::FEMMCorotTruss,
        assembler::ASS,
        geom0::NodalField{GFT},
        u1::NodalField{T},
        Rfield1::NodalField{T},
        dchi::NodalField{TI}
    ) where {ASS<:AbstractSysmatAssembler,GFT<:Number,T<:Number,TI<:Number}

Compute the lumped mass matrix.

This is a general routine for the corotational truss linear-deformation  FEMM.
"""
function mass(
    self::FEMMCorotTruss,
    assembler::ASS,
    geom0::NodalField{GFT},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI}
) where {ASS<:AbstractSysmatAssembler,GFT<:Number,T<:Number,TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums =
        self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, Te = self._F0, self._Ft, self._Te
    elmat, elmatTe = self._elmat, self._elmatTe
    rho = massdensity(self.material)
    A = properties(fes)
    startassembly!(
        assembler,
        size(elmat)..., count(fes),
        nalldofs(dchi),
        nalldofs(dchi),
    )
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i])
        gathervalues_asmat!(u1, edisp1, fes.conn[i])
        ecoords1 .= ecoords0 .+ edisp1
        fill!(elmat, 0.0) # Initialize element matrix
        L0, F0 = local_frame!(Ft, ecoords0)
        Lt, Ft = local_frame!(Ft, ecoords1)
        _transfmat!(Te, Ft)
        local_mass!(elmat, A[i], rho, L0)
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gatherdofnums!(dchi, dofnums, fes.conn[i]) # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)
    end # Loop over elements
    return makematrix!(assembler)
end

function mass(
    self::FEMMCorotTruss,
    geom0::NodalField{GFT},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI}
) where {GFT<:Number,T<:Number,TI<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return mass(self, assembler, geom0, u1, Rfield1, dchi)
end

"""
    gyroscopic(self::FEMMCorotTruss, assembler::ASS, geom0::NodalField{GFT}, u1::NodalField{T}, Rfield1::NodalField{T}, v1::NodalField{T}, dchi::NodalField{TI}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {ASS<:AbstractSysmatAssembler, GFT<:Number, T<:Number, TI<:Number}

Compute the quadratic-inertial-term (gyroscopic) mass matrix

This is a general routine for the abstract linear-deformation  FEMM.
"""
# function gyroscopic(
#     self::FEMMCorotTruss,
#     assembler::ASS,
#     geom0::NodalField{GFT},
#     u1::NodalField{T},
#     Rfield1::NodalField{T},
#     v1::NodalField{T},
#     dchi::NodalField{TI};
#     mass_type = MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA,
# ) where {ASS<:AbstractSysmatAssembler,GFT<:Number,T<:Number,TI<:Number}
#     fes = self.integdomain.fes
#     ecoords0, ecoords1, edisp1, dofnums =
#         self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
#     F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
#     R1I, R1J = self._RI, self._RJ
#     elmat, elmatTe = self._elmat, self._elmatTe
#     dN = self._dN
#     evel1, evel1f = self._evel1, self._evel1f
#     Ge1, Ge2, Ge = self._tempelmat1, self._tempelmat2, self._tempelmat3
#     OmegaTilde = self._elmato
#     OS = self._OS
#     rho = massdensity(self.material)
#     A, I1, I2, I3, J, A2s, A3s, x1x2_vector, dimensions = properties(fes)
#     startassembly!(
#         assembler,
#         size(elmat)..., count(fes),
#         nalldofs(dchi),
#         nalldofs(dchi),
#     )
#     for i = 1:count(fes) # Loop over elements
#         gathervalues_asmat!(geom0, ecoords0, fes.conn[i])
#         gathervalues_asmat!(u1, edisp1, fes.conn[i])
#         ecoords1 .= ecoords0 .+ edisp1
#         R1I[:] .= Rfield1.values[fes.conn[i][1], :]
#         R1J[:] .= Rfield1.values[fes.conn[i][2], :]
#         fill!(elmat, 0.0) # Initialize element matrix
#         L1, Ft, dN = local_frame_and_def!(
#             Ft,
#             dN,
#             F0,
#             FtI,
#             FtJ,
#             ecoords0,
#             x1x2_vector[i],
#             ecoords1,
#             R1I,
#             R1J,
#         )
#         _transfmat!(Te, Ft)
#         L0 = norm(ecoords0[2, :] - ecoords0[1, :])
#         local_mass!(elmat, A[i], I1[i], I2[i], I3[i], rho, L0, mass_type)
#         mul!(elmatTe, elmat, Transpose(Te))
#         mul!(elmat, Te, elmatTe)
#         gathervalues_asmat!(v1, evel1, fes.conn[i])
#         evel1f[:, 1:3] = evel1[:, 1:3] * Ft
#         evel1f[:, 4:6] = evel1[:, 4:6] * Ft
#         Omega =
#             (evel1f[1, 4] + evel1f[2, 4]) / 2 * Ft[:, 1] +
#             (evel1f[1, 3] - evel1f[2, 3]) / L1 * Ft[:, 2] +
#             (evel1f[2, 2] - evel1f[1, 2]) / L1 * Ft[:, 3]
#         skewmat!(OS, Omega)
#         _transfmat!(OmegaTilde, OS)
#         mul!(Ge1, OmegaTilde, elmat)
#         mul!(Ge2, elmat, OmegaTilde)
#         @. Ge = Ge1 - Ge2
#         gatherdofnums!(dchi, dofnums, fes.conn[i]) # degrees of freedom
#         assemble!(assembler, Ge, dofnums, dofnums)
#     end # Loop over elements
#     return makematrix!(assembler)
# end

# function gyroscopic(
#     self::FEMMCorotTruss,
#     geom0::NodalField{GFT},
#     u1::NodalField{T},
#     Rfield1::NodalField{T},
#     v1::NodalField{T},
#     dchi::NodalField{TI};
#     mass_type = MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA,
# ) where {GFT<:Number,T<:Number,TI<:Number}
#     assembler = SysmatAssemblerSparseSymm()
#     return gyroscopic(self, assembler, geom0, u1, Rfield1, v1, dchi; mass_type = mass_type)
# end

"""
    stiffness(self::FEMMCorotTruss, assembler::ASS, geom0::NodalField{GFT}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, GFT<:Number, T<:Number, TI<:Number}

Compute the material stiffness matrix.
"""
function stiffness(
    self::FEMMCorotTruss,
    assembler::ASS,
    geom0::NodalField{GFT},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {ASS<:AbstractSysmatAssembler,GFT<:Number,T<:Number,TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums =
        self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, Te = self._F0, self._Ft, self._Te
    elmat, elmatTe = self._elmat, self._elmatTe
    E = self.material.E
    A = properties(fes)
    startassembly!(
        assembler,
        size(elmat)..., count(fes),
        nalldofs(dchi),
        nalldofs(dchi),
    )
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i])
        gathervalues_asmat!(u1, edisp1, fes.conn[i])
        ecoords1 .= ecoords0 .+ edisp1
        fill!(elmat, 0.0) # Initialize element matrix
        L0, F0 = local_frame!(F0, ecoords0)
        Lt, Ft = local_frame!(Ft, ecoords1)
        _transfmat!(Te, Ft)
        local_stiffness!(elmat, E, A[i], Lt)
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gatherdofnums!(dchi, dofnums, fes.conn[i]) # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)
    end # Loop over elements
    return makematrix!(assembler)
end

function stiffness(
    self::FEMMCorotTruss,
    geom0::NodalField{GFT},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {GFT<:Number,T<:Number,TI<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi)
end


"""
    geostiffness(self::FEMMCorotTruss, assembler::ASS, geom0::NodalField{GFT}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, GFT<:Number, T<:Number, TI<:Number}

Compute the geometric stiffness matrix.
"""
function geostiffness(
    self::FEMMCorotTruss,
    assembler::ASS,
    geom0::NodalField{GFT},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {ASS<:AbstractSysmatAssembler,GFT<:Number,T<:Number,TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums =
        self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, Te = self._F0, self._Ft, self._Te
    elmat, elmatTe = self._elmat, self._elmatTe
    E = self.material.E
    A = properties(fes)
    startassembly!(
        assembler,
        size(elmat)..., count(fes),
        nalldofs(dchi),
        nalldofs(dchi),
    )
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i])
        gathervalues_asmat!(u1, edisp1, fes.conn[i])
        ecoords1 .= ecoords0 .+ edisp1
        fill!(elmat, 0.0) # Initialize element matrix
        L0, F0 = local_frame!(Ft, ecoords0)
        Lt, Ft = local_frame!(Ft, ecoords1)
        _transfmat!(Te, Ft)
        PN = natural_forces!(L0, Lt, E, A[i])
        local_geometric_stiffness!(elmat, PN, Lt)
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gatherdofnums!(dchi, dofnums, fes.conn[i]) # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)
    end # Loop over elements
    return makematrix!(assembler)
end

function geostiffness(
    self::FEMMCorotTruss,
    geom0::NodalField{GFT},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {GFT<:Number,T<:Number,TI<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return geostiffness(self, assembler, geom0, u1, Rfield1, dchi)
end

"""
    restoringforce(self::FEMMCorotTruss, assembler::ASS, geom0::NodalField{GFT}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysvecAssembler, GFT<:Number, T<:Number, TI<:Number}

Compute the vector of the restoring elastic forces
"""
function restoringforce(
    self::FEMMCorotTruss,
    assembler::ASS,
    geom0::NodalField{GFT},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {ASS<:AbstractSysvecAssembler,GFT<:Number,T<:Number,TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums =
        self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, Te = self._F0, self._Ft, self._Te
    elvec = self._elvec
    FL = deepcopy(elvec)
    E = self.material.E
    A = properties(fes)
    startassembly!(assembler, nalldofs(dchi))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i])
        gathervalues_asmat!(u1, edisp1, fes.conn[i])
        ecoords1 .= ecoords0 .+ edisp1
        L0, F0 = local_frame!(Ft, ecoords0)
        Lt, Ft = local_frame!(Ft, ecoords1)
        _transfmat!(Te, Ft)
        PN = natural_forces!(L0, Lt, E, A[i])
        local_forces!(FL, PN)
        mul!(elvec, Te, -FL)
        gatherdofnums!(dchi, dofnums, fes.conn[i]) # degrees of freedom
        assemble!(assembler, elvec, dofnums)
    end # Loop over elements
    return makevector!(assembler)
end

function restoringforce(
    self::FEMMCorotTruss,
    geom0::NodalField{GFT},
    u1::NodalField{T},
    Rfield1::NodalField{T},
    dchi::NodalField{TI},
) where {GFT<:Number,T<:Number,TI<:Number}
    assembler = SysvecAssembler()
    return restoringforce(self, assembler, geom0, u1, Rfield1, dchi)
end

end # module
