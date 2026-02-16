"""
Module for operations on interiors of domains to construct system matrices and
system vectors for linear homogenous shells using the quadrilateral
four-node finite element (Q4RNT). It uses regular membrane and bending stiffness,
and the transverse sheer stiffness is computed with the MITC (DSG) approach.
"""
module FEMMShellQ4RNTModule

using LinearAlgebra: norm, Transpose, mul!, diag, eigen, I, dot, rank
using Statistics: mean
using FinEtools
using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FESetModule: gradN!, nodesperelem, manifdim
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FEMMBaseModule: associategeometry!
import FinEtools.FEMMBaseModule: inspectintegpoints
using FinEtoolsDeforLinear.MatDeforLinearElasticModule:
    tangentmoduli!, update!, thermalstrain!
using FinEtools.MatrixUtilityModule:
    add_btdb_ut_only!, complete_lt!, locjac!, add_nnt_ut_only!, add_btsigma!
using ..FESetShellQ4Module: FESetShellQ4, local_frame!
using ..TransformerModule: TransformerQtEQ, Layup2ElementAngle

const __NN = 4 # number of nodes
const __NDOF = 6 # number of degrees of freedom per node

"""
    mutable struct FEMMShellQ4RNT{ID<:IntegDomain{S} where {S<:FESetQ4}, T<:Real, CS<:CSys{T}, M} <: AbstractFEMM

Type for the finite element modeling machine of the Q4 quadrilateral Flat-Facet
shell with the reduced integration on the shear term and a consistent handling of the
normals. 

Also, the formulation is developed to correctly handle the coupling of twisting
moments and transverse shear (such as in the twisted beam or the Raasch hook
problems) by incorporating "nodal" normals.

The stabilization factor of the shear term of

[2] Mikko Lyly, Rolf Stenberg and Teemu Vihinen, A stable bilinear element for
the Reissner-Mindlin plate model Computer Methods in Applied Mechanics and
Engineering 110 (1993) 343-357 

is incorporated. Refer to expressions (3.12) and (3.13).

The formula for the element to nodal basis transformation is derived similarly 
to the formulation of the robust flat-facet triangle, T3FF.

The following features are incorporated to deal with nodal normals:
- Nodal normals are averages of the normals of elements that meet at a node.
- A crease in the surface is taken into account. In that case the normal are not
  averaged across the crease. At the nodes along the crease every element uses
  the normal to its surface instead of the nodal normal.

Configuration:

These attributes of the FEMM can be set after it's been created.
- `drilling_stiffness_scale`: multiplier of the generalized stiffness coefficient
- `threshold_angle`: angle in degrees. If a nodal normal subtends angle bigger
  then this threshold, the nodal normal at that note is marked as invalid.
- `mult_el_size`: multiplier of the square of the element size, used to control
  transverse shear stiffness.
"""
mutable struct FEMMShellQ4RNT{ID<:IntegDomain{S} where {S<:FESetQ4}, T<:Real, CS<:CSys{T}, M} <: AbstractFEMM
    integdomain::ID # integration domain data
    mcsys::CS # updater of the material orientation matrix
    material::M # material object.
    drilling_stiffness_scale::T
    threshold_angle::T
    mult_el_size::T
    _associatedgeometry::Bool
    _normals::Matrix{T}
    _normal_valid::Vector{Bool}
end

_number_type(femm::FEMMShellQ4RNT{ID, T, CS}) where {ID, T, CS} = T

# Prepare functions to return buffers of various sorts

function _loc(ft::Type{T}) where {T<:Real}
    return fill(zero(ft), 1, 3)
end

function _J(ft::Type{T}) where {T<:Real}
    return fill(zero(ft), 3, 2)
end

function _ecoords(ft::Type{T}) where {T<:Real}
    return fill(zero(ft), __NN, 3)
end

function _edisp(ft::Type{T}) where {T<:Real}
    return fill(zero(ft), __NN * __NDOF)
end

function _ecoords_e(ft::Type{T}) where {T<:Real}
    return fill(zero(ft), __NN, 2)
end

function _edisp_e(ft::Type{T}) where {T<:Real}
    return fill(zero(ft), __NN * __NDOF)
end

function _gradN_e(ft::Type{T}) where {T<:Real}
    return fill(zero(ft), __NN, 2)
end

function _dofnums(it::Type{IT}) where {IT<:Integer}
    return fill(zero(it), 1, __NN * __NDOF)
end

function _Bs(ft::Type{T}) where {T<:Real}
    _Bm = fill(zero(ft), 3, __NN * __NDOF)
    _Bb = fill(zero(ft), 3, __NN * __NDOF)
    _Bs = fill(zero(ft), 2, __NN * __NDOF)
    _DpsBmb = similar(_Bm)
    _DtBs = similar(_Bs)
    return _Bm, _Bb, _Bs, _DpsBmb, _DtBs
end

function _E_G(ft::Type{T}) where {T<:Real}
    return fill(zero(ft), 3, 3)
end

function _A_Es(ft::Type{T}) where {T<:Real}
    return [fill(zero(ft), 3, 3) for i in 1:__NN]
end

function _nvalid()
    return fill(false, __NN)
end

function _T(ft::Type{T}) where {T<:Real}
    return fill(zero(ft), __NN * __NDOF, __NN * __NDOF)
end

function _elmat(ft::Type{T}) where {T<:Real}
    return fill(zero(ft), __NN * __NDOF, __NN * __NDOF)
end

"""
    FEMMShellQ4RNT(
        integdomain::ID,
        mcsys::CSys,
        material::M,
    ) where {IntegDomain{S} where {S<:FESetQ4}, CS<:CSys, M}

Constructor of the Q4RNT shell FEMM.
"""
function FEMMShellQ4RNT(
    integdomain::ID,
    mcsys::CS,
    material::M,
) where {ID<:IntegDomain{S} where {S<:FESetQ4}, T<:Real, CS<:CSys{T}, M}
    _nnmax = 0
    for j in eachindex(integdomain.fes)
        for k in eachindex(integdomain.fes.conn[j])
            _nnmax = max(_nnmax, integdomain.fes.conn[j][k])
        end
    end
    _normals = fill(zero(T), _nnmax, 3)
    _normal_valid = fill(true, _nnmax)
    
    @assert delegateof(integdomain.fes) === FESetShellQ4()

    return FEMMShellQ4RNT(
        integdomain,
        mcsys,
        material,
        # drilling_stiffness_scale::Float64
        # threshold_angle::Float64
        # mult_el_size::Float64
        T(1.0),
        T(30.0),
        T(0.05),
        false,
        _normals,
        _normal_valid,
    )
end

function isoparametric!(E_G, XYZ, J, feid, qpid)
    return _e_g!(E_G, J)
end

function _compute_nodal_normal!(n, mcsys::CSys, XYZ, J, feid, qpid)
    updatecsmat!(mcsys, reshape(XYZ, 1, 3), J, feid, qpid)
    n[:] .= csmat(mcsys)[:, 3]
    return n
end

"""
    FEMMShellQ4RNT(
        integdomain::ID,
        material::M,
    ) where {ID<:IntegDomain{S} where {S<:FESetQ4}, M}

Constructor of the Q4RNT shell FEMM.
"""
function FEMMShellQ4RNT(
    integdomain::ID,
    material::M,
) where {ID<:IntegDomain{S} where {S<:FESetQ4}, M}
    return FEMMShellQ4RNT(integdomain, CSys(3, 3, zero(Float64), isoparametric!), material)
end

"""
    make(integdomain, material)

Make a Q4RNT FEMM from the integration domain,  and a material.
Default isoparametric method for computing the normais is used.
"""
function make(integdomain, material)
    return FEMMShellQ4RNT(integdomain, material)
end

"""
    make(integdomain, mcsys, material)

Make a Q4RNT FEMM from the integration domain, a coordinate system to define the
orientation of the normals, and a material.
"""
function make(integdomain, mcsys, material)
    return FEMMShellQ4RNT(integdomain, mcsys, material)
end

function _e_g!(E_G, J)
    # J[:, 1] is the tangent to the coordinate curve 1
    E_G[:, 1] .= J[:, 1]
    n = sqrt(E_G[1, 1]^2 + E_G[2, 1]^2 + E_G[3, 1]^2)
    E_G[:, 1] ./= n
    # J[:, 2] is the tangent to the coordinate curve 2
    # Now compute the normal
    E_G[1, 3] = -E_G[3, 1] * J[2, 2] + E_G[2, 1] * J[3, 2]
    E_G[2, 3] = E_G[3, 1] * J[1, 2] - E_G[1, 1] * J[3, 2]
    E_G[3, 3] = -E_G[2, 1] * J[1, 2] + E_G[1, 1] * J[2, 2]
    n = sqrt(E_G[1, 3]^2 + E_G[2, 3]^2 + E_G[3, 3]^2)
    E_G[:, 3] ./= n
    E_G[1, 2] = -E_G[3, 3] * E_G[2, 1] + E_G[2, 3] * E_G[3, 1]
    E_G[2, 2] = E_G[3, 3] * E_G[1, 1] - E_G[1, 3] * E_G[3, 1]
    E_G[3, 2] = -E_G[2, 3] * E_G[1, 1] + E_G[1, 3] * E_G[2, 1]
    return E_G
end

function _shell_material_stiffness(material)
    # Compute the two material stiffness matrices: the plane stress matrix, and
    # the transverse shear matrix
    D = fill(0.0, 6, 6)
    t::FFlt, dt::FFlt, loc::FFltMat, label::FInt = 0.0, 0.0, [0.0 0.0 0.0], 0
    tangentmoduli!(material, D, t, dt, loc, label)
    Dps = fill(0.0, 3, 3)
    Dps[1:2, 1:2] =
        D[1:2, 1:2] - (reshape(D[1:2, 3], 2, 1) * reshape(D[3, 1:2], 1, 2)) / D[3, 3]
    ix = [1, 2, 4]
    for i = 1:3
        Dps[3, i] = Dps[i, 3] = D[4, ix[i]]
    end
    Dt = fill(0.0, 2, 2)
    ix = [5, 6]
    for i = 1:2
        Dt[i, i] = D[ix[i], ix[i]]
    end
    return Dps, Dt
end

# struct LocalDerivatives{FT<:Real}
#     G::Matrix{FT}
#     tmp::Vector{FT}
# end
    
# function LocalDerivatives(ft::Type{T}) where {T<:Real}
#     return LocalDerivatives(fill(zero(ft), 2, 2), fill(zero(ft), 2))
# end

# (o::LocalDerivatives)(gradN_e, J, E_G, gradNparams) = let
#     mul!(o.G, transpose(J), J)
#     detG = det(o.G)
    
#     if detG ≈ 0.0
#         error("Singular metric matrix in _gradN_e!")
#         return
#     end
    
#     invG = inv(G)
#     tmp = zeros(2)
    
#     for a in 1:__NN
#         tmp[1] = gradNparams[a, 1]
#         tmp[2] = gradNparams[a, 2]
#         g = J * (invG * tmp)
#         gradN_e[a, 1] = dot(E_G[:, 1], g)
#         gradN_e[a, 2] = dot(E_G[:, 2], g)
#     end
# end

# TODO optimize allocations
function _gradN_e!(gradN_e, J0, E_G, gradNparams)
    G = J0' * J0
    detG = G[1, 1] * G[2, 2] - G[1, 2] * G[2, 1]
    
    if detG ≈ 0.0
        fill!(gradN_e, 0.0)
        return
    end
    
    invG = inv(G)
    E2 = E_G[:, 1:2]
    tmp = zeros(2)
    
    for a in 1:__NN
        tmp[1] = gradNparams[a, 1]
        tmp[2] = gradNparams[a, 2]
        g = J0 * (invG * tmp)
        gradN_e[a, 1] = dot(E2[:, 1], g)
        gradN_e[a, 2] = dot(E2[:, 2], g)
    end
end

struct _NodalTriadsE{FT<:Real}
    r::Vector{FT}
    nk_e::Vector{FT}
    nk::Vector{FT}
    f3_e::Vector{FT}
end

function _NodalTriadsE(ft::Type{T}) where {T<:Real}
    _NodalTriadsE(fill(zero(ft), 3), fill(zero(ft), 3), fill(zero(ft), 3), [zero(ft), zero(ft), zero(ft)+1.0])
end

(o::_NodalTriadsE)(A_Es, nvalid, E_G, normals, normal_valid, c) = begin
    # Components of nodal cartesian ordinate systems such that the third
    # direction is the direction of the nodal normal, and the angle to rotate
    # the element normal into the nodal normal is as short as possible; these
    # components are given on the local element coordinate system basis
    # vectors (basis E). If the angle of rotation is excessively small, the
    # nodal matrix is the identity. 
    # 
    # The array `A_Es` holds the triads that are the nodal-to-element
    # transformation matrices. The array `nvalid` indicates for the three nodes
    # which of the normals is valid (`true`).

    for k in eachindex(A_Es)
        o.nk .= vec(view(normals, c[k], :))
        nvalid[k] = normal_valid[c[k]]
        if nvalid[k]
            # If the nodal normal is valid, pull it back into the element frame.
            mul!(o.nk_e, transpose(E_G), o.nk)
        else
            # Otherwise the element normal replaces the nodal normal.
            o.nk_e .= 0.0
            o.nk_e[3] = 1.0
        end
        cross3!(o.r, o.f3_e, o.nk_e)
        if norm(o.r) > 1.0e-12
            rotmat3!(A_Es[k], o.r)
        else
            # A_Es[k] .= I(3)
            A_Es[k] .= 0.0
            A_Es[k][1, 1] = A_Es[k][2, 2] = A_Es[k][3, 3] = 1.0
        end
    end
    return A_Es, nvalid
end

struct _TransfmatGToA{FT<:Real}
    Tblock::Matrix{FT}
end

function _TransfmatGToA(ft::Type{T}) where {T<:Real}
    _TransfmatGToA(fill(zero(ft), 3, 3))
end

(o::_TransfmatGToA)(T, A_Es, E_G) = begin
    # Global-to-nodal transformation matrix. 

    # The 3x3 blocks consist of the nodal triad expressed on the global basis.
    # The nodal basis vectors in the array `A_Es[i]` are expressed on the
    # element basis, and are then rotated with `E_G` into the global coordinate
    # system.

    # Output
    # - `T` = transformation matrix, input in the global basis, output in the
    #   nodal basis
    T .= zero(eltype(T))
    for i = 1:__NN
        mul!(o.Tblock, transpose(A_Es[i]), transpose(E_G))
        offset = (i - 1) * __NDOF
        r = offset+1:offset+3
        @. T[r, r] = o.Tblock
        r = offset+4:offset+6
        @. T[r, r] = o.Tblock
    end
    return T
end

function _transfmat_a_to_e!(T, A_Es, gradN_e)
    # Nodal-to-element transformation matrix. 
    # - `A_Es` = matrix with the nodal triad vectors in columns, components on
    #   the element basis. Its transpose is the element triad on the nodal
    #   basis vectors.
    # - `gradN_e` = basis function gradients on the element basis.
    # Output
    # - `T` = transformation matrix, input in the nodal basis, output in the
    #   element basis: t_e = T * t_n
    # Rotation degrees of freedom: The drilling rotation of the mid surface
    # produced by the 1/2*(v,x - u,y) effect is linked to the out of plane
    # rotations.

    T .= zero(eltype(T))
    for i = 1:__NN
        roffst = (i - 1) * __NDOF
        iA_E = A_Es[i]
        iA_E_33 = iA_E[3, 3]
        # T[r, r] .= iA_E
        for cl = 1:3
            for rw = 1:3
                T[roffst+rw, roffst+cl] = iA_E[rw, cl]
            end
        end
        for cl = 1:2
            for rw = 1:2
                T[roffst+3+rw, roffst+3+cl] =
                    iA_E[rw, cl] - (1 / iA_E_33) * iA_E[rw, 3] * iA_E[cl, 3]
            end
        end
        m1 = (1 / iA_E_33) * iA_E[1, 3]
        m2 = (1 / iA_E_33) * iA_E[2, 3]
        for j = 1:__NN
            coffst = (j - 1) * __NDOF
            for k = 1:3
                a3 = 1 / 2 * (iA_E[2, k] * gradN_e[j, 1] - iA_E[1, k] * gradN_e[j, 2])
                T[roffst+4, coffst+k] += m1 * a3
                T[roffst+5, coffst+k] += m2 * a3
            end
        end
    end
    return T
end

"""
    associategeometry!(self::FEMMShellQ4RNT,  geom0::NodalField{FFlt})

Associate geometry with the FEMM. 

In this case it means evaluate the nodal normals.
"""
function associategeometry!(self::FEMMShellQ4RNT, geom0::NodalField{FFlt})
    threshold_angle = self.threshold_angle
    FT = _number_type(self)
    J = _J(FT)
    loc = _loc(FT)
    E_G = _E_G(FT)
    normals, normal_valid = self._normals, self._normal_valid
    nnormal = fill(zero(FT), 3)
    nodal_rule = NodalTensorProductRule(2)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain, nodal_rule)
    ecoords = _ecoords(FT)

    normals .= 0.0
    normal_valid .= true

    fes = finite_elements(self)
    for el in eachindex(fes)
        gathervalues_asmat!(geom0, ecoords, fes.conn[el])
        for j in 1:__NN
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianmdim(self.integdomain, J, loc, fes.conn[el], Ns[j], 2)
            _compute_nodal_normal!(nnormal, self.mcsys, loc, J, el, j)
            normals[fes.conn[el][j], :] .+= Jac * nnormal
        end
    end

    # Normalize to unit length
    for j in axes(normals, 1)
        nn = norm(normals[j, :])
        if nn > 0.0
            normals[j, :] ./= nn
        end
    end

    # Now perform a second pass. If the nodal normal differs by a substantial
    # amount from the element normals, zero out the nodal normal to indicate
    # that the nodal normal should not be used at that vertex. 
    ntolerance = 1 - sqrt(1 - sin(threshold_angle / 180 * pi)^2)
    for el in eachindex(self.integdomain.fes)
        gathervalues_asmat!(geom0, ecoords, fes.conn[el])
        for j in 1:__NN
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            _compute_nodal_normal!(nnormal, self.mcsys, loc, J, el, j)
            n = fes.conn[el][j]
            nd = dot(normals[n, :], nnormal)
            if nd < 1 - ntolerance
                normal_valid[n] = false
            end
        end
    end

    # Mark the engine as being associated with geometry
    self._associatedgeometry = true
    return self
end


"""
    num_normals(self::FEMMShellQ4RNT)

Compute the summary of the nodal normals.
"""
function num_normals(self::FEMMShellQ4RNT)
    normals, normal_valid = self._normals, self._normal_valid
    total_normals = length(normal_valid)
    total_invalid_normals = length([v for v in normal_valid if v != true])
    return total_normals, total_invalid_normals
end

struct _Bmmat{FT<:Real}
    tempBm::Matrix{FT}
end

function _Bmmat(ft::Type{T}) where {T<:Real}
    _Bmmat(fill(zero(ft), 3, __NN * __NDOF))
end

(o::_Bmmat)(Bm, gradN, T) = begin
    o.tempBm .= 0.0
    for i in 1:__NN
        off = (i-1)*__NDOF
        o.tempBm[1, off + 1] = gradN[i,1]
        o.tempBm[2, off + 2] = gradN[i,2]
        o.tempBm[3, off + 1] = gradN[i,2]
        o.tempBm[3, off + 2] = gradN[i,1]
    end
    mul!(Bm, o.tempBm, T)
end

struct _Bbmat{FT<:Real}
    tempBb::Matrix{FT}
end

function _Bbmat(ft::Type{T}) where {T<:Real}
    _Bbmat(fill(zero(ft), 3, __NN * __NDOF))
end

(o::_Bbmat)(Bb, gradN, T) = begin
    o.tempBb .= 0.0
    for i in 1:__NN
        off = (i-1)*__NDOF
        o.tempBb[1, off + 5] = gradN[i,1]
        o.tempBb[2, off + 4] = -gradN[i,2]
        o.tempBb[3, off + 4] = -gradN[i,1]
        o.tempBb[3, off + 5] = gradN[i,2]
    end
    mul!(Bb, o.tempBb, T)
end

# TODO eliminate allocations
function _ecoords_e!(ecoords_e, ecoords, E_G)
    centroid = mean(ecoords, dims=1)'
    for j in axes(ecoords_e, 1)
        for k in axes(ecoords_e, 2)
            ecoords_e[j, k] = dot(ecoords[j, :] - centroid, E_G[:, k])
        end
    end
    return ecoords_e
end

struct _Bsmat{FT<:Real}
    tempBs::Matrix{FT}
end

function _Bsmat(ft::Type{T}) where {T<:Real}
    _Bsmat(fill(zero(ft), 2, __NN * __NDOF))
end

(o::_Bsmat)(Bs, ecoords_e, rs, T) = begin
    o.tempBs .= 0.0
    X1, Y1 = ecoords_e[1, 1], ecoords_e[1, 2]
    X2, Y2 = ecoords_e[2, 1], ecoords_e[2, 2]
    X3, Y3 = ecoords_e[3, 1], ecoords_e[3, 2]
    X4, Y4 = ecoords_e[4, 1], ecoords_e[4, 2]
    r, s = rs[1], rs[2]
    J11 = (X1 * (s - 1) / 4 - X2 * (s - 1) / 4 + X3 * (s + 1) / 4 - X4 * (s + 1) / 4)
    J21 = (Y1 * (s - 1) / 4 - Y2 * (s - 1) / 4 + Y3 * (s + 1) / 4 - Y4 * (s + 1) / 4)
    J12 = (X1 * (r - 1) / 4 - X2 * (r + 1) / 4 + X3 * (r + 1) / 4 - X4 * (r - 1) / 4)
    J22 = (Y1 * (r - 1) / 4 - Y2 * (r + 1) / 4 + Y3 * (r + 1) / 4 - Y4 * (r - 1) / 4)
    A = sqrt(J11^2 + J21^2)
    B = sqrt(J12^2 + J22^2)
    ca = J11 / A
    sa = J21 / A    
    cb = J12 / B
    sb = J22 / B
    detJ = J11 * J22 - J12 * J21
    Ax = X1 - X2 - X3 + X4
    Ay = Y1 - Y2 - Y3 + Y4
    Bx = X1 - X2 + X3 - X4
    By = Y1 - Y2 + Y3 - Y4
    Cx = X1 + X2 - X3 - X4
    Cy = Y1 + Y2 - Y3 - Y4
    # gxz
    o.tempBs[1, 3] += (-sa * (r + 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) + sb * (s + 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (16 * detJ)
    o.tempBs[1, 4] += (sa * (Y1 - Y4) * (r + 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) - sb * (Y1 - Y2) * (s + 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[1, 5] += (-sa * (X1 - X4) * (r + 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) + sb * (X1 - X2) * (s + 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[1, 9] += (sa * (r - 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) - sb * (s + 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (16 * detJ)
    o.tempBs[1, 10] += (-sa * (Y2 - Y3) * (r - 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) - sb * (Y1 - Y2) * (s + 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[1, 11] += (sa * (X2 - X3) * (r - 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) + sb * (X1 - X2) * (s + 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[1, 15] += (-sa * (r - 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) + sb * (s - 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (16 * detJ)
    o.tempBs[1, 16] += (-sa * (Y2 - Y3) * (r - 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) - sb * (Y3 - Y4) * (s - 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[1, 17] += (sa * (X2 - X3) * (r - 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) + sb * (X3 - X4) * (s - 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[1, 21] += (sa * (r + 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) - sb * (s - 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (16 * detJ)
    o.tempBs[1, 22] += (sa * (Y1 - Y4) * (r + 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) - sb * (Y3 - Y4) * (s - 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[1, 23] += (-sa * (X1 - X4) * (r + 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) + sb * (X3 - X4) * (s - 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    # gyz
    o.tempBs[2, 3] += (ca * (r + 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) - cb * (s + 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (16 * detJ)
    o.tempBs[2, 4] += (-ca * (Y1 - Y4) * (r + 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) + cb * (Y1 - Y2) * (s + 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[2, 5] += (ca * (X1 - X4) * (r + 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) - cb * (X1 - X2) * (s + 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[2, 9] += (-ca * (r - 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) + cb * (s + 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (16 * detJ)
    o.tempBs[2, 10] += (ca * (Y2 - Y3) * (r - 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) + cb * (Y1 - Y2) * (s + 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[2, 11] += (-ca * (X2 - X3) * (r - 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) - cb * (X1 - X2) * (s + 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[2, 15] += (ca * (r - 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) - cb * (s - 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (16 * detJ)
    o.tempBs[2, 16] += (ca * (Y2 - Y3) * (r - 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) + cb * (Y3 - Y4) * (s - 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[2, 17] += (-ca * (X2 - X3) * (r - 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) - cb * (X3 - X4) * (s - 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[2, 21] += (-ca * (r + 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) + cb * (s - 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (16 * detJ)
    o.tempBs[2, 22] += (-ca * (Y1 - Y4) * (r + 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) + cb * (Y3 - Y4) * (s - 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    o.tempBs[2, 23] += (ca * (X1 - X4) * (r + 1) * sqrt((Ax + Bx * s)^2 + (Ay + By * s)^2) - cb * (X3 - X4) * (s - 1) * sqrt((Bx * r + Cx)^2 + (By * r + Cy)^2)) / (32 * detJ)
    mul!(Bs, o.tempBs, T)
end

function _drilling_penalty_kavg(elmat, normals, normal_valid, conn)
    I3 = Matrix{Float64}(I, 3, 3)
    tangential = Float64[]
    
    for k in 1:__NN
        nnode = Int(conn[k])
        if !normal_valid[nnode]
            continue
        end
        
        nvec = @view normals[nnode, :]
        nn = norm(nvec)
        if nn == 0.0
            continue
        end
        nvec_normalized = nvec / nn
        
        r = ((k - 1) * __NDOF + 4):((k - 1) * __NDOF + 6)
        Krr = @view elmat[r, r]
        P = I3 - nvec_normalized * nvec_normalized'
        Kt = P * Krr * P
        
        push!(tangential, max(0.0, sum(diag(Kt)) / 2.0))
    end
    
    if isempty(tangential)
        return 0.0
    end
    
    return mean(tangential)
end

# TODO optimize allocations
function _add_drilling_stiffness!(elmat, normals, normal_valid, conn, drilling_stiffness_scale)
    if drilling_stiffness_scale != 0.0
        kavg = _drilling_penalty_kavg(elmat, normals, normal_valid, conn) * drilling_stiffness_scale
        if kavg != 0.0
            for k in 1:__NN
                nnode = conn[k]
                if !normal_valid[nnode]
                    continue
                end
                nvec = @view normals[nnode, :]
                if norm(nvec) == 0.0
                    continue
                end
                Ke = kavg .* (nvec * nvec')
                r = ((k-1)*__NDOF+4):((k-1)*__NDOF+6)
                elmat[r, r] .+= Ke
            end
        end
    end
end

function _quad_diameter(ecoords)
    max_dist = 0.0
    for i in 2:__NN
        dist = ((ecoords[i, 1] - ecoords[1, 1])^2 + 
                (ecoords[i, 2] - ecoords[1, 2])^2 +
                (ecoords[i, 3] - ecoords[1, 3])^2)
        max_dist = max(max_dist, dist)
    end
    return sqrt(max_dist)
end

"""
    stiffness(self::FEMMShellQ4RNT, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{TI}, Rfield1::NodalField{TI}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, TI<:Number}

Compute the material stiffness matrix.
"""
function stiffness(
    self::FEMMShellQ4RNT,
    assembler::ASS,
    geom0::NodalField{FFlt},
    u1::NodalField{TI},
    Rfield1::NodalField{TI},
    dchi::NodalField{TI},
) where {ASS<:AbstractSysmatAssembler,TI<:Number}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    normals, normal_valid = self._normals, self._normal_valid
    FT = _number_type(self)
    loc, J = _loc(FT), _J(FT)
    ecoords, ecoords_e, gradN_e = _ecoords(FT), _ecoords_e(FT), _gradN_e(FT)
    dofnums = _dofnums(eltype(dchi.dofnums)) 
    E_G, A_Es, nvalid, T = _E_G(FT), _A_Es(FT), _nvalid(), _T(FT)
    elmat = _elmat(FT)
    full_rule = self.integdomain.integration_rule.rule1
    fi_npts, fi_Ns, fi_gradNparams, fi_w, fi_pc = integrationdata(self.integdomain, full_rule)
    _nodal_triads_e! = _NodalTriadsE(FT)
    _transfmat_g_to_a! = _TransfmatGToA(FT)
    Bm, Bb, Bs, DpsBmb, DtBs = _Bs(FT)
    bmmat! = _Bmmat(FT); bbmat! = _Bbmat(FT); bsmat! = _Bsmat(FT)
    Dps, Dt = _shell_material_stiffness(self.material)
    T = _T(FT); Tae = _T(FT); Tga = _T(FT)
    scf = 5 / 6  # shear correction factor
    Dt .*= scf
    
    drilling_stiffness_scale = self.drilling_stiffness_scale
    mult_el_size = self.mult_el_size
    startassembly!(
        assembler,
        size(elmat)..., count(fes),
        nalldofs(dchi),
        nalldofs(dchi),
    )
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords, fes.conn[i])
        h = _quad_diameter(ecoords)
        fill!(elmat, 0.0) # Initialize element matrix
        for j in 1:fi_npts 
            locjac!(loc, J, ecoords, fi_Ns[j], fi_gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], fi_Ns[j])
            _e_g!(E_G, J)
            _ecoords_e!(ecoords_e, ecoords, E_G)
            _gradN_e!(gradN_e, J, E_G, fi_gradNparams[j])
            t = self.integdomain.otherdimension(loc, fes.conn[i], fi_Ns[j])
            _nodal_triads_e!(A_Es, nvalid, E_G, normals, normal_valid, fes.conn[i])
            _transfmat_g_to_a!(Tga, A_Es, E_G)
            _transfmat_a_to_e!(Tae, A_Es, gradN_e)
            mul!(T, Tae, Tga)
            bmmat!(Bm, gradN_e, T)
            add_btdb_ut_only!(elmat, Bm, t * Jac * fi_w[j], Dps, DpsBmb)
            bbmat!(Bb, gradN_e, T)
            add_btdb_ut_only!(elmat, Bb, (t^3 / 12.0) * Jac * fi_w[j], Dps, DpsBmb)
            bsmat!(Bs, ecoords_e, fi_pc[j, :], T)
            Lylyetal = t^2 / (t^2 + mult_el_size * h^2)
            add_btdb_ut_only!(elmat, Bs,
                t  * Lylyetal * Jac * fi_w[j],
                Dt, DtBs)
        end
        # Complete the elementwise matrix by filling in the lower triangle
        complete_lt!(elmat)
        # Drilling stiffness
        _add_drilling_stiffness!(elmat, normals, normal_valid, fes.conn[i], drilling_stiffness_scale)
        # Assembly
        gatherdofnums!(dchi, dofnums, fes.conn[i])
        assemble!(assembler, elmat, dofnums, dofnums)
    end # Loop over elements
    return makematrix!(assembler)
end

function stiffness(
    self::FEMMShellQ4RNT,
    geom0::NodalField{FFlt},
    u1::NodalField{TI},
    Rfield1::NodalField{TI},
    dchi::NodalField{TI},
) where {TI<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi)
end


"""
    mass(self::FEMMShellQ4RNT,  assembler::A,  geom::NodalField{FFlt}, dchi::NodalField{TI}) where {A<:AbstractSysmatAssembler, TI<:Number}

Compute the diagonal (lumped) mass matrix.

The mass matrix can be expected to be non-singular.
"""
function mass(
    self::FEMMShellQ4RNT,
    assembler::A,
    geom0::NodalField{FFlt},
    dchi::NodalField{TI},
) where {A<:AbstractSysmatAssembler,TI<:Number}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    FT = _number_type(self)
    loc, J = _loc(FT), _J(FT)
    ecoords = _ecoords(FT)
    dofnums = _dofnums(eltype(dchi.dofnums)) 
    elmat = _elmat(FT)
    rho = massdensity(self.material) # mass density
    full_rule = self.integdomain.integration_rule.rule1
    fi_npts, fi_Ns, fi_gradNparams, fi_w, fi_pc = integrationdata(self.integdomain, full_rule)
    startassembly!(
        assembler,
        size(elmat)..., count(fes),
        nalldofs(dchi),
        nalldofs(dchi),
    )
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords, fes.conn[i])
        # Calculate the element mass
        tmass = 0.0
        rmass = 0.0
        for j in 1:fi_npts
            locjac!(loc, J, ecoords, fi_Ns[j], fi_gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], fi_Ns[j])
            t = self.integdomain.otherdimension(loc, fes.conn[i], fi_Ns[j])
            tmass += rho * t * Jac * fi_w[j]
            rmass += rho * t^3 / 12 * Jac * fi_w[j]
        end
        tmass = tmass / __NN
        rmass = rmass / __NN
        # Fill the elementwise matrix in the global basis.
        fill!(elmat, 0.0) # Initialize element matrix
        for k in 1:__NN
            # Translation degrees of freedom
            for d in 1:3
                c = (k - 1) * __NDOF + d
                elmat[c, c] += tmass
            end
            # Bending degrees of freedom
            for d in 4:6
                c = (k - 1) * __NDOF + d
                elmat[c, c] += rmass
            end
        end
        # Assemble
        gatherdofnums!(dchi, dofnums, fes.conn[i])# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function mass(
    self::FEMMShellQ4RNT,
    geom::NodalField{FFlt},
    u::NodalField{TI},
) where {TI<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return mass(self, assembler, geom, u)
end


"""
    inspectintegpoints(self::AbstractFEMMDeforLinear,
      geom::NodalField{FFlt},  u::NodalField{TI},
      dT::NodalField{FFlt},
      felist::FIntVec,
      inspector::F,  idat, quantity=:Cauchy;
      context...) where {TI<:Number, F<:Function}

Inspect integration point quantities.

- `geom` - reference geometry field
- `u` - displacement+rotation field
- `dT` - temperature difference field
- `felist` - indexes of the finite elements that are to be inspected:
     The fes to be included are: `fes[felist]`.
- `context`    - structure: see the update!() method of the material.
- `inspector` - function with the signature
        `idat = inspector(idat, j, conn, x, out, loc);`
   where
    `idat` - a structure or an array that the inspector may
           use to maintain some state,  for instance minimum or maximum of
           stress, `j` is the element number, `conn` is the element connectivity,
           `out` is the output of the update!() method,  `loc` is the location
           of the integration point in the *reference* configuration.
### Return
The updated inspector data is returned.
"""
function inspectintegpoints(
    self::FEMMShellQ4RNT,
    geom0::NodalField{FFlt},
    u::NodalField{TI},
    dT::NodalField{FFlt},
    felist::FIntVec,
    inspector::F,
    idat,
    quantity = :Cauchy;
    context...,
) where {TI<:Number,F<:Function}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    label = self.integdomain.fes.label
    normals, normal_valid = self._normals, self._normal_valid
    FT = _number_type(self)
    loc, J = _loc(FT), _J(FT)
    ecoords, ecoords_e, gradN_e = _ecoords(FT), _ecoords_e(FT), _gradN_e(FT)
    E_G, A_Es, nvalid, T = _E_G(FT), _A_Es(FT), _nvalid(), _T(FT)
    edisp = _edisp(FT)
    _nodal_triads_e! = _NodalTriadsE(FT)
    _transfmat_g_to_a! = _TransfmatGToA(FT)
    T = _T(FT); Tae = _T(FT); Tga = _T(FT)
    full_rule = self.integdomain.integration_rule.rule1
    fi_npts, fi_Ns, fi_gradNparams, fi_w, fi_pc = integrationdata(self.integdomain, full_rule)
    reduced_rule = self.integdomain.integration_rule.rule2
    ri_npts, ri_Ns, ri_gradNparams, ri_w, ri_pc = integrationdata(self.integdomain, reduced_rule)
    Bm, Bb, Bs, DpsBmb, DtBs = _Bs(FT)
    bmmat! = _Bmmat(FT); bbmat! = _Bbmat(FT); bsmat! = _Bsmat(FT)
    lla = Layup2ElementAngle()
    Dps, Dt = _shell_material_stiffness(self.material)
    scf = 5 / 6  # shear correction factor
    Dt .*= scf
    mult_el_size = self.mult_el_size
    out = fill(0.0, 3)
    o2_e = fill(0.0, 2, 2)
    # Sort out  the output requirements
    outputcsys = self.mcsys # default: report the stresses in the material coord system
    for apair in pairs(context)
        sy, val = apair
        if sy == :outputcsys
            outputcsys = val
        end
    end
    BENDING_MOMENT, TRANSVERSE_SHEAR, MEMBRANE_FORCE = 1, 2, 3
    quant = BENDING_MOMENT
    if quantity == :bending || quantity == :moment || quantity == :bending_moment
        quant = BENDING_MOMENT
    end
    if quantity == :transverse_shear || quantity == :transverse || quantity == :shear
        quant = TRANSVERSE_SHEAR
    end
    if quantity == :membrane_force || quantity == :membrane
        quant = MEMBRANE_FORCE
    end
    # Loop over  all the elements and all the quadrature points within them
    for ilist in eachindex(felist) # Loop over elements
        i = felist[ilist]
        gathervalues_asmat!(geom0, ecoords, fes.conn[i])
        gathervalues_asvec!(u, edisp, fes.conn[i])
        if quant == TRANSVERSE_SHEAR
            for j in 1:ri_npts # Loop over quadrature points
                locjac!(loc, J, ecoords, ri_Ns[j], ri_gradNparams[j])
                Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], ri_Ns[j])
                _e_g!(E_G, J)
                _gradN_e!(gradN_e, J, E_G, ri_gradNparams[j])
                t = self.integdomain.otherdimension(loc, fes.conn[i], ri_Ns[j])
                _nodal_triads_e!(A_Es, nvalid, E_G, normals, normal_valid, fes.conn[i])
                _transfmat_g_to_a!(Tga, A_Es, E_G)
                _transfmat_a_to_e!(Tae, A_Es, gradN_e)
                mul!(T, Tae, Tga)
                updatecsmat!(outputcsys, loc, J, i, j)
                if dot(view(csmat(outputcsys), :, 3), view(E_G, :, 3)) < 0.95
                    @warn "Coordinate systems mismatched?"
                end
                ocsm, ocsn = lla(E_G, csmat(outputcsys))
                o2_e[1, 1] = o2_e[2, 2] = ocsm
                o2_e[1, 2] = ocsn
                o2_e[2, 1] = -ocsn
                bsmat!(Bs, gradN_e, ri_Ns[j], T)
                shr = Bs * edisp
                frc = (t^3 / (t^2 + mult_el_size * Jac * ri_w[j])) * Dt * shr
                fo = o2_e' * frc
                out[1:2] .= fo[1], fo[2]
                # Call the inspector
                idat = inspector(idat, i, fes.conn[i], ecoords, out, loc)
            end # Loop over quadrature points
        else # Bending moment or membrane force
            for j in 1:fi_npts # Loop over quadrature points
                locjac!(loc, J, ecoords, fi_Ns[j], fi_gradNparams[j])
                _e_g!(E_G, J)
                _gradN_e!(gradN_e, J, E_G, fi_gradNparams[j])
                t = self.integdomain.otherdimension(loc, fes.conn[i], fi_Ns[j])
                _nodal_triads_e!(A_Es, nvalid, E_G, normals, normal_valid, fes.conn[i])
                _transfmat_g_to_a!(Tga, A_Es, E_G)
                _transfmat_a_to_e!(Tae, A_Es, gradN_e)
                mul!(T, Tae, Tga)
                updatecsmat!(outputcsys, loc, J, i, j)
                if dot(view(csmat(outputcsys), :, 3), view(E_G, :, 3)) < 0.95
                    @warn "Coordinate systems mismatched?"
                end
                ocsm, ocsn = lla(E_G, csmat(outputcsys))
                o2_e[1, 1] = o2_e[2, 2] = ocsm
                o2_e[1, 2] = ocsn
                o2_e[2, 1] = -ocsn
                # Compute the Requested Quantity
                if quant == BENDING_MOMENT
                    bbmat!(Bb, gradN_e, T)
                    kurv = Bb * edisp
                    mom = ((t^3) / 12) * Dps * kurv
                    m = [mom[1] mom[3]; mom[3] mom[2]]
                    mo = o2_e' * m * o2_e
                    out[:] .= mo[1, 1], mo[2, 2], mo[1, 2]
                end
                if quant == MEMBRANE_FORCE
                    bmmat!(Bm, gradN_e, T)
                    strn = Bm * edisp
                    frc = (t) * Dps * strn
                    f = [frc[1] frc[3]; frc[3] frc[2]]
                    fo = o2_e' * f * o2_e
                    out[:] .= fo[1, 1], fo[2, 2], fo[1, 2]
                end
                # Call the inspector
                idat = inspector(idat, i, fes.conn[i], ecoords, out, loc)
            end # Loop over quadrature points
        end # select quantity
    end # Loop over elements
    return idat # return the updated inspector data
end

function inspectintegpoints(
    self::FEMMShellQ4RNT,
    geom::NodalField{FFlt},
    u::NodalField{TI},
    felist::FIntVec,
    inspector::F,
    idat,
    quantity = :Cauchy;
    context...,
) where {TI<:Number,F<:Function}
    dT = NodalField(fill(zero(FFlt), nnodes(geom), 1)) # zero difference in temperature
    return inspectintegpoints(
        self,
        geom,
        u,
        dT,
        felist,
        inspector,
        idat,
        quantity;
        context...,
    )
end


end # module
