module FEMMShellT3FFCompModule

using LinearAlgebra: norm, Transpose, mul!, diag, eigen, I, dot, rank
using Statistics: mean
using FinEtools
import FinEtools.FESetModule: gradN!, nodesperelem, manifdim
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FEMMBaseModule: associategeometry!
import FinEtools.FEMMBaseModule: inspectintegpoints
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!, update!, thermalstrain!
using FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, locjac!, add_nnt_ut_only!, add_btsigma!, add_b1tdb2!
using ..FESetShellT3Module: FESetShellT3
using ..TransformerModule: QTEQTransformer, QEQTTransformer, Layup2ElementAngle
using ..CompositeLayupModule: CompositeLayup, thickness, laminate_stiffnesses!, laminate_transverse_stiffness!, plane_stress_T_matrix!, transverse_shear_T_matrix!, laminate_inertia!


const __nn = 3 # number of nodes
const __ndof = 6 # number of degrees of freedom per node

# Formulation for the transverse shear stiffness which averages the
# strain-displacement matrix.
const __TRANSV_SHEAR_FORMULATION_AVERAGE_B = 0
# Formulation for the transverse shear stiffness which averages the
# stiffness matrix.
const __TRANSV_SHEAR_FORMULATION_AVERAGE_K = 1


"""
    FEMMShellT3FFComp{S<:FESetT3, F<:Function} <: AbstractFEMM

Class for the finite element modeling machine of the T3 triangular Flat-Facet
shell with the Discrete Shear Gap technology and a consistent handling of the
normals. This formulation is suitable for modelling of COMPOSITE (layered) materials. 

For details for the homogeneous-shell refer to [`FEMMShellT3FF`](@ref).
"""
mutable struct FEMMShellT3FFComp{S<:FESetT3, F<:Function} <: AbstractFEMM
    integdomain::IntegDomain{S, F} # integration domain data
    # Definitions of layups
    layup_groups::Vector{Tuple{CompositeLayup, Vector{FInt}}} # layups: vector of pairs of the composite layup and the group of elements using that layup.
    # Configuration parameters
    transv_shear_formulation::FInt
    drilling_stiffness_scale::Float64
    threshold_angle::Float64
    mult_el_size::Float64
    # Private data
    _associatedgeometry::Bool
    _normals::FFltMat
    _normal_valid::Vector{Bool}
    _layup_group_lookup::FIntVec # look up the layup group for each element
    # The attributes below are buffers used in various operations.
    _loc::FFltMat
    _J0::FFltMat
    _ecoords::FFltMat
    _edisp::FFltVec
    _ecoords_e::FFltMat
    _edisp_e::FFltVec
    _dofnums::FIntMat
    _E_G::FFltMat
    _A_Es::Vector{FFltMat} # transformation nodal-element matrices  
    _nvalid::Vector{Bool}
    _T::FFltMat  # element transformation matrix
    _elmat::FFltMat
    _gradN_e::FFltMat
    _Bm::FFltMat
    _Bb::FFltMat
    _Bs::FFltMat
    _DpsBmb::FFltMat
    _DtBs::FFltMat
end


"""
    FEMMShellT3FFComp(integdomain::IntegDomain{S, F}, layup::CompositeLayup) where {S<:FESetT3, F<:Function}

Constructor of the T3FFComp shell FEMM. All elements use a single layup.
"""
function FEMMShellT3FFComp(integdomain::IntegDomain{S,F}, layup::CompositeLayup) where {S<:FESetT3,F<:Function}
    _nnmax = 0
    for j in 1:count(integdomain.fes)
        for k in eachindex(integdomain.fes.conn[j])
            _nnmax = max(_nnmax, integdomain.fes.conn[j][k])
        end
    end
    # When there is only a single layup, all elements belong to a single group
    layup_groups = [(layup, collect(1:count(integdomain.fes)))]
    # Establish the mapping from the elements to the layup group
    _layup_group_lookup = fill(zero(FInt), count(integdomain.fes))
    for j in 1:length(layup_groups)
        eset = layup_groups[j][2]
        for i in eset # Loop over elements in the layup group
            _layup_group_lookup[i] = j
        end
    end
    # Allocate space for the normals
    _normals = fill(0.0, _nnmax, 3)
    _normal_valid = fill(true, _nnmax)
    # Alocate the buffers
    _loc = fill(0.0, 1, 3)
    _J0 = fill(0.0, 3, 2)
    _ecoords = fill(0.0, __nn, 3)
    _edisp = fill(0.0, __nn * __ndof)
    _ecoords_e = fill(0.0, __nn, 2)
    _edisp_e = fill(0.0, __nn * __ndof)
    _dofnums = zeros(FInt, 1, __nn * __ndof)
    _E_G = fill(0.0, 3, 3)
    _A_Es = [fill(0.0, 3, 3), fill(0.0, 3, 3), fill(0.0, 3, 3)]
    _nvalid = fill(false, 3)
    _T = fill(0.0, __nn * __ndof, __nn * __ndof)
    _elmat = fill(0.0, __nn * __ndof, __nn * __ndof)
    _gradN_e = fill(0.0, __nn, 2)
    _Bm = fill(0.0, 3, __nn * __ndof)
    _Bb = fill(0.0, 3, __nn * __ndof)
    _Bs = fill(0.0, 2, __nn * __ndof)
    _DpsBmb = similar(_Bm)
    _DtBs = similar(_Bs)

    @assert delegateof(integdomain.fes) === FESetShellT3()

    return FEMMShellT3FFComp(integdomain, layup_groups,
        __TRANSV_SHEAR_FORMULATION_AVERAGE_B, 1.0, 30.0, 5 / 12 / 1.5,
        false, _normals, _normal_valid,
        _layup_group_lookup,
        _loc, _J0,
        _ecoords, _edisp, _ecoords_e, _edisp_e,
        _dofnums,
        _E_G, _A_Es, _nvalid, _T,
        _elmat,
        _gradN_e,
        _Bm, _Bb, _Bs, _DpsBmb, _DtBs)
end

"""
    make(integdomain, layup)

Make a T3FFComp FEMM from the integration domain and a composite layup.
"""
function make(integdomain, layup)
    return FEMMShellT3FFComp(integdomain, layup)
end

function _compute_nodal_normal!(n, mcsys::CSys, XYZ, J0::FFltMat, labl::FInt)
    updatecsmat!(mcsys, reshape(XYZ, 1, 3), J0, labl);
    n[:] .= mcsys.csmat[:, 3]
    return n
end

function _compute_J!(J0, ecoords)
    # Compute the Jacobian matrix: the Jacobian matrix is constant within the
    # triangle
    # x, y, z = ecoords[2, :].-ecoords[1, :]
    # J0[:, 1] .= (x, y, z)
    # x, y, z = ecoords[3, :].-ecoords[1, :]
    # J0[:, 2] .= (x, y, z)
    for j in 1:3 
        J0[j, 1] = ecoords[2, j] - ecoords[1, j]
        J0[j, 2] = ecoords[3, j] - ecoords[1, j]  
    end
    return J0
end

function _e_g!(E_G, J0)
    # J0[:, 1] is the tangent to the coordinate curve 1
    for j in 1:3
        E_G[j,1]  = J0[j, 1]
    end
    n = sqrt(E_G[1, 1]^2 + E_G[2, 1]^2 + E_G[3, 1]^2)
    for j in 1:3
        E_G[j,1]  /= n
    end
    # J0[:, 2] is the tangent to the coordinate curve 2
    # Now compute the normal
    E_G[1, 3] = -E_G[3, 1]*J0[2, 2]+E_G[2, 1]*J0[3, 2]
    E_G[2, 3] =  E_G[3, 1]*J0[1, 2]-E_G[1, 1]*J0[3, 2]
    E_G[3, 3] = -E_G[2, 1]*J0[1, 2]+E_G[1, 1]*J0[2, 2]
    n = sqrt(E_G[1, 3]^2 + E_G[2, 3]^2 + E_G[3, 3]^2);
    for j in 1:3
        E_G[j, 3]  /= n
    end
    E_G[1, 2] = -E_G[3, 3]*E_G[2, 1]+E_G[2, 3]*E_G[3, 1]
    E_G[2, 2] =  E_G[3, 3]*E_G[1, 1]-E_G[1, 3]*E_G[3, 1]
    E_G[3, 2] = -E_G[2, 3]*E_G[1, 1]+E_G[1, 3]*E_G[2, 1]
    return E_G
end

function _ecoords_e!(ecoords_e, J0, E_G)
    ecoords_e[1, 1] = ecoords_e[1, 2] = 0.0
    ecoords_e[2, 1] = dot(view(J0, :, 1), view(E_G, :, 1))
    ecoords_e[2, 2] = dot(view(J0, :, 1), view(E_G, :, 2))
    ecoords_e[3, 1] = dot(view(J0, :, 2), view(E_G, :, 1))
    ecoords_e[3, 2] = dot(view(J0, :, 2), view(E_G, :, 2))
    return ecoords_e
end

function _centroid!(centroid, ecoords)
    centroid[1]  = (ecoords[1, 1]+ecoords[2, 1]+ecoords[3, 1])/3
    centroid[2]  = (ecoords[1, 2]+ecoords[2, 2]+ecoords[3, 2])/3
    centroid[3]  = (ecoords[1, 3]+ecoords[2, 3]+ecoords[3, 3])/3
    return centroid
end
   
struct NodalTriadsE
    r::FFltVec
    nk_e::FFltVec
    nk::FFltVec
    f3_e::FFltVec
end

function NodalTriadsE()
    NodalTriadsE(fill(0.0, 3), fill(0.0, 3), fill(0.0, 3), [0.0, 0.0, 1.0])
end
   
(o::NodalTriadsE)(A_Es, nvalid, E_G, normals, normal_valid, c) = begin
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

    for k in 1:length(A_Es)
        o.nk .= vec(view(normals, c[k], :))
        nvalid[k] = normal_valid[c[k]]
        if nvalid[k] 
            # If the nodal normal is valid, pull it back into the element frame.
            mul!(o.nk_e, transpose(E_G), o.nk)
        else
            # Otherwise the element normal replaces the nodal normal.
            o.nk_e .= 0.0; o.nk_e[3] = 1.0
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

struct TransfmatGToA
    Tblock::FFltMat
end

function TransfmatGToA()
    TransfmatGToA(fill(0.0, 3, 3))
end

(o::TransfmatGToA)(T, A_Es, E_G) = begin
    # Global-to-nodal transformation matrix. 

    # The 3x3 blocks consist of the nodal triad expressed on the global basis.
    # The nodal basis vectors in the array `A_Es[i]` are expressed on the
    # element basis, and are then rotated with `E_G` into the global coordinate
    # system.

    # Output
    # - `T` = transformation matrix, input in the global basis, output in the
    #   nodal basis
    T .= 0.0
    for i in 1:__nn
        mul!(o.Tblock,  transpose(A_Es[i]), transpose(E_G))
        offset = (i-1)*__ndof
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
    
    T .= 0.0
    for i in 1:__nn
        roffst = (i-1)*__ndof
        iA_E = A_Es[i] 
        iA_E_33 = iA_E[3, 3]
        # T[r, r] .= iA_E
        for cl in 1:3
            for rw in 1:3
                T[roffst+rw, roffst+cl] = iA_E[rw, cl]
            end
        end
        for cl in 1:2
            for rw in 1:2
                T[roffst+3+rw, roffst+3+cl] = iA_E[rw, cl] - (1/iA_E_33)*iA_E[rw, 3]*iA_E[cl, 3]
            end
        end
        m1 = (1/iA_E_33)*iA_E[1, 3]
        m2 = (1/iA_E_33)*iA_E[2, 3]
        for j in 1:__nn
            coffst = (j-1)*__ndof
            for k in 1:3
                a3 = 1/2 * (iA_E[2, k] * gradN_e[j, 1] - iA_E[1, k] * gradN_e[j, 2])
                T[roffst+4, coffst+k] += m1 * a3
                T[roffst+5, coffst+k] += m2 * a3
            end
        end
    end
    return T
end

function _gradN_e_Ae!(gradN_e, ecoords_e)
    # Compute the gradient with respect to the element coordinates, three
    # gradients in two coordinates, and the area of the triangle
    # a, b = ecoords_e[2, :] .- ecoords_e[1, :]
    # c, d = ecoords_e[3, :] .- ecoords_e[1, :]
    # J = (a*d - b*c)
    # gradN_e[:] .= ((b-d)/J, d/J, -b/J, (c-a)/J, -c/J, a/J)
     a = ecoords_e[2, 1] - ecoords_e[1, 1]
     b = ecoords_e[2, 2] - ecoords_e[1, 2]
     c = ecoords_e[3, 1] - ecoords_e[1, 1]
     d = ecoords_e[3, 2] - ecoords_e[1, 2]
     J = (a*d - b*c)
     gradN_e[1, 1] = (b-d)/J
     gradN_e[2, 1] = d/J
     gradN_e[3, 1] = -b/J
     gradN_e[1, 2] = (c-a)/J
     gradN_e[2, 2] = -c/J
     gradN_e[3, 2] = a/J
    return gradN_e, J/2
end

function _add_Bsmat_o!(Bs, ecoords_e, Ae, ordering)
    # Compute the linear transverse shear strain-displacement matrix for one
    # particular ordering of the nodes. The computed entries are ADDED. The
    # matrix is NOT zeroed out initially. That is the responsibility of the
    # caller.

    # Orientation 
    s, p, q = ordering
    a = ecoords_e[p, 1] - ecoords_e[s, 1]
    b = ecoords_e[p, 2] - ecoords_e[s, 2]
    c = ecoords_e[q, 1] - ecoords_e[s, 1]
    d = ecoords_e[q, 2] - ecoords_e[s, 2]
    m = (1/2/Ae) # multiplier

    # The first node in the triangle 
    # Node s
    co = (s - 1) * 6 # column offset
    Bs[1, co+3] += m*(b-d);                             Bs[1, co+5] += m*(Ae) 
    Bs[2, co+3] += m*(c-a); Bs[2, co+4] += m*(-Ae); 
    # The other two nodes
    # Node p
    co = (p - 1) * 6 # column offset
    Bs[1, co+3] += m*(d);   Bs[1, co+4] += m*(-b*d/2);  Bs[1, co+5] += m*(a*d/2) 
    Bs[2, co+3] += m*(-c);  Bs[2, co+4] += m*(b*c/2);   Bs[2, co+5] += m*(-a*c/2) 
    # Node q
    co = (q - 1) * 6 # column offset
    Bs[1, co+3] += m*(-b);  Bs[1, co+4] += m*(b*d/2);   Bs[1, co+5] += m*(-b*c/2) 
    Bs[2, co+3] += m*(a);   Bs[2, co+4] += m*(-a*d/2);  Bs[2, co+5] += m*(a*c/2) 
    
    return Bs
end

function _Bsmat!(Bs, ecoords_e, Ae)
    # Compute the linear transverse shear strain-displacement matrix.
    Bs .= 0.0 # Zero out initially, then add the three contributions  
    _add_Bsmat_o!(Bs, ecoords_e, Ae, (1, 2, 3))
    _add_Bsmat_o!(Bs, ecoords_e, Ae, (2, 3, 1))
    _add_Bsmat_o!(Bs, ecoords_e, Ae, (3, 1, 2))
    # Now averaging of the three contributions
    Bs .*= (1/3)
    return Bs
end

function _Bmmat!(Bm, gradN)
    # Compute the linear membrane strain-displacement matrix.
    fill!(Bm, 0.0)
    for i in 1:__nn
        Bm[1,6*(i-1)+1] = gradN[i,1];
                                             Bm[2,6*(i-1)+2] = gradN[i,2];
        Bm[3,6*(i-1)+1] = gradN[i,2];        Bm[3,6*(i-1)+2] = gradN[i,1];
    end
end

function _Bbmat!(Bb, gradN)
    # Compute the linear, displacement independent, curvature-displacement/rotation
    # matrix for a shell quadrilateral element with nfens=3 nodes. Displacements and
    # rotations are in a local coordinate system.
    fill!(Bb, 0.0)
    for i in 1:__nn
                                              Bb[1,6*(i-1)+5] = gradN[i,1];
        Bb[2,6*(i-1)+4] = -gradN[i,2];
        Bb[3,6*(i-1)+4] = -gradN[i,1];        Bb[3,6*(i-1)+5] = gradN[i,2];
    end
end

"""
    associategeometry!(self::FEMMShellT3FFComp,  geom::NodalField{FFlt})

Associate geometry with the FEMM. 

In this case it means evaluate the nodal normals.
"""
function associategeometry!(self::FEMMShellT3FFComp,  geom::NodalField{FFlt})
    threshold_angle = self.threshold_angle
    J0 = self._J0
    E_G = self._E_G
    normals, normal_valid = self._normals, self._normal_valid
    nnormal = fill(0.0, 3)

    # Compute the normals at the nodes
    for lg in self.layup_groups
        layup = lg[1]
        eset = lg[2]
        for el in eset
            i, j, k = self.integdomain.fes.conn[el]
            J0[:, 1] = geom.values[j, :] - geom.values[i, :]
            J0[:, 2] = geom.values[k, :] - geom.values[i, :]
            for n in [i, j, k]
                _compute_nodal_normal!(nnormal, layup.csys, geom.values[n, :], J0, self.integdomain.fes.label[el])
                normals[n, :] .+= nnormal
            end
        end
    end
    
    # Normalize to unit length
    for j in 1:size(normals, 1)
        nn = norm(normals[j, :])
        if nn > 0.0
            normals[j, :] ./= nn
        end
    end

    # Now perform a second pass. If the nodal normal differs by a substantial
    # amount from the element normals, zero out the nodal normal to indicate
    # that the nodal normal should not be used at that vertex. 
    ntolerance = 1 - sqrt(1 - sin(threshold_angle/180*pi)^2)
    for el in 1:count(self.integdomain.fes)
        i, j, k = self.integdomain.fes.conn[el]
        _compute_J!(J0, geom.values[[i, j, k], :])
        for n in [i, j, k]
            _e_g!(E_G, J0)
            nd = dot(normals[n, :], E_G[:, 3])
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
    num_normals(self::FEMMShellT3FFComp)

Compute the summary of the nodal normals.
"""
function num_normals(self::FEMMShellT3FFComp)
    normals, normal_valid = self._normals, self._normal_valid
    total_normals = length(normal_valid)
    total_invalid_normals = length([v for v in normal_valid if v != true])
    return total_normals, total_invalid_normals
end

"""
    stiffness(self::FEMMShellT3FFComp, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{TI}, Rfield1::NodalField{TI}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, TI<:Number}

Compute the material stiffness matrix.
"""
function stiffness(self::FEMMShellT3FFComp, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{TI}, Rfield1::NodalField{TI}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, TI<:Number}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    label = self.integdomain.fes.label
    normals, normal_valid = self._normals, self._normal_valid
    centroid, J0 = self._loc, self._J0
    ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
    ecoords_e, gradN_e = self._ecoords_e, self._gradN_e 
    E_G, A_Es, nvalid, T = self._E_G, self._A_Es, self._nvalid, self._T
    elmat = self._elmat
    transformwith = QTEQTransformer(elmat)
    lla = Layup2ElementAngle()
    _nodal_triads_e! = NodalTriadsE()
    _transfmat_g_to_a! = TransfmatGToA()
    Bm, Bb, Bs, DpsBmb, DtBs = self._Bm, self._Bb, self._Bs, self._DpsBmb, self._DtBs
    A, B, C, BT = zeros(3, 3), zeros(3, 3), zeros(3, 3), zeros(3, 3)
    H = zeros(2, 2)
    sA, sB, sC = zeros(3, 3), zeros(3, 3), zeros(3, 3)
    sH = zeros(2, 2)
    Tps, Tts = zeros(3, 3), zeros(2, 2)
    tps! = QEQTTransformer(Tps)
    tts! = QEQTTransformer(Tts)
    drilling_stiffness_scale = self.drilling_stiffness_scale
    transv_shear_formulation = self.transv_shear_formulation
    mult_el_size = self.mult_el_size
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), dchi.nfreedofs, dchi.nfreedofs);
    for lg in self.layup_groups
        layup = lg[1]
        eset = lg[2]
        t = thickness(layup)
        laminate_stiffnesses!(layup, A, B, C)
        laminate_transverse_stiffness!(layup, H)
        for i in eset # Loop over elements in the layup group
            gathervalues_asmat!(geom0, ecoords, fes.conn[i]);
            _centroid!(centroid, ecoords)
            _compute_J!(J0, ecoords)
            _e_g!(E_G, J0)
            _ecoords_e!(ecoords_e, J0, E_G)
            gradN_e, Ae = _gradN_e_Ae!(gradN_e, ecoords_e)
            # Working copies to be transformed
            sA[:] .= A[:]; sB[:] .= B[:]; sC[:] .= C[:];     sH[:] .= H[:]
            # Transform the laminate stiffnesses
            updatecsmat!(layup.csys, reshape(centroid, 1, 3), J0, -1);
            m, n = lla(E_G, layup.csys.csmat) 
            plane_stress_T_matrix!(Tps, m, -n)
            tps!(sA, Tps); tps!(sB, Tps); tps!(sC, Tps); 
            transverse_shear_T_matrix!(Tts, m, n)
            tts!(sH, Tts)
            # Construct the Stiffness Matrix
            fill!(elmat,  0.0); # Initialize element matrix
            _Bmmat!(Bm, gradN_e)
            add_btdb_ut_only!(elmat, Bm, Ae, sA, DpsBmb)
            _Bbmat!(Bb, gradN_e)
            add_btdb_ut_only!(elmat, Bb, Ae, sC, DpsBmb)
            add_b1tdb2!(elmat, Bm, Bb, Ae, sB, DpsBmb)
            add_b1tdb2!(elmat, Bb, Bm, Ae, sB, DpsBmb)
            # he = sqrt(2*Ae) # we avoid taking the square root here, replacing
            # he^2 with 2*Ae
            if transv_shear_formulation == __TRANSV_SHEAR_FORMULATION_AVERAGE_K
                for  o in [(1, 2, 3), (2, 3, 1), (3, 1, 2)]
                    Bs .= 0.0; _add_Bsmat_o!(Bs, ecoords_e, Ae, o)
                    add_btdb_ut_only!(elmat, Bs, (t^2/(t^2+mult_el_size*2*Ae))*Ae/3, sH, DtBs)
                end
            else
                _Bsmat!(Bs, ecoords_e, Ae)
                add_btdb_ut_only!(elmat, Bs, (t^2/(t^2+mult_el_size*2*Ae))*Ae, sH, DtBs)
            end
            # Complete the elementwise matrix by filling in the lower triangle
            complete_lt!(elmat)
            # Transformation from the nodal (A) to the element (E) basis
            _nodal_triads_e!(A_Es, nvalid, E_G, normals, normal_valid, fes.conn[i])
            _transfmat_a_to_e!(T, A_Es, gradN_e)
            transformwith(elmat, T)
            # Bending diagonal stiffness coefficients
            kavg = mean((elmat[4, 4], elmat[10, 10], elmat[16, 16],
                elmat[5, 5], elmat[11, 11], elmat[17, 17])) * drilling_stiffness_scale
            # Add the artificial drilling stiffness in the nodal basis, but only
            # if the nodal normal is valid. 
            nvalid[1] && (elmat[6,6] += kavg)
            nvalid[2] && (elmat[12,12] += kavg)
            nvalid[3] && (elmat[18,18] += kavg)
            # Transform from global (G) into the nodal (A) basis
            _transfmat_g_to_a!(T, A_Es, E_G)
            transformwith(elmat, T)
            # Assembly
            gatherdofnums!(dchi, dofnums, fes.conn[i]); 
            assemble!(assembler, elmat, dofnums, dofnums); 
        end # Loop over elements
    end # Loop over layup groups
    return makematrix!(assembler);
end

function stiffness(self::FEMMShellT3FFComp, geom0::NodalField{FFlt}, u1::NodalField{TI}, Rfield1::NodalField{TI}, dchi::NodalField{TI}) where {TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi);
end


"""
    mass(self::FEMMShellT3FFComp,  assembler::A,  geom::NodalField{FFlt}, dchi::NodalField{TI}) where {A<:AbstractSysmatAssembler, TI<:Number}

Compute the diagonal (lumped) mass matrix

The mass matrix can be expected to be non-singular.
"""
function mass(self::FEMMShellT3FFComp,  assembler::A,  geom0::NodalField{FFlt}, dchi::NodalField{TI}) where {A<:AbstractSysmatAssembler, TI<:Number}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    label = self.integdomain.fes.label
    normals, normal_valid = self._normals, self._normal_valid
    centroid, J0 = self._loc, self._J0
    ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
    ecoords_e, gradN_e = self._ecoords_e, self._gradN_e 
    E_G, A_Es, nvalid, T = self._E_G, self._A_Es, self._nvalid, self._T
    elmat = self._elmat
    transformwith = QTEQTransformer(elmat)
    _nodal_triads_e! = NodalTriadsE()
    _transfmat_g_to_a! = TransfmatGToA()
    npe = nodesperelem(fes)
    ndn = ndofs(dchi)
    drilling_stiffness_scale = self.drilling_stiffness_scale
    startassembly!(assembler,  size(elmat,1),  size(elmat,2),  count(fes), dchi.nfreedofs,  dchi.nfreedofs);
    for lg in self.layup_groups
            layup = lg[1]
            eset = lg[2]
            mass_density, moment_of_inertia_density = laminate_inertia!(layup)
            for i in eset # Loop over elements in the layup group
                gathervalues_asmat!(geom0, ecoords, fes.conn[i]);
                _centroid!(centroid, ecoords)
                _compute_J!(J0, ecoords)
                _e_g!(E_G, J0)
                _ecoords_e!(ecoords_e, J0, E_G)
                gradN_e, Ae = _gradN_e_Ae!(gradN_e, ecoords_e)
                # Compute the translational and rotational masses corresponding
                # to nodes

                tfactor = mass_density*(Ae/3);
                rfactor = moment_of_inertia_density*(Ae/3);
                # end # Loop over quadrature points

                # Now treat the transformation from the element to the nodal
                # triad
                _nodal_triads_e!(A_Es, nvalid, E_G, normals, normal_valid, fes.conn[i])
                # Fill the elementwise matrix in the nodal basis: the matrix is
                # lumped, hence diagonal, and no transformation from element to
                # nodal is needed.
                fill!(elmat,  0.0); # Initialize element matrix
                for k in 1:npe
                    # Translation degrees of freedom
                    for d in 1:3
                        c = (k - 1) * __ndof + d
                        elmat[c, c] += tfactor
                    end
                    # Bending degrees of freedom
                    for d in 4:5
                        c = (k - 1) * __ndof + d
                        elmat[c, c] += rfactor
                    end
                    # Drilling rotations
                    if nvalid[k]
                        d = 6
                        c = (k - 1) * __ndof + d
                        elmat[c, c] += rfactor * drilling_stiffness_scale / 1e3
                    end
                end
                # Transform into global coordinates
                _transfmat_g_to_a!(T, A_Es, E_G)
                transformwith(elmat, T)
                # Assemble
                gatherdofnums!(dchi,  dofnums,  fes.conn[i]);# retrieve degrees of freedom
                assemble!(assembler,  elmat,  dofnums,  dofnums);# assemble symmetric matrix
            end # Loop over elements
        end # Loop over layup groups
    return makematrix!(assembler);
end

function mass(self::FEMMShellT3FFComp,  geom::NodalField{FFlt},  u::NodalField{TI}) where {TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return mass(self, assembler, geom, u);
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
function inspectintegpoints(self::FEMMShellT3FFComp, geom0::NodalField{FFlt},  u::NodalField{TI}, dT::NodalField{FFlt}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {TI<:Number, F<:Function}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    label = self.integdomain.fes.label
    normals, normal_valid = self._normals, self._normal_valid
    centroid, J0 = self._loc, self._J0
    ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
    ecoords_e, gradN_e = self._ecoords_e, self._gradN_e 
    E_G, A_Es, nvalid, T = self._E_G, self._A_Es, self._nvalid, self._T
    elmat = self._elmat
    transformwith = QTEQTransformer(elmat)
    _nodal_triads_e! = NodalTriadsE()
    _transfmat_g_to_a! = TransfmatGToA()
    Bm, Bb, Bs, DpsBmb, DtBs = self._Bm, self._Bb, self._Bs, self._DpsBmb, self._DtBs
    A, B, C, BT = zeros(3, 3), zeros(3, 3), zeros(3, 3), zeros(3, 3)
    H = zeros(2, 2)
    sA, sB, sC = zeros(3, 3), zeros(3, 3), zeros(3, 3)
    sH = zeros(2, 2)
    Tps, Tts = zeros(3, 3), zeros(2, 2)
    tps! = QEQTTransformer(Tps)
    tts! = QEQTTransformer(Tts)
    lla = Layup2ElementAngle()
    mult_el_size = self.mult_el_size
    edisp_e = deepcopy(edisp)
    edisp_n = deepcopy(edisp_e)
    out = fill(0.0, 3)
    o2_e = fill(0.0, 2, 2)
    # Sort out  the output requirements
    outputcsys = self.layup_groups[1][1].csys; # default: report the stresses in the material coord system
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
    for ilist = 1:length(felist) # Loop over elements
        i = felist[ilist];
        # Look up the layup
        layup = self.layup_groups[self._layup_group_lookup[i]][1]
        t = thickness(layup)
        laminate_stiffnesses!(layup, A, B, C)
        laminate_transverse_stiffness!(layup, H)
        gathervalues_asmat!(geom0, ecoords, fes.conn[i]);
        gathervalues_asvec!(u, edisp, fes.conn[i]);
        _centroid!(centroid, ecoords)
        _compute_J!(J0, ecoords)
        _e_g!(E_G, J0)
        _ecoords_e!(ecoords_e, J0, E_G)
        gradN_e, Ae = _gradN_e_Ae!(gradN_e, ecoords_e)
        t = self.integdomain.otherdimension(centroid, fes.conn[i], [1.0/3 1.0/3])
        # Establish nodal triads
        _nodal_triads_e!(A_Es, nvalid, E_G, normals, normal_valid, fes.conn[i])
        # Transform from global into nodal coordinates
        _transfmat_g_to_a!(T, A_Es, E_G)
        mul!(edisp_n, T, edisp)
        # Now treat the transformation from the nodal to the element triad
        _transfmat_a_to_e!(T, A_Es, gradN_e)
         # Transform the nodal vector into the elementwise coordinates
        mul!(edisp_e, T, edisp_n)
        updatecsmat!(outputcsys, centroid, J0, fes.label[i]);
        if dot(view(outputcsys.csmat, :, 3), view(E_G, :, 3)) < 0.95
            @warn "Coordinate systems mismatched?"
        end
        # Established the stiffness matrices
        sA[:] .= A[:]; sB[:] .= B[:]; sC[:] .= C[:];     sH[:] .= H[:]
        # Transform the laminate stiffnesses
        updatecsmat!(layup.csys, centroid, J0, fes.label[i]);
        m, n = lla(E_G, layup.csys.csmat)
        plane_stress_T_matrix!(Tps, m, -n)
        tps!(sA, Tps); tps!(sB, Tps); tps!(sC, Tps); 
        transverse_shear_T_matrix!(Tts, m, n)
        tts!(sH, Tts)
        # The output coordinate system
        ocsm, ocsn = lla(E_G, outputcsys.csmat)
        o2_e[1, 1] = o2_e[2, 2]  = ocsm
        o2_e[1, 2] = ocsn; o2_e[2, 1] = -ocsn
        # Compute the Requested Quantity
        _Bbmat!(Bb, gradN_e)
        _Bmmat!(Bm, gradN_e)
        kurv = Bb * edisp_e
        memstr = Bm * edisp_e
        if quant == BENDING_MOMENT
            mom = sB * memstr + sC * kurv
            m = [mom[1] mom[3]; mom[3] mom[2]]
            mo = o2_e' * m * o2_e
            out[:] .= mo[1, 1], mo[2, 2], mo[1, 2]
        end
        if quant == MEMBRANE_FORCE
            frc = sA * memstr + sB * kurv
            f = [frc[1] frc[3]; frc[3] frc[2]]
            fo = o2_e' * f * o2_e
            out[:] .= fo[1, 1], fo[2, 2], fo[1, 2]
            # @infiltrate
        end
        if quant == TRANSVERSE_SHEAR
            _Bsmat!(Bs, ecoords_e, Ae)
            he = sqrt(2*Ae)
            shrstr = Bs * edisp_e
            frc = (t^2/(t^2+mult_el_size*2*Ae))*sH * shrstr
            fo = o2_e' * frc
            out[1:2] .= fo[1], fo[2]
        end 
        # Call the inspector
        idat = inspector(idat, i, fes.conn[i], ecoords, out, centroid);
        # end # Loop over quadrature points
    end # Loop over elements
    return idat; # return the updated inspector data
end

function inspectintegpoints(self::FEMMShellT3FFComp, geom::NodalField{FFlt},  u::NodalField{TI}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {TI<:Number, F<:Function}
    dT = NodalField(fill(zero(FFlt), nnodes(geom), 1)) # zero difference in temperature
    return inspectintegpoints(self, geom, u, dT, felist, inspector, idat, quantity; context...);
end

end # module
