module FEMMShellT3FFModule

using LinearAlgebra: norm, Transpose, mul!, diag, eigen, I, dot, rank
using Statistics: mean
using FinEtools
import FinEtools.FESetModule: gradN!, nodesperelem, manifdim
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FEMMBaseModule: associategeometry!
import FinEtools.FEMMBaseModule: inspectintegpoints
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!, update!, thermalstrain!
using FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, locjac!, add_nnt_ut_only!, add_btsigma!
using ..FESetShellT3Module: FESetShellT3
using ..TransformerModule: Transformer

using Infiltrator

const __nn = 3 # number of nodes
const __ndof = 6 # number of degrees of freedom per node

# Formulation for the transverse shear stiffness which averages the
# strain-displacement matrix.
const __TRANSV_SHEAR_FORMULATION_AVERAGE_B = 0
# Formulation for the transverse shear stiffness which averages the
# stiffness matrix.
const __TRANSV_SHEAR_FORMULATION_AVERAGE_K = 1


"""
    FEMMShellT3FF{S<:AbstractFESet, F<:Function} <: AbstractFEMM

Class for the finite element modeling machine of the Flat-Facet shell with the
Discrete Shear Gap technology and a consistent handling of the normals. 

With averaging of the transverse strain-displacement matrix or averaging of the
transverse shear stiffness matrix to provide isotropic transverse shear
response. 

Also, the formulation is developed to correctly handle the coupling of twisting
moments and transverse shear (such as in the twisted beam or the Raasch hook
problems) by incorporating "nodal" normals.

Some of the programming developed consistently with the paper

[1] Cui et al, Analysis of plates and shells using an edge-based smoothed finite
element method, Comput Mech (2010) 45:141–156 DOI 10.1007/s00466-009-0429-9

The stabilization factor of the shear term of

[2] Mikko Lyly, Rolf Stenberg and Teemu Vihinen, A stable bilinear element for
the Reissner-Mindlin plate model Computer Methods in Applied Mechanics and
Engineering 110 (1993) 343-357 

is incorporated. Refer to expressions (3.12) and (3.13).

The treatment of the transformation between the element and nodal coordinates
is carried out using a clean alternative to the publication

[3] Finite Elements in Analysis and Design 30 (1998) 235—242
The treatment of shell normals in ﬁnite element analysis
Richard H. MacNeal, Charles T. Wilson, Robert L. Harder, Claus C. Hoﬀ
The MacNeal-Schwendler Corporation, 815 Colorado Blvd., Los Angeles, CA 90041, USA

The formula for the element to nodal basis transformation is derived from the
expression

    [theta]_n = [A]_E^T [theta]_e + [A_3,1:3]_E^T [alpha_3]_e

by disconnecting the drilling degree of freedom from the bending degrees of
freedom in the nodal basis.

The following features are incorporated to deal with nodal normals:
- Nodal normals are averages of the normals of elements that meet at a node.
- A crease in the surface is take into account. In that case the normal are not
  averaged across the crease. At the nodes along the crease every element uses
  the normal to its surface instead of the nodal normal.
"""
mutable struct FEMMShellT3FF{S<:AbstractFESet, F<:Function, M} <: AbstractFEMM
    integdomain::IntegDomain{S, F} # integration domain data
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object.
    transv_shear_formulation::FInt
    drilling_stiffness_scale::Float64
    threshold_angle::Float64
    mult_el_size::Float64
    _associatedgeometry::Bool
    _normals::FFltMat
    _normal_valid::Vector{Bool}
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
    FEMMShellT3FF(integdomain::IntegDomain{S, F}, mcsys::CSys, material::M) where {S<:AbstractFESet, F<:Function, M}

Constructor of the T3FF shell FEMM.
"""
function FEMMShellT3FF(integdomain::IntegDomain{S, F}, mcsys::CSys, material::M) where {S<:AbstractFESet, F<:Function, M}
    _nnmax = 0
    for j in 1:count(integdomain.fes)
        for k in eachindex(integdomain.fes.conn[j])
            _nnmax = max(_nnmax, integdomain.fes.conn[j][k])
        end
    end
    _normals = fill(0.0, _nnmax, 3)
    _normal_valid = fill(true, _nnmax)
    # Alocate the buffers
    _loc = fill(0.0, 1, 3)
    _J0 = fill(0.0, 3, 2)
    _ecoords = fill(0.0, __nn, 3)
    _edisp = fill(0.0, __nn*__ndof); 
    _ecoords_e = fill(0.0, __nn, 2) 
    _edisp_e = fill(0.0, __nn*__ndof); 
    _dofnums = zeros(FInt, 1, __nn*__ndof); 
    _E_G = fill(0.0, 3, 3); 
    _A_Es = [fill(0.0, 3, 3), fill(0.0, 3, 3), fill(0.0, 3, 3)]; 
    _nvalid = fill(false, 3)
    _T = fill(0.0, __nn*__ndof, __nn*__ndof)
    _elmat = fill(0.0, __nn*__ndof, __nn*__ndof);   
    _gradN_e = fill(0.0, __nn, 2)
    _Bm = fill(0.0, 3, __nn*__ndof)
    _Bb = fill(0.0, 3, __nn*__ndof)
    _Bs = fill(0.0, 2, __nn*__ndof)
    _DpsBmb = similar(_Bm)
    _DtBs = similar(_Bs)

    return FEMMShellT3FF(integdomain, mcsys, 
        material, 
        __TRANSV_SHEAR_FORMULATION_AVERAGE_B, 1.0, 30.0, 5/12/1.5,
        false,
        _normals, _normal_valid,
        _loc, _J0,
        _ecoords, _edisp, _ecoords_e, _edisp_e,
        _dofnums, 
        _E_G, _A_Es, _nvalid, _T,
        _elmat, 
        _gradN_e,
        _Bm, _Bb, _Bs, _DpsBmb, _DtBs)
end

function isoparametric!(E_G::FFltMat, XYZ::FFltMat, J0::FFltMat, fe_label::FInt)
    # This is the tangent to the coordinate curve 1
    a = @view J0[:, 1]
    L0 = norm(a);
    E_G[:,1] = a/L0;
    # This is the tangent to the coordinate curve 2
    b = J0[:, 2]
    #     E_G(:,3)=skewmat(E_G(:,1))*b;
    b = @view J0[:, 2]
    # Now compute the normal
    E_G[:, 3] .= (-E_G[3,1]*b[2]+E_G[2,1]*b[3],
                   E_G[3,1]*b[1]-E_G[1,1]*b[3],
                  -E_G[2,1]*b[1]+E_G[1,1]*b[2]);
    E_G[:, 3] /= norm(@view E_G[:, 3]);
    #     E_G(:,2)=skewmat(E_G(:,3))*E_G(:,1);
    E_G[:, 2] .= (-E_G[3,3]*E_G[2,1]+E_G[2,3]*E_G[3,1],
                   E_G[3,3]*E_G[1,1]-E_G[1,3]*E_G[3,1],
                  -E_G[2,3]*E_G[1,1]+E_G[1,3]*E_G[2,1]);
    return E_G
end

function _compute_nodal_normal!(n, mcsys::CSys, XYZ, J0::FFltMat, labl::FInt)
    updatecsmat!(mcsys, reshape(XYZ, 1, 3), J0, labl);
    n[:] .= mcsys.csmat[:, 3]
    return n
end

"""
    make(integdomain, mcsys, material)

Make a T3FF FEMM from the integration domain,  and a material.
Default isoparametric method for computing the normals is used.
"""
function FEMMShellT3FF(integdomain::IntegDomain{S, F}, material::M) where {S<:AbstractFESet, F<:Function, M}
    return FEMMShellT3FF(integdomain, CSys(3, 3, isoparametric!), material)
end


"""
    make(integdomain, mcsys, material)

Make a T3FF FEMM from the integration domain, a coordinate system to define the
orientation of the normals, and a material.
"""
function make(integdomain, mcsys, material)
    return FEMMShellT3FF(integdomain, mcsys, material)
end

function _compute_J!(J0, ecoords)
    # Compute the Jacobian matrix: the Jacobian matrix is constant within the
    # triangle
    x, y, z = ecoords[2, :].-ecoords[1, :]
    J0[:, 1] .= (x, y, z)
    x, y, z = ecoords[3, :].-ecoords[1, :]
    J0[:, 2] .= (x, y, z)
    return J0
end

function _e_g!(E_G, J0)
    # This is the tangent to the coordinate curve 1
    a = @view J0[:, 1]
    L0 = norm(a);
    E_G[:,1] = a/L0;
    # This is the tangent to the coordinate curve 2
    b = @view J0[:, 2]
    # Now compute the normal
    E_G[:, 3] .= (-E_G[3, 1]*b[2]+E_G[2, 1]*b[3],
                   E_G[3, 1]*b[1]-E_G[1, 1]*b[3],
                  -E_G[2, 1]*b[1]+E_G[1, 1]*b[2]);
    E_G[:, 3] /= norm(@view E_G[:, 3]);
    E_G[:, 2] .= (-E_G[3, 3]*E_G[2, 1]+E_G[2, 3]*E_G[3, 1],
                   E_G[3, 3]*E_G[1, 1]-E_G[1, 3]*E_G[3, 1],
                  -E_G[2, 3]*E_G[1, 1]+E_G[1, 3]*E_G[2, 1]);
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
    centroid[:] .= (
        (ecoords[1, 1]+ecoords[2, 1]+ecoords[3, 1])/3,
        (ecoords[1, 2]+ecoords[2, 2]+ecoords[3, 2])/3,
        (ecoords[1, 3]+ecoords[2, 3]+ecoords[3, 3])/3
        )
end
    
function _shell_material_stiffness(material)
    # Compute the two material stiffness matrices: the plane stress matrix, and
    # the transverse shear matrix
    D = fill(0.0, 6, 6)
    t::FFlt, dt::FFlt, loc::FFltMat, label::FInt = 0.0, 0.0, [0.0 0.0 0.0], 0
    tangentmoduli!(material,  D,  t, dt, loc, label)
    Dps = fill(0.0, 3, 3)
    Dps[1:2, 1:2] = D[1:2, 1:2] -  (reshape(D[1:2,3], 2, 1) * reshape(D[3,1:2], 1, 2))/D[3, 3]
    ix=[1, 2, 4];
    for i = 1:3
        Dps[3,i] = Dps[i,3] = D[4, ix[i]];
    end
    Dt = fill(0.0, 2, 2)
    ix=[5, 6];
    for i = 1:2
        Dt[i,i] = D[ix[i], ix[i]];
    end
    return Dps, Dt
end

function _nodal_triads_e!(A_Es, nvalid, E_G, normals, normal_valid, c)
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

    # TO DO Get rid of the temporaries
    r = fill(0.0, 3)
    nk_e = fill(0.0, 3)
    nk = fill(0.0, 3)
    f3_e = [0.0, 0.0, 1.0]
    for k in 1:length(A_Es)
        nk .= vec(view(normals, c[k], :))
        nvalid[k] = normal_valid[c[k]]
        if nvalid[k] 
            # If the nodal normal is valid, pull it back into the element frame.
            nk_e .= E_G'*nk
        else
            # Otherwise the element normal replaces the nodal normal.
            nk_e .= 0.0; nk_e[3] = 1.0
        end
        cross3!(r, f3_e, nk_e)
        if norm(r) > 1.0e-12
            rotmat3!(A_Es[k], r) 
        else
            A_Es[k] .= I(3)
        end
    end
    return A_Es, nvalid
end

function _transfmat_g_to_a!(T, A_Es, E_G)
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
        Tblock = transpose(A_Es[i]) * transpose(E_G) # TO DO remove temporary
        offset = (i-1)*__ndof
        r = offset+1:offset+3
        @. T[r, r] = Tblock
        r = offset+4:offset+6
        @. T[r, r] = Tblock
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
    
    # TO DO avoid a temporary
    T .= 0.0
    for i in 1:__nn
        roffst = (i-1)*__ndof
        iA_E = A_Es[i] 
        iA_E_33 = iA_E[3, 3]
        r = roffst+1:roffst+3
        T[r, r] .= iA_E
        r = roffst+4:roffst+5
        # TO DO avoid temporary
        T[r, r] .= iA_E[1:2, 1:2] - (1/iA_E_33).*(vec(iA_E[1:2, 3])*vec(iA_E[3, 1:2])') 
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
    a, b = ecoords_e[2, :] .- ecoords_e[1, :]
    c, d = ecoords_e[3, :] .- ecoords_e[1, :]
    J = (a*d - b*c)
    gradN_e[:] .= ((b-d)/J, d/J, -b/J, (c-a)/J, -c/J, a/J)
    return gradN_e, J/2
end

function _add_Bsmat_o!(Bs, ecoords_e, Ae, ordering)
    # Compute the linear transverse shear strain-displacement matrix for one
    # particular ordering of the nodes. The computed entries are ADDED. The
    # matrix is NOT zeroed out initially. That is the responsibility of the
    # caller.

    # Orientation 
    s, p, q = ordering
    a, b = ecoords_e[p, :] .- ecoords_e[s, :]
    c, d = ecoords_e[q, :] .- ecoords_e[s, :]
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
    associategeometry!(self::FEMMShellT3FF,  geom::NodalField{FFlt})

Associate geometry with the FEMM. 

In this case it means evaluate the nodal normals.
"""
function associategeometry!(self::FEMMShellT3FF,  geom::NodalField{FFlt})
    threshold_angle = self.threshold_angle
    J0 = self._J0
    E_G = self._E_G
    normals, normal_valid = self._normals, self._normal_valid
    nnormal = fill(0.0, 3)

    # Compute the normals at the nodes
    for el in 1:count(self.integdomain.fes)
        i, j, k = self.integdomain.fes.conn[el]
        J0[:, 1] = geom.values[j, :] - geom.values[i, :]
        J0[:, 2] = geom.values[k, :] - geom.values[i, :]
        for n in [i, j, k]
            _compute_nodal_normal!(nnormal, self.mcsys, geom.values[n, :], J0, self.integdomain.fes.label[el])
            normals[n, :] .+= nnormal
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
    num_normals(self::FEMMShellT3FF)

Compute the summary of the nodal normals.
"""
function num_normals(self::FEMMShellT3FF)
    normals, normal_valid = self._normals, self._normal_valid
    total_normals = length(normal_valid)
    total_invalid_normals = length([v for v in normal_valid if v != true])
    return total_normals, total_invalid_normals
end

"""
    stiffness(self::FEMMShellT3FF, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{TI}, Rfield1::NodalField{TI}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, TI<:Number}

Compute the material stiffness matrix.
"""
function stiffness(self::FEMMShellT3FF, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{TI}, Rfield1::NodalField{TI}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, TI<:Number}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    label = self.integdomain.fes.label
    normals, normal_valid = self._normals, self._normal_valid
    centroid, J0 = self._loc, self._J0
    ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
    ecoords_e, gradN_e = self._ecoords_e, self._gradN_e 
    E_G, A_Es, nvalid, T = self._E_G, self._A_Es, self._nvalid, self._T
    elmat = self._elmat
    transformwith = Transformer(elmat)
    Bm, Bb, Bs, DpsBmb, DtBs = self._Bm, self._Bb, self._Bs, self._DpsBmb, self._DtBs
    Dps, Dt = _shell_material_stiffness(self.material)
    scf = 5/6;  # shear correction factor
    Dt .*= scf
    drilling_stiffness_scale = self.drilling_stiffness_scale
    transv_shear_formulation = self.transv_shear_formulation
    mult_el_size = self.mult_el_size
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), dchi.nfreedofs, dchi.nfreedofs);
    for i in 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords, fes.conn[i]);
        _centroid!(centroid, ecoords)
        _compute_J!(J0, ecoords)
        _e_g!(E_G, J0)
        _ecoords_e!(ecoords_e, J0, E_G)
        gradN_e, Ae = _gradN_e_Ae!(gradN_e, ecoords_e)
        t = self.integdomain.otherdimension(centroid, fes.conn[i], [1.0/3 1.0/3])
        # Construct the Stiffness Matrix
        fill!(elmat,  0.0); # Initialize element matrix
        _Bmmat!(Bm, gradN_e)
        add_btdb_ut_only!(elmat, Bm, t*Ae, Dps, DpsBmb)
        _Bbmat!(Bb, gradN_e)
        add_btdb_ut_only!(elmat, Bb, (t^3)/12*Ae, Dps, DpsBmb)
        he = sqrt(2*Ae)
        if transv_shear_formulation == __TRANSV_SHEAR_FORMULATION_AVERAGE_K
            Bs .= 0.0; _add_Bsmat_o!(Bs, ecoords_e, Ae, (1, 2, 3))
            add_btdb_ut_only!(elmat, Bs, (t^3/(t^2+mult_el_size*he^2))*Ae/3, Dt, DtBs)
            Bs .= 0.0; _add_Bsmat_o!(Bs, ecoords_e, Ae, (2, 3, 1))
            add_btdb_ut_only!(elmat, Bs, (t^3/(t^2+mult_el_size*he^2))*Ae/3, Dt, DtBs)
            Bs .= 0.0; _add_Bsmat_o!(Bs, ecoords_e, Ae, (3, 1, 2))
            add_btdb_ut_only!(elmat, Bs, (t^3/(t^2+mult_el_size*he^2))*Ae/3, Dt, DtBs)
        else
            _Bsmat!(Bs, ecoords_e, Ae)
            add_btdb_ut_only!(elmat, Bs, (t^3/(t^2+mult_el_size*he^2))*Ae, Dt, DtBs)
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
        # Add the artificial drilling stiffness in the nodal basis, but only if
        # the nodal normal is valid. 
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
    return makematrix!(assembler);
end

function stiffness(self::FEMMShellT3FF, geom0::NodalField{FFlt}, u1::NodalField{TI}, Rfield1::NodalField{TI}, dchi::NodalField{TI}) where {TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi);
end


"""
    mass(self::FEMMShellT3FF,  assembler::A,  geom::NodalField{FFlt}, dchi::NodalField{TI}) where {A<:AbstractSysmatAssembler, TI<:Number}

Compute the diagonal (lumped) mass matrix

The mass matrix can be expected to be non-singular.
"""
function mass(self::FEMMShellT3FF,  assembler::A,  geom0::NodalField{FFlt}, dchi::NodalField{TI}) where {A<:AbstractSysmatAssembler, TI<:Number}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    label = self.integdomain.fes.label
    normals, normal_valid = self._normals, self._normal_valid
    centroid, J0 = self._loc, self._J0
    ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
    ecoords_e, gradN_e = self._ecoords_e, self._gradN_e 
    E_G, A_Es, nvalid, T = self._E_G, self._A_Es, self._nvalid, self._T
    elmat = self._elmat
    transformwith = Transformer(elmat)
    rho::FFlt = massdensity(self.material); # mass density
    npe = nodesperelem(fes)
    ndn = ndofs(dchi)
    drilling_stiffness_scale = self.drilling_stiffness_scale
    startassembly!(assembler,  size(elmat,1),  size(elmat,2),  count(fes), dchi.nfreedofs,  dchi.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords, fes.conn[i]);
        _centroid!(centroid, ecoords)
        _compute_J!(J0, ecoords)
        _e_g!(E_G, J0)
        _ecoords_e!(ecoords_e, J0, E_G)
        gradN_e, Ae = _gradN_e_Ae!(gradN_e, ecoords_e)
        # Compute the translational and rotational masses corresponding to nodes
        t = self.integdomain.otherdimension(centroid, fes.conn[i], [1.0/3 1.0/3])
        tfactor = rho*(t*Ae)/3;
        rfactor = rho*(t^3/12*Ae)/3;
        # end # Loop over quadrature points
        # Now treat the transformation from the element to the nodal triad
        _nodal_triads_e!(A_Es, nvalid, E_G, normals, normal_valid, fes.conn[i])
        # Fill the elementwise matrix in the nodal basis: the matrix is lumped,
        # hence diagonal, and no transformation from element to nodal is needed.
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
    return makematrix!(assembler);
end

function mass(self::FEMMShellT3FF,  geom::NodalField{FFlt},  u::NodalField{TI}) where {TI<:Number}
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
function inspectintegpoints(self::FEMMShellT3FF, geom0::NodalField{FFlt},  u::NodalField{TI}, dT::NodalField{FFlt}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {TI<:Number, F<:Function}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    label = self.integdomain.fes.label
    normals, normal_valid = self._normals, self._normal_valid
    centroid, J0 = self._loc, self._J0
    ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
    ecoords_e, gradN_e = self._ecoords_e, self._gradN_e 
    E_G, A_Es, nvalid, T = self._E_G, self._A_Es, self._nvalid, self._T
    elmat = self._elmat
    transformwith = Transformer(elmat)
    Bm, Bb, Bs, DpsBmb, DtBs = self._Bm, self._Bb, self._Bs, self._DpsBmb, self._DtBs
    Dps, Dt = _shell_material_stiffness(self.material)
    scf = 5/6;  # shear correction factor
    Dt .*= scf
    mult_el_size = self.mult_el_size
    edisp_e = deepcopy(edisp)
    edisp_n = deepcopy(edisp_e)
    out = fill(0.0, 3)
    # Sort out  the output requirements
    outputcsys = self.mcsys; # default: report the stresses in the material coord system
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
        o_e = E_G' * outputcsys.csmat
        o2_e = o_e[1:2, 1:2]
        # Compute the Requested Quantity
        if quant == BENDING_MOMENT
            _Bbmat!(Bb, gradN_e)
            kurv = Bb * edisp_e
            mom = ((t^3)/12)*Dps * kurv
            m = [mom[1] mom[3]; mom[3] mom[2]]
            mo = o2_e' * m * o2_e
            out[:] .= mo[1, 1], mo[2, 2], mo[1, 2]
        end
        if quant == MEMBRANE_FORCE
            _Bmmat!(Bm, gradN_e)
            strn = Bm * edisp_e
            frc = (t)*Dps * strn
            f = [frc[1] frc[3]; frc[3] frc[2]]
            fo = o2_e' * f * o2_e
            # @infiltrate
            out[:] .= fo[1, 1], fo[2, 2], fo[1, 2]
        end
        if quant == TRANSVERSE_SHEAR
            _Bsmat!(Bs, ecoords_e, Ae)
            he = sqrt(2*Ae)
            shr = Bs * edisp_e
            frc = ((t^3/(t^2+mult_el_size*he^2)))*Dt * shr
            fo = o2_e' * frc
            out[1:2] .= fo[1], fo[2]
        end 
        # Call the inspector
        idat = inspector(idat, i, fes.conn[i], ecoords, out, centroid);
        # end # Loop over quadrature points
    end # Loop over elements
    return idat; # return the updated inspector data
end

function inspectintegpoints(self::FEMMShellT3FF, geom::NodalField{FFlt},  u::NodalField{TI}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {TI<:Number, F<:Function}
    dT = NodalField(fill(zero(FFlt), nnodes(geom), 1)) # zero difference in temperature
    return inspectintegpoints(self, geom, u, dT, felist, inspector, idat, quantity; context...);
end

# function _resultant_check(self::FEMMShellT3FF, geom0::NodalField{FFlt}, u1::NodalField{TI}, Rfield1::NodalField{TI}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, TI<:Number}
#     @assert self._associatedgeometry == true
#     fes = self.integdomain.fes
#     label = self.integdomain.fes.label
#     normals, normal_valid = self._normals, self._normal_valid
#     centroid, J0 = self._loc, self._J0
#     ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
#     ecoords_e, gradN_e = self._ecoords_e, self._gradN_e 
#     E_G, A_Es, nvalid, T = self._E_G, self._A_Es, self._nvalid, self._T
#     elmat = self._elmat
#     transformwith = Transformer(elmat)
#     Bm, Bb, Bs, DpsBmb, DtBs = self._Bm, self._Bb, self._Bs, self._DpsBmb, self._DtBs
#     Dps, Dt = _shell_material_stiffness(self.material)
#     scf = 5/6;  # shear correction factor
#     Dt .*= scf
#     drilling_stiffness_scale = self.drilling_stiffness_scale
#     mult_el_size = self.mult_el_size
#     edisp = fill(0.0, 18)
#     nodal_moments = 0.0 .* deepcopy(normals)
#     nodal_rotations = 0.0 .* deepcopy(normals)
#     for i in 1:count(fes) # Loop over elements
#         gathervalues_asmat!(geom0, ecoords, fes.conn[i]);
#         gathervalues_asvec!(dchi, edisp, fes.conn[i]);
#         _centroid!(centroid, ecoords)
#         _compute_J!(J0, ecoords)
#         _e_g!(E_G, J0)
#         _ecoords_e!(ecoords_e, J0, E_G)
#         gradN_e, Ae = _gradN_e_Ae!(gradN_e, ecoords_e)
#         t = self.integdomain.otherdimension(centroid, fes.conn[i], [1.0/3 1.0/3])
#         # Construct the Stiffness Matrix
#         fill!(elmat,  0.0); # Initialize element matrix
#         _Bmmat!(Bm, gradN_e)
#         add_btdb_ut_only!(elmat, Bm, t*Ae, Dps, DpsBmb)
#         _Bbmat!(Bb, gradN_e)
#         add_btdb_ut_only!(elmat, Bb, (t^3)/12*Ae, Dps, DpsBmb)
#         he = sqrt(2*Ae)
#         _Bsmat!(Bs, ecoords_e, Ae)
#         add_btdb_ut_only!(elmat, Bs, (t^3/(t^2+mult_el_size*he^2))*Ae, Dt, DtBs)
#         # Complete the elementwise matrix by filling in the lower triangle
#         complete_lt!(elmat)
#         # Now treat the transformation from the element to the nodal triad
#         _nodal_triads_e!(A_Es, nvalid, E_G, normals, normal_valid, fes.conn[i])
#         _transfmat_a_to_e!(T, A_Es, gradN_e)
#         # Transform the elementwise matrix into the nodal coordinates
#         transformwith(elmat, T)
#         # Bending diagonal stiffness coefficients
#         kavg = mean((elmat[4, 4], elmat[10, 10], elmat[16, 16],
#             elmat[5, 5], elmat[11, 11], elmat[17, 17])) * drilling_stiffness_scale
#         # Add the artificial drilling stiffness in the nodal cordinates, but
#         # only if the nodal normal is valid. 
#         nvalid[1] && (elmat[6,6] += kavg)
#         nvalid[2] && (elmat[12,12] += kavg)
#         nvalid[3] && (elmat[18,18] += kavg)
#         # Transform from nodal into global coordinates
#         _transfmat_g_to_a!(T, A_Es, E_G)
#         transformwith(elmat, T)
#         # Compute the resultants
#         eforc = elmat * edisp
#         nodal_moments[fes.conn[i][1], :] += eforc[4:6]
#         nodal_moments[fes.conn[i][2], :] += eforc[10:12]
#         nodal_moments[fes.conn[i][3], :] += eforc[16:18]
#         nodal_rotations[fes.conn[i][1], :] += edisp[4:6]
#         nodal_rotations[fes.conn[i][2], :] += edisp[10:12]
#         nodal_rotations[fes.conn[i][3], :] += edisp[16:18]
#     end # Loop over elements
#     nodal_moment_magnitude = 0.0
#     nodal_moment_normal_magnitude = 0.0
#     for k in 1:size(nodal_moments, 1)
#         nodal_moment_magnitude += norm(nodal_moments[k, :])
#         nodal_moment_normal_magnitude += abs(dot(vec(normals[k, :]), vec(nodal_moments[k, :])))
#     end
#     nodal_moment_magnitude /= size(nodal_moments, 1)
#     nodal_moment_normal_magnitude /= size(nodal_moments, 1)
#     @show nodal_moment_magnitude,     nodal_moment_normal_magnitude
#     nodal_rotation_magnitude = 0.0
#     nodal_rotation_normal_magnitude = 0.0
#     for k in 1:size(nodal_rotations, 1)
#         nodal_rotation_magnitude += norm(nodal_rotations[k, :])
#         nodal_rotation_normal_magnitude += abs(dot(vec(normals[k, :]), vec(nodal_rotations[k, :])))
#     end
#     nodal_rotation_magnitude /= size(nodal_rotations, 1)
#     nodal_rotation_normal_magnitude /= size(nodal_rotations, 1)
#     @show nodal_rotation_magnitude,     nodal_rotation_normal_magnitude
#     return true
# end

end # module

# @show nvalid
# @show [v for v in abs.(eigen(elmat).values) if v > 1.0e-3]


"""
    _Bsmat!(Bs, gradN, N)

Compute the linear transverse shear strain-displacement matrix.

International Journal of Applied Mechanics
Vol. 9, No. 4 (2017) 1750055 (30 pages)
c⃝ World Scientific Publishing Europe Ltd.
DOI: 10.1142/S1758825117500557
A Central Point-Based Discrete Shear Gap Method
for Plates and Shells Analysis Using
Triangular Elements
X. Y. Cui∗ and L. Tian
State Key Laboratory of Advanced Design and
Manufacturing for Vehicle Body Hunan University
Changsha 410082, P. R. China
∗cuixy@hnu.edu.cn
Received 18 February 2017

This matrix computation is precisely the same as the original averaging
I implemented, even though this is based on the central point.

"""
# function _Bsmat!(Bs, ecoords_e)
#     x = @view ecoords_e[:, 1]
#     y = @view ecoords_e[:, 2]
#     Vx, Vy = x[2] .- x[1], y[2] .- y[1]
#     Wx, Wy = x[3] .- x[1], y[3] .- y[1]
#     J = (Vx*Wy - Vy*Wx)
#     m = ((y[2]-y[3])/J, (y[3]-y[1])/J, (y[1]-y[2])/J)
#     n = (-(x[2]-x[3])/J, -(x[3]-x[1])/J, -(x[1]-x[2])/J)
#     xo = (x[1]+x[2]+x[3])/3
#     yo = (y[1]+y[2]+y[3])/3
#     a = x[1]-xo
#     b = y[1]-yo
#     c = x[2]-xo
#     d = y[2]-yo
#     e = x[3]-xo
#     f = y[3]-yo
#     fill!(Bs, 0.0)
#     s = 6*(1-1)
#     Bs[1, s+3] = 2/3*m[1]-1/3*m[2]-1/3*m[3]
#     Bs[1, s+4] = -2/3*b*m[1]-1/6*d*m[2]-1/6*f*m[3]
#     Bs[1, s+5] = +2/3*a*m[1]+1/6*c*m[2]+1/6*e*m[3]
#     Bs[2, s+3] = 2/3*n[1]-1/3*n[2]-1/3*n[3]
#     Bs[2, s+4] = -2/3*b*n[1]-1/6*d*n[2]-1/6*f*n[3]
#     Bs[2, s+5] = +2/3*a*n[1]+1/6*c*n[2]+1/6*e*n[3]
#     s = 6*(2-1)
#     Bs[1, s+3] = -1/3*m[1]+2/3*m[2]-1/3*m[3]
#     Bs[1, s+4] = -1/6*b*m[1]-2/3*d*m[2]-1/6*f*m[3]
#     Bs[1, s+5] = +1/6*a*m[1]+2/3*c*m[2]+1/6*e*m[3]
#     Bs[2, s+3] = -1/3*n[1]+2/3*n[2]-1/3*n[3]
#     Bs[2, s+4] = -1/6*b*n[1]-2/3*d*n[2]-1/6*f*n[3]
#     Bs[2, s+5] = +1/6*a*n[1]+2/3*c*n[2]+1/6*e*n[3]
#     s = 6*(3-1)
#     Bs[1, s+3] = -1/3*m[1]-1/3*m[2]+2/3*m[3]
#     Bs[1, s+4] = -1/6*b*m[1]-1/6*d*m[2]-2/3*f*m[3]
#     Bs[1, s+5] = +1/6*a*m[1]+1/6*c*m[2]+2/3*e*m[3]
#     Bs[2, s+3] = -1/3*n[1]-1/3*n[2]+2/3*n[3]
#     Bs[2, s+4] = -1/6*b*n[1]-1/6*d*n[2]-2/3*f*n[3]
#     Bs[2, s+5] = +1/6*a*n[1]+1/6*c*n[2]+2/3*e*n[3]
#     return Bs
# end