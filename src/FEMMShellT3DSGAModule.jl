module FEMMShellT3DSGAModule

using LinearAlgebra: norm, Transpose, mul!, diag, eigen, I, dot
using Statistics: mean
using FinEtools
import FinEtools.FESetModule: gradN!, nodesperelem, manifdim
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FEMMBaseModule: associategeometry!
import FinEtools.FEMMBaseModule: inspectintegpoints
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!, update!, thermalstrain!
using FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, locjac!, add_nnt_ut_only!, add_btsigma!
using ..FESetShellT3Module: FESetShellT3, local_frame!
using ..TransformerModule: Transformer

using Infiltrator

const __nn = 3 # number of nodes
const __ndof = 6 # number of degrees of freedom per node

"""
    FEMMShellT3DSGA{S<:AbstractFESet, F<:Function} <: AbstractFEMM

Class for Discrete Shear Gap shell finite element modeling machine. With
averaging of the transverse strain-displacement matrix to provide isotropic
transverse shear response. Also, the formulation is developed to correctly
handle the coupling of twisting moments and transverse shear (such as in the
twisted beam or the Raasch hook problems) by incorporating "nodal" normals.

Programming developed consistently with the paper

[1] Cui et al, Analysis of plates and shells using an edge-based smoothed finite
element method, Comput Mech (2010) 45:141–156 DOI 10.1007/s00466-009-0429-9

In the below reference, the sign next to Ae in equation (44) is wrong:

[2] A superconvergent alpha finite element method (S a FEM) for static and
free vibration analysis of shell structures
Chai et al. (2017).

The stabilization factor of the shear term of

Mikko Lyly, Rolf Stenberg and Teemu Vihinen
A stable bilinear element for the
Reissner-Mindlin plate model
Computer Methods in Applied Mechanics and Engineering 110 (1993) 343-357 

is incorporated. Refer to expressions (3.12) and (3.13).

The treatment of the transformation between the element and nodal coordinates
is carried out using a clean alternative to the publication

[3] Finite Elements in Analysis and Design 30 (1998) 235—242
The treatment of shell normals in ﬁnite element analysis
Richard H. MacNeal, Charles T. Wilson, Robert L. Harder, Claus C. Hoﬀ
The MacNeal-Schwendler Corporation, 815 Colorado Blvd., Los Angeles, CA 90041, USA

The formula for the element to nodal basis transformation is
    [theta]_n = [n]_e^T [theta]_e + [n_3,1:3]_e^T [alpha_3]_e

The following features are incorporated to deal with nodal normals:
- Nodal normals are averages of the normals of elements that meet at a node.
- A crease in the surface is take into account. In that case the normal are not
  averaged across the crease. At the nodes along the crease every element uses
  the normal to its surface instead of the nodal normal.
"""
mutable struct FEMMShellT3DSGA{S<:AbstractFESet, F<:Function, M} <: AbstractFEMM
    integdomain::IntegDomain{S, F} # integration domain data
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object.
    drilling_stiffness_scale::Float64
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
    _e_g::FFltMat
    _n_e::Vector{FFltMat} # transformation nodal-element matrices  
    _nvalid::Vector{Bool}
    _Te::FFltMat
    _elmat::FFltMat
    _gradN_e::FFltMat
    _Bm::FFltMat
    _Bb::FFltMat
    _Bs::FFltMat
    _DpsBmb::FFltMat
    _DtBs::FFltMat
end

function FEMMShellT3DSGA(integdomain::IntegDomain{S, F}, material::M) where {S<:AbstractFESet, F<:Function, M}
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
    _e_g = fill(0.0, 3, 3); 
    _n_e = [fill(0.0, 3, 3), fill(0.0, 3, 3), fill(0.0, 3, 3)]; 
    _nvalid = fill(false, 3)
    _Te = fill(0.0, __nn*__ndof, __nn*__ndof)
    _elmat = fill(0.0, __nn*__ndof, __nn*__ndof);   
    _gradN_e = fill(0.0, __nn, 2)
    _Bm = fill(0.0, 3, __nn*__ndof)
    _Bb = fill(0.0, 3, __nn*__ndof)
    _Bs = fill(0.0, 2, __nn*__ndof)
    _DpsBmb = similar(_Bm)
    _DtBs = similar(_Bs)
    
    return FEMMShellT3DSGA(integdomain, CSys(3), 
        material, 0.1,
        false,
        _normals, _normal_valid,
        _loc, _J0,
        _ecoords, _edisp, _ecoords_e, _edisp_e,
        _dofnums, 
        _e_g, _n_e, _nvalid, _Te,
        _elmat, 
        _gradN_e,
        _Bm, _Bb, _Bs, _DpsBmb, _DtBs)
end

function make(integdomain, material)
    return FEMMShellT3DSGA(integdomain, material)
end

function _compute_J_loc!(J0, loc, ecoords)
    # Compute the Jacobian matrix: the Jacobian matrix is constant within the
    # triangle
    x, y, z = ecoords[2, :].-ecoords[1, :]
    J0[:, 1] .= (x, y, z)
    x, y, z = ecoords[3, :].-ecoords[1, :]
    J0[:, 2] .= (x, y, z)
    loc[1] = (ecoords[1, 1]+ecoords[2, 1]+ecoords[3, 1])/3
    loc[2] = (ecoords[1, 2]+ecoords[2, 2]+ecoords[3, 2])/3
    loc[3] = (ecoords[1, 3]+ecoords[2, 3]+ecoords[3, 3])/3
    return J0, loc
end

function _local_coordinates!(ecoords_e, J0, e_g)
    ecoords_e[1, 1] = ecoords_e[1, 2] = 0.0
    ecoords_e[2, 1] = dot(view(J0, :, 1), view(e_g, :, 1))
    ecoords_e[2, 2] = dot(view(J0, :, 1), view(e_g, :, 2))
    ecoords_e[3, 1] = dot(view(J0, :, 2), view(e_g, :, 1))
    ecoords_e[3, 2] = dot(view(J0, :, 2), view(e_g, :, 2))
    return ecoords_e
end

function _centroid(ecoords)
    (
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

function _nodal_triads_e!(n_e, nvalid, e_g, normals, normal_valid, c)
    # Components of nodal cartesian ordinate systems such that the third
    # direction is the direction of the nodal normal, and the angle to rotate
    # the element normal into the nodal normal is as short as possible; these
    # components are given on the local element coordinate system basis
    # vectors. If the angle of rotation is excessively small, the nodal matrix
    # is the identity. 
    # 
    # These triads are the nodal-to-element transformation
    # matrices.

    # TO DO Get rid of the temporaries
    r = fill(0.0, 3)
    nk_e = fill(0.0, 3)
    nk = fill(0.0, 3)
    f3_e = [0.0, 0.0, 1.0]
    for k in 1:length(n_e)
        nk .= vec(view(normals, c[k], :))
        nvalid[k] = normal_valid[c[k]]
        if nvalid[k] 
            # If the nodal normal is valid, pull it back into the element frame.
            nk_e .= e_g'*nk
        else
            # Otherwise the element normal replaces the nodal normal.
            nk_e .= 0.0; nk_e[3] = 1.0
        end
        cross3!(r, f3_e, nk_e)
        if norm(r) > 1.0e-6
            rotmat3!(n_e[k], r) 
        else
            n_e[k] .= I(3)
        end
    end
    return n_e, nvalid
end

function _transfmat_n_to_g!(Te, n_e, e_g)
    # Nodal-to-global transformation matrix. 
    # The 3x3 blocks consist of the nodal
    # triad expressed on the global basis. The nodal basis vectors in `n_e
    # [i]` are expressed on the element basis, and are then rotated with `e_g`
    # into the global coordinate system.
    # Output
    # - `Te` = transformation matrix, input in the nodal basis, output in the
    #   global basis
    Te .= 0.0
    for i in 1:__nn
        Teblock = e_g * n_e[i] # TO DO remove temporary
        offset = (i-1)*__ndof
        r = offset+1:offset+3
        @. Te[r, r] = Teblock
        r = offset+4:offset+6
        @. Te[r, r] = Teblock
    end
    return Te
end

function _transfmat_e_to_n!(Te, n_e, gradN_e)
    # Element-to-nodal transformation matrix. 
    # - `n_e` = matrix with the nodal triad vectors in columns, components on
    #   the element basis. Its transpose is the element triad on the nodal
    #   basis vectors.
    # - `gradN_e` = basis function gradients on the element basis.
    # Output
    # - `Te` = transformation matrix, input in the element basis, output in the
    #   nodal basis

    # TO DO avoid a temporary
    Te .= 0.0
    # Translation and rotation degrees of freedom
    for i in 1:__nn
        roffset = (i-1)*__ndof
        n_eT = n_e[i]'
        r = roffset+1:roffset+3
        Te[r, r] .= n_eT 
        r = roffset+4:roffset+6
        Te[r, r] .= n_eT 
    end
    # The drilling rotation of the mid surface produced by the 1/2*(v,x - u,y)
    # effect is linked to the out of plane rotations.
    for i in 1:__nn
        coffst = (i-1)*__ndof
        n_eT = n_e[i]'
        a1 = n_e[i][1, 3]
        a2 = n_e[i][2, 3]
        a3 = n_e[i][3, 3]
        for j in 1:__nn
            roffst = (j-1)*__ndof
            Te[roffst+1, coffst+4] = (-a1 * 1/2 * gradN_e[j, 2])
            Te[roffst+2, coffst+4] = (+a1 * 1/2 * gradN_e[j, 1])
            Te[roffst+1, coffst+5] = (-a2 * 1/2 * gradN_e[j, 2])
            Te[roffst+2, coffst+5] = (+a2 * 1/2 * gradN_e[j, 1])
            Te[roffst+1, coffst+6] = (-a3 * 1/2 * gradN_e[j, 2])
            Te[roffst+2, coffst+6] = (+a3 * 1/2 * gradN_e[j, 1])
        end
    end
    return Te
end

function _gradN_e_Ae!(gradN_e, ecoords_e)
    # Compute the gradient with respect to the element ordinates, three
    # gradients in two coordinates, and the area of the triangle
    a, b = ecoords_e[2, :] .- ecoords_e[1, :]
    c, d = ecoords_e[3, :] .- ecoords_e[1, :]
    J = (a*d - b*c)
    gradN_e[:] .= ((b-d)/J, d/J, -b/J, (c-a)/J, -c/J, a/J)
    return gradN_e, J/2
end

"""
    _Bsmat!(Bs, gradN, N)

Compute the linear transverse shear strain-displacement matrix.
"""
function _Bsmat!(Bs, ecoords_e)
    Bs .= 0.0 

    # Orientation 1
    s, p, q = 3, 1, 2
    a, b = ecoords_e[p, :] .- ecoords_e[s, :]
    c, d = ecoords_e[q, :] .- ecoords_e[s, :]
    Ae = (a*d - b*c)/2
    m = (1/2/Ae) * (1/3) # multiplier

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
    
    # Orientation 2
    s, p, q = 1, 2, 3
    a, b = ecoords_e[p, :] .- ecoords_e[s, :]
    c, d = ecoords_e[q, :] .- ecoords_e[s, :]

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


    # Orientation 3
    s, p, q = 2, 3, 1
    a, b = ecoords_e[p, :] .- ecoords_e[s, :]
    c, d = ecoords_e[q, :] .- ecoords_e[s, :]
    
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


"""
    _Bmmat!(Bm, gradN)

Compute the linear membrane strain-displacement matrix.
"""
function _Bmmat!(Bm, gradN)
    fill!(Bm, 0.0)
    for i in 1:__nn
        Bm[1,6*(i-1)+1] = gradN[i,1];
                                             Bm[2,6*(i-1)+2] = gradN[i,2];
        Bm[3,6*(i-1)+1] = gradN[i,2];        Bm[3,6*(i-1)+2] = gradN[i,1];
    end
end

"""
    _Bbmat!(Bb, gradN)

Compute the linear, displacement independent, curvature-displacement/rotation matrix for a shell quadrilateral element with nfens=3 nodes. Displacements and rotations are in a local coordinate system.
"""
function _Bbmat!(Bb, gradN)
    fill!(Bb, 0.0)
    for i in 1:__nn
                                              Bb[1,6*(i-1)+5] = gradN[i,1];
        Bb[2,6*(i-1)+4] = -gradN[i,2];
        Bb[3,6*(i-1)+4] = -gradN[i,1];        Bb[3,6*(i-1)+5] = gradN[i,2];
    end
end

function associategeometry!(self::FEMMShellT3DSGA,  geom::NodalField{FFlt})
    J0 = self._J0
    e_g = self._e_g
    normals, normal_valid = self._normals, self._normal_valid
    # Compute the normals at the nodes
    for el in 1:count(self.integdomain.fes)
        i, j, k = self.integdomain.fes.conn[el]
        J0[:, 1] = geom.values[j, :] - geom.values[i, :]
        J0[:, 2] = geom.values[k, :] - geom.values[i, :]
        local_frame!(delegateof(self.integdomain.fes), e_g, J0)
        normals[i, :] .+= e_g[:, 3]
        normals[j, :] .+= e_g[:, 3]
        normals[k, :] .+= e_g[:, 3]
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
    # that the nodal normal should not be used at that vertex. TO DO: Should we
    # make this a configurable parameter?
    ntolerance = 1 - sqrt(1 - sin(30/180*pi)^2)
    for el in 1:count(self.integdomain.fes)
        i, j, k = self.integdomain.fes.conn[el]
        J0[:, 1] = geom.values[j, :] - geom.values[i, :]
        J0[:, 2] = geom.values[k, :] - geom.values[i, :]
        local_frame!(delegateof(self.integdomain.fes), e_g, J0)
        en = view(e_g, :, 3)
        nd = dot(normals[i, :], en)
        if nd < 1 - ntolerance
            normal_valid[i] = false
        end
        nd = dot(normals[j, :], en)
        if nd < 1 - ntolerance
            normal_valid[j] = false
        end
        nd = dot(normals[k, :], en)
        if nd < 1 - ntolerance
            normal_valid[k] = false
        end
    end
    nz = 0
    for j in 1:size(normals, 1)
        if !normal_valid[j]
            nz += 1
        end
    end
    # @show size(normals, 1) , nz
    self._associatedgeometry = true
    return self
end

"""
    stiffness(self::FEMMShellT3DSGA, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the material stiffness matrix.
"""
function stiffness(self::FEMMShellT3DSGA, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    normals, normal_valid = self._normals, self._normal_valid
    loc, J0 = self._loc, self._J0
    ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
    ecoords_e, gradN_e = self._ecoords_e, self._gradN_e 
    e_g, n_e, nvalid, Te = self._e_g, self._n_e, self._nvalid, self._Te
    elmat = self._elmat
    transformwith = Transformer(elmat)
    Bm, Bb, Bs, DpsBmb, DtBs = self._Bm, self._Bb, self._Bs, self._DpsBmb, self._DtBs
    Dps, Dt = _shell_material_stiffness(self.material)
    scf=5/6;  # shear correction factor
    Dt .*= scf
    drilling_stiffness_scale = self.drilling_stiffness_scale
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), dchi.nfreedofs, dchi.nfreedofs);
    for i in 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords, fes.conn[i]);
        _compute_J_loc!(J0, loc, ecoords)
        local_frame!(delegateof(fes), e_g, J0)
        _local_coordinates!(ecoords_e, J0, e_g)
        gradN_e, Ae = _gradN_e_Ae!(gradN_e, ecoords_e)
        t = self.integdomain.otherdimension(loc, fes.conn[i], [1.0/3 1.0/3])
        # Construct the Stiffness Matrix
        fill!(elmat,  0.0); # Initialize element matrix
        _Bmmat!(Bm, gradN_e)
        _Bbmat!(Bb, gradN_e)
        _Bsmat!(Bs, ecoords_e)
        add_btdb_ut_only!(elmat, Bm, t*Ae, Dps, DpsBmb)
        add_btdb_ut_only!(elmat, Bb, (t^3)/12*Ae, Dps, DpsBmb)
        he = sqrt(2*Ae)
        add_btdb_ut_only!(elmat, Bs, (t^3/(t^2+0.2*he^2))*Ae, Dt, DtBs)
        # Complete the elementwise matrix by filling in the lower triangle
        complete_lt!(elmat)
        # Now treat the transformation from the element to the nodal triad
        _nodal_triads_e!(n_e, nvalid, e_g, normals, normal_valid, fes.conn[i])
        _transfmat_e_to_n!(Te, n_e, gradN_e)
        # Transform the elementwise matrix into the nodal coordinates
        transformwith(elmat, Te)
        # Bending diagonal stiffness coefficients
        kavg = mean((elmat[4, 4], elmat[10, 10], elmat[16, 16],
            elmat[5, 5], elmat[11, 11], elmat[17, 17])) * drilling_stiffness_scale
        # Add the artificial drilling stiffness in the nodal cordinates, but
        # only if the nodal normal is valid. 
        nvalid[1] && (elmat[6,6] += kavg)
        nvalid[2] && (elmat[12,12] += kavg)
        nvalid[3] && (elmat[18,18] += kavg)
        # Transform from nodal into global coordinates
        _transfmat_n_to_g!(Te, n_e, e_g)
        transformwith(elmat, Te)
        # Assembly
        gatherdofnums!(dchi, dofnums, fes.conn[i]); 
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function stiffness(self::FEMMShellT3DSGA, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi);
end


"""
    mass(self::FEMMShellT3DSGA,  assembler::A,  geom::NodalField{FFlt}, dchi::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the consistent mass matrix

This is a general routine for the shell FEMM.
"""
function mass(self::FEMMShellT3DSGA,  assembler::A,  geom0::NodalField{FFlt}, dchi::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    normals, normal_valid = self._normals, self._normal_valid
    loc, J0 = self._loc, self._J0
    ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
    ecoords_e, gradN_e = self._ecoords_e, self._gradN_e 
    e_g, n_e, nvalid, Te = self._e_g, self._n_e, self._nvalid, self._Te
    elmat = self._elmat
    transformwith = Transformer(elmat)
    rho::FFlt = massdensity(self.material); # mass density
    npe = nodesperelem(fes)
    ndn = ndofs(dchi)
    startassembly!(assembler,  size(elmat,1),  size(elmat,2),  count(fes), dchi.nfreedofs,  dchi.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords, fes.conn[i]);
        _compute_J_loc!(J0, loc, ecoords)
        local_frame!(delegateof(fes), e_g, J0)
        _local_coordinates!(ecoords_e, J0, e_g)
        gradN_e, Ae = _gradN_e_Ae!(gradN_e, ecoords_e)
        # Compute the translational and rotational masses corresponding to nodes
        t = self.integdomain.otherdimension(loc, fes.conn[i], [1.0/3 1.0/3])
        tfactor = rho*(t*Ae)/3;
        rfactor = rho*(t^3/12*Ae)/3;
        # end # Loop over quadrature points
        # Now treat the transformation from the element to the nodal triad
        _nodal_triads_e!(n_e, nvalid, e_g, normals, normal_valid, fes.conn[i])
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
                elmat[c, c] += rfactor / 1e6
            end
        end
        # Transform into global coordinates
        _transfmat_n_to_g!(Te, n_e, e_g)
        transformwith(elmat, Te)
        # Assemble
        gatherdofnums!(dchi,  dofnums,  fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler,  elmat,  dofnums,  dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function mass(self::FEMMShellT3DSGA,  geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return mass(self, assembler, geom, u);
end


"""
    inspectintegpoints(self::AbstractFEMMDeforLinear,
      geom::NodalField{FFlt},  u::NodalField{T},
      dT::NodalField{FFlt},
      felist::FIntVec,
      inspector::F,  idat, quantity=:Cauchy;
      context...) where {T<:Number, F<:Function}

Inspect integration point quantities.

- `geom` - reference geometry field
- `u` - displacement field
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
function inspectintegpoints(self::FEMMShellT3DSGA, geom0::NodalField{FFlt},  u::NodalField{T}, dT::NodalField{FFlt}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {T<:Number, F<:Function}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    normals, normal_valid = self._normals, self._normal_valid
    loc, J0 = self._loc, self._J0
    ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
    ecoords_e, gradN_e = self._ecoords_e, self._gradN_e 
    e_g, n_e, nvalid, Te = self._e_g, self._n_e, self._nvalid, self._Te
    elmat = self._elmat
    transformwith = Transformer(elmat)
    Bm, Bb, Bs, DpsBmb, DtBs = self._Bm, self._Bb, self._Bs, self._DpsBmb, self._DtBs
    Dps, Dt = _shell_material_stiffness(self.material)
    scf=5/6;  # shear correction factor
    Dt .*= scf
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
        _compute_J_loc!(J0, loc, ecoords)
        local_frame!(delegateof(fes), e_g, J0)
        _local_coordinates!(ecoords_e, J0, e_g)
        gradN_e, Ae = _gradN_e_Ae!(gradN_e, ecoords_e)
        t = self.integdomain.otherdimension(loc, fes.conn[i], [1.0/3 1.0/3])
        # Establish nodal triads
        _nodal_triads_e!(n_e, nvalid, e_g, normals, normal_valid, fes.conn[i])
        # Transform from global into nodal coordinates
        _transfmat_n_to_g!(Te, n_e, e_g)
        mul!(edisp_n, Te', edisp)
        # Now treat the transformation from the nodal to the element triad
        _transfmat_e_to_n!(Te, n_e, gradN_e)
         # Transform the nodal vector into the elementwise coordinates
        mul!(edisp_e, Te', edisp_n)
        updatecsmat!(outputcsys, loc, J0, fes.label[i]);
        if dot(view(outputcsys.csmat, :, 3), view(e_g, :, 3)) < 0.95
            @warn "Ordinate Systems Mismatched?"
        end
        o_e = e_g' * outputcsys.csmat
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
            _Bsmat!(Bs, ecoords_e)
            he = sqrt(2*Ae)
            shr = Bs * edisp_e
            frc = ((t^3/(t^2+0.2*he^2)))*Dt * shr
            fo = o2_e' * frc
            out[1:2] .= fo[1], fo[2]
        end
        # Call the inspector
        idat = inspector(idat, i, fes.conn[i], ecoords, out, loc);
        # end # Loop over quadrature points
    end # Loop over elements
    return idat; # return the updated inspector data
end

function inspectintegpoints(self::FEMMShellT3DSGA, geom::NodalField{FFlt},  u::NodalField{T}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {T<:Number, F<:Function}
    dT = NodalField(fill(zero(FFlt), nnodes(geom), 1)) # zero difference in temperature
    return inspectintegpoints(self, geom, u, dT, felist, inspector, idat, quantity; context...);
end


end # module
