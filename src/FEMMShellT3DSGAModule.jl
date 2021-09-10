module FEMMShellT3DSGAModule

using LinearAlgebra: norm, Transpose, mul!, diag, eigen, I, dot
using Statistics: mean
using FinEtools
import FinEtools.FESetModule: gradN!, nodesperelem, manifdim
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FEMMBaseModule: associategeometry!
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!, update!, thermalstrain!
using FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, locjac!, add_nnt_ut_only!, add_btsigma!
using ..FESetShellT3Module: FESetShellT3, local_frame!
using ..TransformerModule: Transformer

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
element method, Comput Mech (2010) 45:141–156
DOI 10.1007/s00466-009-0429-9

In this reference, the sign next to Ae in equation (44) is wrong:
[2] A superconvergent alpha finite element method (S a FEM) for static and
free vibration analysis of shell structures
Chai et al. (2017).

The treatment of the transformation between the element and nodal coordinates
is carried out using an alternative to the publication

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
    material::M # material object.
    _associatedgeometry::Bool
    _normals::FFltMat
    _normal_valid::Vector{Bool}
    # The attributes below are buffers used in various operations.
    _loc::FFltMat
    _J::FFltMat
    _J0::FFltMat
    _ecoords::FFltMat
    _edisp::FFltMat
    _ecoords_e::FFltMat
    _dofnums::FIntMat
    _e_g::FFltMat
    _n_e::Vector{FFltMat} # transformation nodal-element matrices  
    _nvalid::Vector{Bool}
    _Te::FFltMat
    _elmat::FFltMat
    _loc_e::FFltMat
    _J_e::FFltMat
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
    # Alocate_e the buffers2
    _loc = fill(0.0, 1, 3)
    _J = fill(0.0, 3, 2)
    _J0 = fill(0.0, 3, 2)
    _ecoords = fill(0.0, __nn, 3)
    _edisp = fill(0.0, __nn, 6); 
    _ecoords_e = fill(0.0, __nn, 2) 
    _dofnums = zeros(FInt, 1, __nn*__ndof); 
    _e_g = fill(0.0, 3, 3); 
    _n_e = [fill(0.0, 3, 3), fill(0.0, 3, 3), fill(0.0, 3, 3)]; 
    _nvalid = fill(false, 3)
    _Te = fill(0.0, __nn*__ndof, __nn*__ndof)
    _elmat = fill(0.0, __nn*__ndof, __nn*__ndof);   
    _loc_e = fill(0.0, 1, 2)
    _J_e = fill(0.0, 2, 2)
    _gradN_e = fill(0.0, __nn, 2)
    _Bm = fill(0.0, 3, __nn*__ndof)
    _Bb = fill(0.0, 3, __nn*__ndof)
    _Bs = fill(0.0, 2, __nn*__ndof)
    _DpsBmb = similar(_Bm)
    _DtBs = similar(_Bs)
    
    return FEMMShellT3DSGA(integdomain, material,
        false,
        _normals, _normal_valid,
        _loc, _J, _J0,
        _ecoords, _edisp, _ecoords_e,
        _dofnums, 
        _e_g, _n_e, _nvalid, _Te,
        _elmat, 
        _loc_e, _J_e, _gradN_e,
        _Bm, _Bb, _Bs, _DpsBmb, _DtBs)
end

function make(integdomain, material)
    return FEMMShellT3DSGA(integdomain, material)
end

function _compute_J0!(J0, ecoords)
    # Compute the Jacobian: the Jacobian is constant within the triangle
    x, y, z = ecoords[2, :].-ecoords[1, :]
    J0[:, 1] .= (x, y, z)
    x, y, z = ecoords[3, :].-ecoords[1, :]
    J0[:, 2] .= (x, y, z)
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
        Te[r, r] .= n_eT # this needs to be inverse
        r = roffset+4:roffset+6
        Te[r, r] .= n_eT # this needs to be inverse
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

    # The first node in the triangle 
    m = (1/2/Ae) * (1/3) # multiplier
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
    m = (1/2/Ae) * (1/3) # multiplier
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
    m = (1/2/Ae) * (1/3) # multiplier
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
    for i in 1:__nn
        Bm[1,6*(i-1)+1] = gradN[i,1];
        Bm[2,6*(i-1)+2] = gradN[i,2];
        Bm[3,6*(i-1)+1] = gradN[i,2];
        Bm[3,6*(i-1)+2] = gradN[i,1];
    end
end

"""
    _Bbmat!(Bb, gradN)

Compute the linear, displacement independent, curvature-displacement/rotation matrix for a shell quadrilateral element with nfens=3 nodes. Displacements and rotations are in a local coordinate system.
"""
function _Bbmat!(Bb, gradN)
    for i in 1:__nn
        Bb[1,6*(i-1)+5] = gradN[i,1];
        Bb[2,6*(i-1)+4] = -gradN[i,2];
        Bb[3,6*(i-1)+4] = -gradN[i,1];
        Bb[3,6*(i-1)+5] = gradN[i,2];
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
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    normals, normal_valid = self._normals, self._normal_valid
    loc, J, J0 = self._loc, self._J, self._J0
    ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
    ecoords_e = self._ecoords_e
    e_g, n_e, nvalid, Te = self._e_g, self._n_e, self._nvalid, self._Te
    elmat = self._elmat
    transformwith = Transformer(elmat)
    loc_e, J_e, gradN_e = self._loc_e, self._J_e, self._gradN_e 
    Bm, Bb, Bs, DpsBmb, DtBs = self._Bm, self._Bb, self._Bs, self._DpsBmb, self._DtBs
    Br = fill(0.0, 1, size(elmat, 1))
    Dps, Dt = _shell_material_stiffness(self.material)
    scf=5/6;  # shear correction factor
    Dt .*= scf
    drilling_stiffness_scale = 1.0e0
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), dchi.nfreedofs, dchi.nfreedofs);
    for i in 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords, fes.conn[i]);
        _compute_J0!(J0, ecoords)
        local_frame!(delegateof(fes), e_g, J0)
        fill!(elmat,  0.0); # Initialize element matrix
        for j in 1:npts
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            t = self.integdomain.otherdimension(loc, fes.conn[i], Ns[j])
            mul!(ecoords_e, ecoords, view(e_g, :, 1:2))
            locjac!(loc_e, J_e, ecoords_e, Ns[j], gradNparams[j])
            gradN!(fes, gradN_e, gradNparams[j], J_e);
            _Bmmat!(Bm, gradN_e)
            _Bbmat!(Bb, gradN_e)
            _Bsmat!(Bs, ecoords_e)
            add_btdb_ut_only!(elmat, Bm, t*Jac*w[j], Dps, DpsBmb)
            add_btdb_ut_only!(elmat, Bb, (t^3)/12*Jac*w[j], Dps, DpsBmb)
            # TO DO The stabilization expression has a significant effect
            # (at least for the pinched cylinder). What is the recommended
            # multiplier of he^2?
            he = sqrt(Jac)
            add_btdb_ut_only!(elmat, Bs, (t^3/(t^2+0.2*he^2))*Jac*w[j], Dt, DtBs)
        end
        # Complete the elementwise matrix by filling in the lower triangle
        complete_lt!(elmat)
        # Now treat the transformation from the element to the nodal triad
        _nodal_triads_e!(n_e, nvalid, e_g, normals, normal_valid, fes.conn[i])
        _transfmat_e_to_n!(Te, n_e, gradN_e)
         # Transform the elementwise matrix into the nodal coordinates
        transformwith(elmat, Te)
        # Bending diagonal stiffness coefficients
        kavg4 = mean((elmat[4, 4], elmat[10, 10], elmat[16, 16]))
        kavg5 = mean((elmat[5, 5], elmat[11, 11], elmat[17, 17]))
        kavg = (kavg4 + kavg5) * drilling_stiffness_scale
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
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    normals, normal_valid = self._normals, self._normal_valid
    loc, J, J0 = self._loc, self._J, self._J0
    ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
    ecoords_e = self._ecoords_e
    e_g, n_e, nvalid, Te = self._e_g, self._n_e, self._nvalid, self._Te
    elmat = self._elmat
    transformwith = Transformer(elmat)
    loc_e, J_e, gradN_e = self._loc_e, self._J_e, self._gradN_e 
    rho::FFlt = massdensity(self.material); # mass density
    npe = nodesperelem(fes)
    tmss = fill(0.0, npe);# basis f. matrix -- buffer
    rmss = fill(0.0, npe);# basis f. matrix -- buffer
    ndn = ndofs(dchi)
    startassembly!(assembler,  size(elmat,1),  size(elmat,2),  count(fes), dchi.nfreedofs,  dchi.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords, fes.conn[i]);
        _compute_J0!(J0, ecoords)
        local_frame!(delegateof(fes), e_g, J0)
        fill!(tmss, 0.0)
        fill!(rmss, 0.0)
        # Compute the translational and rotational masses corresponding to nodes
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            t = self.integdomain.otherdimension(loc, fes.conn[i], Ns[j])
            mul!(ecoords_e, ecoords, view(e_g, :, 1:2))
            tfactor = rho*(t*Jac*w[j]);
            rfactor = rho*(t^3/12*Jac*w[j]);
            for k in 1:npe
                tmss[k] += tfactor*(Ns[j][k])
                rmss[k] += rfactor*(Ns[j][k])
            end
        end # Loop over quadrature points
        # Now treat the transformation from the element to the nodal triad
        _nodal_triads_e!(n_e, nvalid, e_g, normals, normal_valid, fes.conn[i])
        # Fill the elementwise matrix in the nodal basis: the matrix is lumped,
        # hence diagonal, and no transformation from element to nodal is needed.
        fill!(elmat,  0.0); # Initialize element matrix
        for k in 1:npe
            # Translation degrees of freedom
            for d in 1:3
                c = (k - 1) * __ndof + d
                elmat[c, c] += tmss[k]
            end
            # Bending degrees of freedom
            for d in 4:5
                c = (k - 1) * __ndof + d
                elmat[c, c] += rmss[k]
            end
            # Drilling rotations
            if nvalid[k]
                d = 6
                c = (k - 1) * __ndof + d
                elmat[c, c] += rmss[k] / 1e6
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
function inspectintegpoints(self::FEMMShellT3DSGA, geom::NodalField{FFlt},  u::NodalField{T}, dT::NodalField{FFlt}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {T<:Number, F<:Function}
    fes = self.integdomain.fes
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    normals, normal_valid = self._normals, self._normal_valid
    loc, J, J0 = self._loc, self._J, self._J0
    ecoords, edisp, dofnums = self._ecoords, self._edisp, self._dofnums
    ecoords_e = self._ecoords_e
    e_g, n_e, n, ln, nvalid, Te = self._e_g, self._n_e, self._n, self._nv_e, self._nvalid, self._Te
    elmat = self._elmat
    loc_e, J_e, gradN_e = self._loc_e, self._J_e, self._gradN_e 
    Bm, Bb, Bs, DpsBmb, DtBs = self._Bm, self._Bb, self._Bs, self._DpsBmb, self._DtBs
    Br = fill(0.0, 1, size(elmat, 1))
    Dps, Dt = _shell_material_stiffness(self.material)
    scf=5/6;  # shear correction factor
    Dt .*= scf
    # Sort out  the output requirements
    outputcsys = self.mcsys; # default: report the stresses in the material coord system
    for apair in pairs(context)
        sy, val = apair
        if sy == :outputcsys
            outputcsys = val
        end
    end
    # t= 0.0
    # dt = 0.0
    # dTe = fill(zero(FFlt), nodesperelem(fes)) # nodal temperatures -- buffer
    # ue = fill(zero(FFlt), size(elmat, 1)); # array of node displacements -- buffer
    # nne = nodesperelem(fes); # number of nodes for element
    # sdim = ndofs(geom);            # number of space dimensions
    # xe = fill(zero(FFlt), nne, sdim); # array of node coordinates -- buffer
    # qpdT = 0.0; # node temperature increment
    # qpstrain = fill(zero(FFlt), nstressstrain(self.mr), 1); # total strain -- buffer
    # qpthstrain = fill(zero(FFlt), nthermstrain(self.mr)); # thermal strain -- buffer
    # qpstress = fill(zero(FFlt), nstressstrain(self.mr)); # stress -- buffer
    # out1 = fill(zero(FFlt), nstressstrain(self.mr)); # stress -- buffer
    # out =  fill(zero(FFlt), nstressstrain(self.mr));# output -- buffer
    # Loop over  all the elements and all the quadrature points within them
    for ilist = 1:length(felist) # Loop over elements
        i = felist[ilist];
        gathervalues_asmat!(geom0, ecoords, fes.conn[i]);
        gathervalues_asmat!(u1, edisp, fes.conn[i]);
        ecoords1 .= ecoords .+ edisp
        _compute_J0!(J0, ecoords)
        local_frame!(delegateof(fes), e_g, J0)
        fill!(elmat,  0.0); # Initialize element matrix
        for j in 1:npts
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            t = self.integdomain.otherdimension(loc, fes.conn[i], Ns[j])
            mul!(ecoords_e, ecoords, view(e_g, :, 1:2))
            locjac!(loc_e, J_e, ecoords_e, Ns[j], gradNparams[j])
            gradN!(fes, gradN_e, gradNparams[j], J_e);
            _Bmmat!(Bm, gradN_e)
            _Bbmat!(Bb, gradN_e)
            _Bsmat!(Bs, ecoords_e)
            # add_btdb_ut_only!(elmat, Bm, t*Jac*w[j], Dps, DpsBmb)
            # add_btdb_ut_only!(elmat, Bb, (t^3)/12*Jac*w[j], Dps, DpsBmb)
            # TO DO The stabilization expression has a significant effect
            # (at least for the pinched cylinder). What is the recommended
            # multiplier of he^2?
            he = sqrt(Jac)
            # add_btdb_ut_only!(elmat, Bs, (t^3/(t^2+0.2*he^2))*Jac*w[j], Dt, DtBs)
            updatecsmat!(outputcsys, loc, J, fes.label[i]);
            # Call the inspector
            #     idat = inspector(idat, i, fes.conn[i], ecoords, out, loc);
        end # Loop over quadrature points

        # gathervalues_asmat!(geom, ecoords, fes.conn[i]);
        # gathervalues_asvec!(u, ue, fes.conn[i]);# retrieve element displacements
        # # gathervalues_asvec!(dT, dTe, fes.conn[i]);# retrieve element temp. increments
        # for j = 1:npts # Loop over quadrature points
        #     locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
        #     Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
        #     updatecsmat!(self.mcsys, loc, J, fes.label[i]);
        #     At_mul_B!(csmatTJ,  self.mcsys.csmat,  J); # local Jacobian matrix
        #     gradN!(fes, gradN, gradNparams[j], csmatTJ);
        #     Blmat!(self.mr, B, Ns[j], gradN, loc, self.mcsys.csmat);
        #     updatecsmat!(outputcsys, loc, J, fes.label[i]);
        #     # Quadrature point quantities
        #     A_mul_B!(qpstrain, B, ue); # strain in material coordinates
        #     qpdT = dot(vec(dTe), vec(Ns[j]));# Quadrature point temperature increment
        #     thermalstrain!(self.material, qpthstrain, qpdT)
        #     # Material updates the state and returns the output
        #     out = update!(self.material, qpstress, out, vec(qpstrain), qpthstrain, t, dt, loc, fes.label[i], quantity)
        #     if (quantity == :Cauchy)   # Transform stress tensor,  if that is "out"
        #         (length(out1) >= length(out)) || (out1 = zeros(length(out)))
        #         rotstressvec!(self.mr, out1, out, transpose(self.mcsys.csmat))# To global coord sys
        #         rotstressvec!(self.mr, out, out1, outputcsys.csmat)# To output coord sys
        #     end
        #     # Call the inspector
        #     idat = inspector(idat, i, fes.conn[i], ecoords, out, loc);
        # end # Loop over quadrature points
    end # Loop over elements
    return idat; # return the updated inspector data
end

function inspectintegpoints(self::FEMMShellT3DSGA, geom::NodalField{FFlt},  u::NodalField{T}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {T<:Number, F<:Function}
    dT = NodalField(fill(zero(FFlt), nnodes(geom), 1)) # zero difference in temperature
    return inspectintegpoints(self, geom, u, dT, felist, inspector, idat, quantity; context...);
end


end # module
