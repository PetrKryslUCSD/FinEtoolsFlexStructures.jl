module FEMMShellDSG3IFModule

using LinearAlgebra: norm, Transpose, mul!, diag, eigen, I
using Statistics: mean
using FinEtools
import FinEtools.FESetModule: gradN!, nodesperelem, manifdim
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FEMMBaseModule: associategeometry!
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!, update!, thermalstrain!
using FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, locjac!, add_nnt_ut_only!, add_btsigma!
using ..FESetShellT3Module: FESetShellT3, local_frame!

const __nn = 3 # number of nodes
const __ndof = 6 # number of degrees of freedom per node

"""
    FEMMShellDSG3IF{S<:AbstractFESet, F<:Function} <: AbstractFEMM

Class for Discrete Shear Gap shell finite element modeling machine. With
averaging of the transverse strain-displacement matrix to provide isotropic
transverse shear response. Also, the formulation is developed to correctly
handle the coupling of twisting moments and transverse shear (such as in the
twisted beam or the Raasch hook problems).

Programming developed consistently with the paper
[1] Cui et al, Analysis of plates and shells using an edge-based smoothed finite
element method, Comput Mech (2010) 45:141–156
DOI 10.1007/s00466-009-0429-9

In this reference, the sign next to Ae in equation (44) is wrong.
[2] A superconvergent alpha finite element method (S a FEM) for static and
free vibration analysis of shell structures
Chai et al. (2017).

TO DO: 
- Take into account the possibility of a crease in the surface.
In that case the normal should not be averaged across the crease.
Along the crease every element should use the normal to its surface.
"""
mutable struct FEMMShellDSG3IF{S<:AbstractFESet, F<:Function, M} <: AbstractFEMM
    integdomain::IntegDomain{S, F} # integration domain data
    material::M # material object
    _associatedgeometry::Bool
    _normals::FFltMat
    # The attributes below are buffers used in various operations.
    _loc::FFltMat
    _J::FFltMat
    _J0::FFltMat
    _ecoords0::FFltMat
    _ecoords1::FFltMat
    _edisp1::FFltMat
    _evel1::FFltMat
    _evel1f::FFltMat
    _lecoords0::FFltMat
    _dofnums::FIntMat
    _F0::FFltMat
    _Tn::Vector{FFltMat} # transformation nodal-global matrices  
    _lTn::Vector{FFltMat} # transformation nodal-element matrices  
    _n::Vector{FFltVec} # nodal normals
    _ln::Vector{FFltVec} # nodal normals
    _Te::FFltMat
    _elmat::FFltMat
    _elmatTe::FFltMat
    _lloc::FFltMat
    _lJ::FFltMat
    _lgradN::FFltMat
    _Bm::FFltMat
    _Bb::FFltMat
    _Bs::FFltMat
    _DpsBmb::FFltMat
    _DtBs::FFltMat
end

function FEMMShellDSG3IF(integdomain::IntegDomain{S, F}, material::M) where {S<:AbstractFESet, F<:Function, M}
    _nnmax = 0
    for j in 1:count(integdomain.fes)
        for k in eachindex(integdomain.fes.conn[j])
            _nnmax = max(_nnmax, integdomain.fes.conn[j][k])
        end
    end
    _normals = fill(0.0, _nnmax, 3)
    # Allocate the buffers2
    _loc = fill(0.0, 1, 3)
    _J = fill(0.0, 3, 2)
    _J0 = fill(0.0, 3, 2)
    _ecoords0 = fill(0.0, __nn, 3); 
    _ecoords1 = fill(0.0, __nn, 3)
    _edisp1 = fill(0.0, __nn, 3); 
    _evel1 = fill(0.0, __nn, __ndof); 
    _evel1f = fill(0.0, __nn, __ndof)
    _lecoords0 = fill(0.0, __nn, 2) 
    _dofnums = zeros(FInt, 1, __nn*__ndof); 
    _F0 = fill(0.0, 3, 3); 
    _Tn = [fill(0.0, 3, 3), fill(0.0, 3, 3), fill(0.0, 3, 3)]; 
    _lTn = [fill(0.0, 3, 3), fill(0.0, 3, 3), fill(0.0, 3, 3)]; 
    _n = [fill(0.0, 3), fill(0.0, 3), fill(0.0, 3)]; 
    _ln = [fill(0.0, 3), fill(0.0, 3), fill(0.0, 3)]; 
    _Te = fill(0.0, __nn*__ndof, __nn*__ndof)
    _elmat = fill(0.0, __nn*__ndof, __nn*__ndof);    
    _elmatTe = fill(0.0, __nn*__ndof, __nn*__ndof); 
    _lloc = fill(0.0, 1, 2)
    _lJ = fill(0.0, 2, 2)
    _lgradN = fill(0.0, __nn, 2)
    _Bm = fill(0.0, 3, __nn*__ndof)
    _Bb = fill(0.0, 3, __nn*__ndof)
    _Bs = fill(0.0, 2, __nn*__ndof)
    _DpsBmb = similar(_Bm)
    _DtBs = similar(_Bs)
    
    return FEMMShellDSG3IF(integdomain, material,
        false,
        _normals,
        _loc, _J, _J0,
        _ecoords0, _ecoords1, _edisp1, _evel1, _evel1f, _lecoords0,
        _dofnums, 
        _F0, _Tn, _lTn, _n, _ln, _Te,
        _elmat, _elmatTe, 
        _lloc, _lJ, _lgradN,
        _Bm, _Bb, _Bs, _DpsBmb, _DtBs)
end

function make(integdomain, material)
    return FEMMShellDSG3IF(integdomain, material)
end

function _compute_J0!(J0, ecoords)
    x, y, z = ecoords[2, :].-ecoords[1, :]
    J0[:, 1] .= (x, y, z)
    x, y, z = ecoords[3, :].-ecoords[1, :]
    J0[:, 2] .= (x, y, z)
end
    
function _shell_material_stiffness(material)
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

function _gather_normals!(n, normals, c)
    for k in 1:length(n)
        n[k] .= vec(view(normals, c[k], :))
    end
    return n
end

function _compute_normals_e!(ln, n, F0)
    # The nodal normal is projected onto the vectors of the element cartesian
    # coordinate system. The normals are expressed using components on the
    # element basis.
    for k in 1:length(ln)
        ln[k] .= F0'*vec(view(n[k], :))
    end
    return ln
end

function _compute_nodal_triads_e!(lTn, ln)
    # Components of a cartesian ordinate system such that the third direction is
    # the direction of the nodal normal, and the angle to rotate one into the
    # other is as short as possible; these components are given on the local
    # element coordinate system basis vectors. 
    # These triads are the element-to-nodal transformation matrices.

    # TO DO Get rid of the temporaries
    r = similar(ln[1])
    f3_e = [0.0, 0.0, 1.0]
    for k in 1:length(lTn)
        cross3!(r, f3_e, ln[k])
        if norm(r) > 1.0e-6
            rotmat3!(lTn[k], r) 
        else
            lTn[k] .= I(3)
        end
    end
    return lTn
end

function _compute_nodal_triads_g!(Tn, lTn, F0)
    # The nodal triads expressed on the element basis are pushed forward with the
    # element rotation matrix.
    for k in 1:length(Tn)
        mul!(Tn[k], F0, lTn[k])
    end
    return Tn
end
 
function _g_n_transfmat!(Te, Tn, F0)
    # Global-to-nodal transformation matrix
    Te .= 0.0
    for i in 1:__nn
        offset = (i-1)*__ndof
        r = offset+1:offset+3
        @. Te[r, r] = Tn[i]
        r = offset+4:offset+6
        @. Te[r, r] = Tn[i]
    end
    return Te
end

function _n_e_transfmat!(Te, lTn, F0, lgradN)
    # Nodal-to-element transformation matrix. 
    # lTn = matrix with the nodal triad vectors in columns, components 
    # on the element basis
    # TO DO avoid a temporary
    Te .= 0.0
    # Translation degrees of freedom
    for i in 1:__nn
        roffset = (i-1)*__ndof
        r = roffset+1:roffset+3
        Te[r, r] .= lTn[i]' # this needs to be inverse
    end
    # Rotation degrees of freedom. The drilling rotation of the mid surface
    # produced by the 1/2*(v,x - u,y) effect is linked to the out of plane
    # rotations.
    for i in 1:__nn
        roffset = (i-1)*__ndof
        invlTn = inv(lTn[i][1:2, 1:2]') 
        Te[roffset.+(4:5), roffset.+(4:5)] .= invlTn
        Te[roffset+6, roffset+6] = 1/lTn[i][3, 3]
        r = invlTn * vec(lTn[i][1:2, 3]) # TO DO avoid the temporary
        for j in 1:__nn
            coffset = (j-1)*__ndof
            Te[coffset+1, roffset+4] = (-r[1] * 1/2 * lgradN[j, 2])
            Te[coffset+2, roffset+4] = (+r[1] * 1/2 * lgradN[j, 1])
            Te[coffset+1, roffset+5] = (-r[2] * 1/2 * lgradN[j, 2])
            Te[coffset+2, roffset+5] = (+r[2] * 1/2 * lgradN[j, 1])
        end
    end
    return Te
end


"""
    _Bsmat!(Bs, gradN, N)

Compute the linear transverse shear strain-displacement matrix.
"""
function _Bsmat!(Bs, lecoords)
    Bs .= 0.0 

    # Orientation 1
    s, p, q = 3, 1, 2
    a, b = lecoords[p, :] .- lecoords[s, :]
    c, d = lecoords[q, :] .- lecoords[s, :]
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
    a, b = lecoords[p, :] .- lecoords[s, :]
    c, d = lecoords[q, :] .- lecoords[s, :]

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
    a, b = lecoords[p, :] .- lecoords[s, :]
    c, d = lecoords[q, :] .- lecoords[s, :]
    
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

"""
    _Brmat!(Bb, gradN, N:)

Compute the linear, displacement independent, curvature-displacement/rotation matrix for a shell quadrilateral element with nfens=3 nodes. Displacements and rotations are in a local coordinate system.
"""
function _Brmat!(Br, gradN, N)
    for i in 1:__nn
        Br[1,6*(i-1)+1] = -1/2*gradN[i,2];
        Br[1,6*(i-1)+2] =  1/2*gradN[i,1];
        Br[1,6*(i-1)+6] =  -N[i];
    end
end

function associategeometry!(self::FEMMShellDSG3IF,  geom::NodalField{FFlt})
    J0 = self._J0
    F0 = self._F0
    normals = self._normals
    # Compute the normals at the nodes
    for el in 1:count(self.integdomain.fes)
        i, j, k = self.integdomain.fes.conn[el]
        J0[:, 1] = geom.values[j, :] - geom.values[i, :]
        J0[:, 2] = geom.values[k, :] - geom.values[i, :]
        local_frame!(delegateof(self.integdomain.fes), F0, J0)
        normals[i, :] .+= F0[:, 3]
        normals[j, :] .+= F0[:, 3]
        normals[k, :] .+= F0[:, 3]
    end
    for j in 1:size(normals, 1)
        nn = norm(normals[j, :])
        if nn > 0.0
            normals[j, :] ./= nn
        end
    end
    self._associatedgeometry = true
    return self
end

"""
    stiffness(self::FEMMShellDSG3IF, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the material stiffness matrix.
"""
function stiffness(self::FEMMShellDSG3IF, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    @assert self._associatedgeometry == true
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    normals = self._normals
    loc, J, J0 = self._loc, self._J, self._J0
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    lecoords0 = self._lecoords0
    F0, Tn, lTn, n, ln, Te = self._F0, self._Tn, self._lTn, self._n, self._ln, self._Te
    elmat, elmatTe = self._elmat, self._elmatTe
    lloc, lJ, lgradN = self._lloc, self._lJ, self._lgradN 
    Bm, Bb, Bs, DpsBmb, DtBs = self._Bm, self._Bb, self._Bs, self._DpsBmb, self._DtBs
    Br = fill(0.0, 1, size(elmat, 1))
    Dps, Dt = _shell_material_stiffness(self.material)
    scf=5/6;  # shear correction factor
    Dt .*= scf
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), dchi.nfreedofs, dchi.nfreedofs);
    for i in 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        _compute_J0!(J0, ecoords0)
        local_frame!(delegateof(fes), F0, J0)
        fill!(elmat,  0.0); # Initialize element matrix
        for j in 1:npts
            locjac!(loc, J, ecoords0, Ns[j], gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            t = self.integdomain.otherdimension(loc, fes.conn[i], Ns[j])
            mul!(lecoords0, ecoords0, view(F0, :, 1:2))
            locjac!(lloc, lJ, lecoords0, Ns[j], gradNparams[j])
            gradN!(fes, lgradN, gradNparams[j], lJ);
            _Bmmat!(Bm, lgradN)
            _Bbmat!(Bb, lgradN)
            _Bsmat!(Bs, lecoords0)
            add_btdb_ut_only!(elmat, Bm, t*Jac*w[j], Dps, DpsBmb)
            add_btdb_ut_only!(elmat, Bb, (t^3)/12*Jac*w[j], Dps, DpsBmb)
            # TO DO The stabilization expression has a significant effect
            # (at least for the pinched cylinder). What is the recommended
            # multiplier of he^2?
            he = sqrt(Jac)
            add_btdb_ut_only!(elmat, Bs, (t^3/(t^2+0.2*he^2))*Jac*w[j], Dt, DtBs)
            # add_btdb_ut_only!(elmat, Bs, t*Jac*w[j], Dt, DtBs)
        end
        # Complete the elementwise matrix by filling in the lower triangle
        complete_lt!(elmat)
        # Bending diagonal stiffness coefficients
        kavg4 = mean((elmat[4, 4], elmat[10, 10], elmat[16, 16]))
        kavg5 = mean((elmat[5, 5], elmat[11, 11], elmat[17, 17]))
        kavg = (kavg4 + kavg5) / 1e0
        # Add the artificial drilling stiffness
        elmat[6,6] += kavg
        elmat[12,12] += kavg
        elmat[18,18] += kavg
        # Now treat the transformation from the nodal to the element triad
        _gather_normals!(n, normals, fes.conn[i])
        _compute_normals_e!(ln, n, F0)
        _compute_nodal_triads_e!(lTn, ln)
        _n_e_transfmat!(Te, lTn, F0, lgradN)
        # Transform the elementwise matrix into the nodal coordinates
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        # Transform into global coordinates
        _compute_nodal_triads_g!(Tn, lTn, F0)
        _g_n_transfmat!(Te, Tn, F0)
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        # Assembly
        gatherdofnums!(dchi, dofnums, fes.conn[i]); 
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function stiffness(self::FEMMShellDSG3IF, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi);
end


"""
    mass(self::FEMMShellDSG3IF,  assembler::A,  geom::NodalField{FFlt}, dchi::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the consistent mass matrix

This is a general routine for the shell FEMM.
"""
function mass(self::FEMMShellDSG3IF,  assembler::A,  geom0::NodalField{FFlt}, dchi::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    loc, J, J0 = self._loc, self._J, self._J0
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    lecoords0 = self._lecoords0
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    
    elmat, elmatTe = self._elmat, self._elmatTe
    lloc, lJ, lgradN = self._lloc, self._lJ, self._lgradN 
    rho::FFlt = massdensity(self.material); # mass density
    npe = nodesperelem(fes)
    tmss = fill(0.0, npe);# basis f. matrix -- buffer
    rmss = fill(0.0, npe);# basis f. matrix -- buffer
    ndn = ndofs(dchi)
    startassembly!(assembler,  size(elmat,1),  size(elmat,2),  count(fes), dchi.nfreedofs,  dchi.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        _compute_J0!(J0, ecoords0)
        local_frame!(delegateof(fes), Ft, J0)
        fill!(tmss, 0.0)
        fill!(rmss, 0.0)
        # Compute the translational and rotational masses corresponding to nodes
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords0, Ns[j], gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            t = self.integdomain.otherdimension(loc, fes.conn[i], Ns[j])
            mul!(lecoords0, ecoords0, view(Ft, :, 1:2))
            tfactor = rho*(t*Jac*w[j]);
            rfactor = rho*(t^3/12*Jac*w[j]);
            for k in 1:npe
                tmss[k] += tfactor*(Ns[j][k])
                rmss[k] += rfactor*(Ns[j][k])
            end
        end # Loop over quadrature points
        fill!(elmat,  0.0); # Initialize element matrix
        for k in 1:npe
            for d in 1:3
                c = (k - 1) * __ndof + d
                elmat[c, c] += tmss[k]
            end
            for d in 4:5
                c = (k - 1) * __ndof + d
                elmat[c, c] += rmss[k]
            end
            d = 6
            c = (k - 1) * __ndof + d
            elmat[c, c] += rmss[k] / 1e6
        end
        # Transformation into global ordinates
        _transfmat!(Te, Ft)
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        # Assemble
        gatherdofnums!(dchi,  dofnums,  fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler,  elmat,  dofnums,  dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function mass(self::FEMMShellDSG3IF,  geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return mass(self, assembler, geom, u);
end

end # module
