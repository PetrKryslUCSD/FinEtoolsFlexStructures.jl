module FEMMShellDSG3Module

using LinearAlgebra: norm, Transpose, mul!
using FinEtools
import FinEtools.FESetModule: gradN!, nodesperelem, manifdim
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!, update!, thermalstrain!
using FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, locjac!, add_nnt_ut_only!, add_btsigma!
using ..FESetShellDSG3Module: FESetShellDSG3, local_frame!

const __nn = 3 # number of nodes
const __ndof = 6 # number of degrees of freedom per node

"""
    FEMMShellDSG3{S<:AbstractFESet, F<:Function} <: AbstractFEMM

Class for Discrete Shear Gap shell finite element modeling machine.

Programming developed consistently with the paper
[1] A superconvergent alpha finite element method (S a FEM) for static and
free vibration analysis of shell structures
Chai et al. (2017).
"""
mutable struct FEMMShellDSG3{S<:AbstractFESet, F<:Function, M} <: AbstractFEMM
    integdomain::IntegDomain{S, F} # integration domain data
    material::M # material object
    # The attributes below are buffers used in various operations.
    _loc::FFltMat
    _J::FFltMat
    _ecoords0::FFltMat
    _ecoords1::FFltMat
    _edisp1::FFltMat
    _evel1::FFltMat
    _evel1f::FFltMat
    _lecoords0::FFltMat
    _dofnums::FIntMat
    _F0::FFltMat
    _Ft::FFltMat
    _FtI::FFltMat
    _FtJ::FFltMat
    _Te::FFltMat
    _tempelmat1::FFltMat
    _tempelmat2::FFltMat
    _tempelmat3::FFltMat
    _elmat::FFltMat
    _elmatTe::FFltMat
    _elmato::FFltMat
    _elvec::FFltVec
    _elvecf::FFltVec
    _lloc::FFltMat
    _lJ::FFltMat
    _lgradN::FFltMat
    _Bm::FFltMat
    _Bb::FFltMat
    _Bs::FFltMat
    _DpsBmb::FFltMat
    _DtBs::FFltMat
    _LF::FFltVec
    _RI::FFltMat
    _RJ::FFltMat
    _OS::FFltMat
end

function FEMMShellDSG3(integdomain::IntegDomain{S, F}, material::M) where {S<:AbstractFESet, F<:Function, M}
    _loc = fill(0.0, 1, 3)
    _J = fill(0.0, 3, 2)
    _ecoords0 = fill(0.0, __nn, 3); 
    _ecoords1 = fill(0.0, __nn, 3)
    _edisp1 = fill(0.0, __nn, 3); 
    _evel1 = fill(0.0, __nn, __ndof); 
    _evel1f = fill(0.0, __nn, __ndof)
    _lecoords0 = fill(0.0, __nn, 2) 
    _dofnums = zeros(FInt, 1, __nn*__ndof); 
    _F0 = fill(0.0, 3, 3); 
    _Ft = fill(0.0, 3, 3); 
    _FtI = fill(0.0, 3, 3); 
    _FtJ = fill(0.0, 3, 3)
    _Te = fill(0.0, __nn*__ndof, __nn*__ndof)
    _tempelmat1 = fill(0.0, __nn*__ndof, __nn*__ndof); 
    _tempelmat2 = fill(0.0, __nn*__ndof, __nn*__ndof); 
    _tempelmat3 = fill(0.0, __nn*__ndof, __nn*__ndof)
    _elmat = fill(0.0, __nn*__ndof, __nn*__ndof);    
    _elmatTe = fill(0.0, __nn*__ndof, __nn*__ndof);    
    _elmato = fill(0.0, __nn*__ndof, __nn*__ndof)
    _elvec = fill(0.0, __nn*__ndof);    
    _elvecf = fill(0.0, __nn*__ndof)
    _lloc = fill(0.0, 1, 2)
    _lJ = fill(0.0, 2, 2)
    _lgradN = fill(0.0, __nn, 2)
    _Bm = fill(0.0, 3, __nn*__ndof)
    _Bb = fill(0.0, 3, __nn*__ndof)
    _Bs = fill(0.0, 2, __nn*__ndof)
    _DpsBmb = similar(_Bm)
    _DtBs = similar(_Bs)
    _LF = fill(0.0, __nn*__ndof)
    _RI = fill(0.0, 3, 3);    
    _RJ = fill(0.0, 3, 3);    
    _OS = fill(0.0, 3, 3)
    return FEMMShellDSG3(integdomain, material,
        _loc, _J,
        _ecoords0, _ecoords1, _edisp1, _evel1, _evel1f, _lecoords0,
        _dofnums, 
        _F0, _Ft, _FtI, _FtJ, _Te,
        _tempelmat1, _tempelmat2, _tempelmat3, _elmat, _elmatTe, _elmato, 
        _elvec, _elvecf, 
        _lloc, _lJ, _lgradN,
        _Bm, _Bb, _Bs, _DpsBmb, _DtBs, _LF, 
        _RI, _RJ, _OS)
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

function _transfmat!(Te, Ft)
    for i in 1:2*__nn
        r = (i-1)*__nn .+ (1:3)
        @. Te[r, r] = Ft
    end
    return Te
end

"""
    _Bsmat!(Bs, gradN, N)

Compute the linear transverse shear strain-displacement matrix.
"""
function _Bsmat!(Bs, gradN, N)
    for i in 1:__nn
        Bs[1,6*(i-1)+3] = gradN[i,1];
        Bs[1,6*(i-1)+5] = N[i];
        Bs[2,6*(i-1)+3] = gradN[i,2];
        Bs[2,6*(i-1)+4] = -N[i];
    end
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
    stiffness(self::FEMMShellDSG3, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the material stiffness matrix.
"""
function stiffness(self::FEMMShellDSG3, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    loc, J = self._loc, self._J
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    lecoords0 = self._lecoords0
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    lloc, lJ, lgradN = self._lloc, self._lJ, self._lgradN 
    Bm, Bb, Bs, DpsBmb, DtBs = self._Bm, self._Bb, self._Bs, self._DpsBmb, self._DtBs
   Dps, Dt = _shell_material_stiffness(self.material)
    scf=5/6;  # shear correction factor
    Dt .*= scf
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), dchi.nfreedofs, dchi.nfreedofs);
    for i in 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        fill!(elmat,  0.0); # Initialize element matrix
        for j in 1:npts
            locjac!(loc, J, ecoords0, Ns[j], gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            t = self.integdomain.otherdimension(loc, fes.conn[i], Ns[j])
            local_frame!(delegateof(fes), Ft, J)
            mul!(lecoords0, ecoords0, view(Ft, :, 1:2))
            locjac!(lloc, lJ, lecoords0, Ns[j], gradNparams[j])
            gradN!(fes, lgradN, gradNparams[j], lJ);
            _Bmmat!(Bm, lgradN)
            _Bbmat!(Bb, lgradN)
            _Bsmat!(Bs, lgradN, Ns[j])
            add_btdb_ut_only!(elmat, Bm, t*Jac*w[j], Dps, DpsBmb)
            add_btdb_ut_only!(elmat, Bb, (t^3)/12*Jac*w[j], Dps, DpsBmb)
            add_btdb_ut_only!(elmat, Bs, t*Jac*w[j], Dt, DtBs)
            complete_lt!(elmat)
            _transfmat!(Te, Ft)
            mul!(elmatTe, elmat, Transpose(Te))
            mul!(elmat, Te, elmatTe)
        end
        gatherdofnums!(dchi, dofnums, fes.conn[i]); 
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function stiffness(self::FEMMShellDSG3, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi);
end


    # % Loop over all integration points: bending and membrane contributions
    # for j=1:npts_per_gcell_bm
    #     Nder = Ndermat_param (gcells(i), pc_bm(j,:));
    #     ts=Nder'*x; % rows: tangent vectors to parametric curves
    #     es=local_basis(self,ts); % calculate the local Cartesian basis
    #     xl=x*es(:,1:2);% local coordinates
    #     [Nspatialder,detJ] = Ndermat_spatial (gcells(i), Nder, xl);
    #     D = tangent_moduli (mat, struct ('ms',matstates{i,j}));% note: the material must be plane stress
    #     h =thickness(gcells(i), pc_bm(j,:));
    #     Bb = Bblmat (gcells(i), Nspatialder);
    #     Bm = Bmlmat (gcells(i), Nspatialder);
    #     Kl = ((h^3/15)*Bb'*D*Bb + h*Bm'*D*Bm) * detJ * w_bm(j);
    #     T = zeros(2*dim*nfens); for q=1:2*nfens, s=3*(q-1)+1; T(s:s+2,s:s+2)=es; end
    #     Ke = Ke + T*Kl*T';
    #     k=k+Kl;
    # end
    # Kbavg=sum(diag(Ke(4:6:2*dim*nfens,4:6:2*dim*nfens))+diag(Ke(5:6:2*dim*nfens,5:6:2*dim*nfens)))/nfens/2;
    # % Loop over all integration points: out of plane shear
    # for j=1:npts_per_gcell_s
    #     Nder = Ndermat_param (gcells(i), pc_s(j,:));
    #     ts=Nder'*x; % rows: tangent vectors to parametric curves
    #     es=local_basis(self,ts); % calculate the local Cartesian basis
    #     xl=x*es(:,1:2);% local coordinates
    #     [Nspatialder,detJ] = Ndermat_spatial (gcells(i), Nder, xl);
    #     h =thickness(gcells(i), pc_s(j,:));
    #     G=get(mat,'E')/2/(1.0+get(mat,'nu')); % shear modulus: raw -- need a better way to account for anisotropic materials
    #     Ds=G*scf*eye(2);  % shear stiffness
    #     Bs = Bslmat (gcells(i), Nspatialder, Nmat(gcells(i), pc_s(j,:)));
    #     Kl = (Bs'*Ds*Bs) * h * detJ * w_s(j);
    #     T = zeros(2*dim*nfens); for q=1:2*nfens, s=3*(q-1)+1; T(s:s+2,s:s+2)=es; end
#         % add a very small stiffness (compared to bending) to stabilize
#         % twisting deformation about the normal to the surface
#         Kl(6:6:2*dim*nfens,6:6:2*dim*nfens)=Kl(6:6:2*dim*nfens,6:6:2*dim*nfens)+10e-7*Kbavg;
#         Ke = Ke + T*Kl*T';
#         k=k+Kl;
#     end
#     ems(i) = set(ems(i), 'mat', Ke);
#     ems(i) = set(ems(i), 'eqnums', gather(u, conn, 'eqnums'));
# end
# return;


"""
    distribloads_global(self::FEMMShellDSG3, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}, fi) where {T<:Number}
    
Compute the load vector due to distributed loads.

Compute the global load vector corresponding to applied distributed
load. Here it means force per unit length of the beam,
in the configuration u1,Rfield1.

Note: the force intensity must be uniform across the entire element.
Note: the force intensity is given in the global coordinates.
"""
function distribloads_global(self::FEMMShellDSG3, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}, fi) where {ASS<:AbstractSysvecAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    aN, dN, DN, PN = self._aN, self._dN, self._DN, self._PN
    elvec, elvecf = self._elvec, self._elvecf
    Lforce = fill(0.0, 3)
    ignore = fill(0.0 , 0, 0)
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)
    A, I2, I3, J, x1x2_vector = fes.A, fes.I2, fes.I3, fes.J, fes.x1x2_vector
    startassembly!(assembler, dchi.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        R1I[:] .= Rfield1.values[fes.conn[i][1], :];
        R1J[:] .= Rfield1.values[fes.conn[i][2], :];
        fill!(elmat,  0.0); # Initialize element matrix
        L1, Ft, dN = local_frame_and_def!(Ft, dN, F0, FtI, FtJ, ecoords0, x1x2_vector[i], ecoords1, R1I, R1J);
        _transfmat!(Te, Ft)
        L0 = norm(ecoords0[2,:]-ecoords0[1,:]); 
        force = updateforce!(fi, ignore, ignore, fes.label[i]); # retrieve the applied load
        Lforce = Ft' * force # local distributed load components
        elvecf[1] = Lforce[1]*L0/2;
        elvecf[2] = Lforce[2]*L0/2;
        elvecf[3] = Lforce[3]*L0/2;
        elvecf[4] = 0;
        elvecf[5] = -Lforce[3]*L0^2/12;
        elvecf[6] = +Lforce[2]*L0^2/12;
        elvecf[7] = Lforce[1]*L0/2;
        elvecf[8] = Lforce[2]*L0/2;
        elvecf[9] = Lforce[3]*L0/2;
        elvecf[10] = 0;
        elvecf[11] = +Lforce[3]*L0^2/12;
        elvecf[12] = -Lforce[2]*L0^2/12;
        mul!(elvec, Te, elvecf)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elvec, dofnums); 
    end # Loop over elements
    return makevector!(assembler);
end

function distribloads_global(self::FEMMShellDSG3, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}, fi) where {T<:Number, TI<:Number}
    assembler = SysvecAssembler();
    return distribloads_global(self, assembler, geom0, u1, Rfield1, dchi, fi);
end

"""
    mass(self::FEMMShellDSG3,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the consistent mass matrix

This is a general routine for the abstract linear-deformation  FEMM.
"""
function mass(self::FEMMShellDSG3, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    dN = self._dN
    rho = massdensity(self.material)
    A, I1, I2, I3, x1x2_vector = fes.A, fes.I1, fes.I2, fes.I3, fes.x1x2_vector
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), dchi.nfreedofs, dchi.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        R1I[:] .= Rfield1.values[fes.conn[i][1], :];
        R1J[:] .= Rfield1.values[fes.conn[i][2], :];
        fill!(elmat,  0.0); # Initialize element matrix
        L1, Ft, dN = local_frame_and_def!(Ft, dN, F0, FtI, FtJ, ecoords0, x1x2_vector[i], ecoords1, R1I, R1J);
        _transfmat!(Te, Ft)
        L0 = norm(ecoords0[2,:]-ecoords0[1,:]); 
        local_mass!(elmat, A[i], I1[i], I2[i], I3[i], rho, L0, mass_type);
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function mass(self::FEMMShellDSG3, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return mass(self, assembler, geom0, u1, Rfield1, dchi; mass_type = mass_type);
end

end # module
