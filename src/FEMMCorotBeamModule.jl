module FEMMCorotBeamModule

using LinearAlgebra: norm, Transpose, mul!
using FinEtools
using FinEtools.IntegDomainModule: IntegDomain
import FinEtoolsDeforLinear.MatDeforElastIsoModule: MatDeforElastIso
using ..FESetCorotBeamModule: FESetL2CorotBeam, local_frame_and_def!, local_mass!, local_stiffness!, natural_forces!, local_geometric_stiffness!, local_forces!, MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA


"""
    FEMMDeforLinear{MR<:AbstractDeforModelRed,  S<:AbstractFESet, F<:Function, M<:AbstractMatDeforLinearElastic} <: FEMMCorotBeam

Class for linear deformation finite element modeling machine.
"""
mutable struct FEMMCorotBeam{S<:AbstractFESet, F<:Function} <: AbstractFEMM
    integdomain::IntegDomain{S, F} # integration domain data
    material::MatDeforElastIso # material object
    # The attributes below are buffers used in various operations.
    _ecoords0::FFltMat
    _ecoords1::FFltMat
    _edisp1::FFltMat
    _evel1::FFltMat
    _evel1f::FFltMat
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
    _aN::FFltMat
    _dN::FFltVec
    _DN::FFltMat
    _PN::FFltVec
    _LF::FFltVec
    _RI::FFltMat
    _RJ::FFltMat
    _OS::FFltMat
end

function FEMMCorotBeam(integdomain::IntegDomain{S, F}, material::MatDeforElastIso) where {S<:FESetL2CorotBeam, F<:Function}
    _ecoords0 = fill(0.0, 2, 3); 
    _ecoords1 = fill(0.0, 2, 3)
    _edisp1 = fill(0.0, 2, 3); 
    _evel1 = fill(0.0, 2, 6); 
    _evel1f = fill(0.0, 2, 6)
    _dofnums = zeros(FInt, 1, 12); 
    _F0 = fill(0.0, 3, 3); 
    _Ft = fill(0.0, 3, 3); 
    _FtI = fill(0.0, 3, 3); 
    _FtJ = fill(0.0, 3, 3)
    _Te = fill(0.0, 12, 12)
    _tempelmat1 = fill(0.0, 12, 12); 
    _tempelmat2 = fill(0.0, 12, 12); 
    _tempelmat3 = fill(0.0, 12, 12)
    _elmat = fill(0.0, 12, 12);    
    _elmatTe = fill(0.0, 12, 12);    
    _elmato = fill(0.0, 12, 12)
    _elvec = fill(0.0, 12);    
    _elvecf = fill(0.0, 12)
    _aN = fill(0.0, 6, 12)
    _dN = fill(0.0, 6)
    _DN = fill(0.0, 6, 6)
    _PN = fill(0.0, 6)
    _LF = fill(0.0, 12)
    _RI = fill(0.0, 3, 3);    
    _RJ = fill(0.0, 3, 3);    
    _OS = fill(0.0, 3, 3)
    return FEMMCorotBeam(integdomain, material,
     _ecoords0, _ecoords1, _edisp1, _evel1, _evel1f, 
     _dofnums, 
     _F0, _Ft, _FtI, _FtJ, _Te,
     _tempelmat1, _tempelmat2, _tempelmat3, _elmat, _elmatTe, _elmato, 
     _elvec, _elvecf, 
     _aN, _dN, _DN, _PN, _LF, 
     _RI, _RJ, _OS)
end

function _transfmat!(Te, Ft)
    Te[1:3, 1:3] = Te[4:6, 4:6] = Te[7:9, 7:9] = Te[10:12, 10:12] = Ft
    return Te
end


"""
    mass(self::FEMMCorotBeam,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the consistent mass matrix

This is a general routine for the abstract linear-deformation  FEMM.
"""
function mass(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
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

function mass(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return mass(self, assembler, geom0, u1, Rfield1, dchi; mass_type = mass_type);
end

"""
    gyroscopic(self::FEMMCorotBeam,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the quadratic-inertial-term (gyroscopic) mass matrix

This is a general routine for the abstract linear-deformation  FEMM.
"""
function gyroscopic(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, v1::NodalField{T}, dchi::NodalField{TI}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    dN = self._dN
    evel1, evel1f = self._evel1, self._evel1f
    Ge1, Ge2, Ge = self._tempelmat1, self._tempelmat2, self._tempelmat3
    OmegaTilde = self._elmato
    OS = self._OS
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
        gathervalues_asmat!(v1, evel1, fes.conn[i]);
        evel1f[:, 1:3] = evel1[:,1:3]*Ft
        evel1f[:, 4:6] = evel1[:,4:6]*Ft
        Omega = (evel1f[1,4]+evel1f[2,4])/2*Ft[:,1] + (evel1f[1,3]-evel1f[2,3])/L1*Ft[:,2] + (evel1f[2,2]-evel1f[1,2])/L1*Ft[:,3];
        skewmat!(OS, Omega);
        _transfmat!(OmegaTilde, OS)
        mul!(Ge1, OmegaTilde, elmat)
        mul!(Ge2, elmat, OmegaTilde)
        @. Ge = Ge1 - Ge2;
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, Ge, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function gyroscopic(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, v1::NodalField{T}, dchi::NodalField{TI}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return gyroscopic(self, assembler, geom0, u1, Rfield1, v1, dchi; mass_type = mass_type);
end

"""
    stiffness(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the material stiffness matrix.
"""
function stiffness(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    aN, dN, DN = self._aN, self._dN, self._DN
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)::Float64
    A, I2, I3, J, A2s, A3s, x1x2_vector = fes.A, fes.I2, fes.I3, fes.J, fes.A2s, fes.A3s, fes.x1x2_vector
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
        local_stiffness!(elmat, E, G, A[i], I2[i], I3[i], J[i], A2s[i], A3s[i], L1, aN, DN);
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function stiffness(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi);
end


"""
    geostiffness(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the geometric stiffness matrix.
"""
function geostiffness(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    aN, dN, DN, PN = self._aN, self._dN, self._DN, self._PN
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)
    A, I2, I3, J, A2s, A3s, x1x2_vector = fes.A, fes.I2, fes.I3, fes.J, fes.A2s, fes.A3s, fes.x1x2_vector
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
        natural_forces!(PN, E, G, A[i], I2[i], I3[i], J[i], A2s[i], A3s[i], L1, dN, DN)
        local_geometric_stiffness!(elmat, A[i], I2[i], I3[i], PN, L1)
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function geostiffness(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return geostiffness(self, assembler, geom0, u1, Rfield1, dchi);
end

"""
    restoringforce(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the vector of the restoring elastic forces
"""
function restoringforce(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysvecAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    aN, dN, DN, PN, LF = self._aN, self._dN, self._DN, self._PN, self._LF
    elvec = self._elvec
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)
    A, I2, I3, J, A2s, A3s, x1x2_vector = fes.A, fes.I2, fes.I3, fes.J, fes.A2s, fes.A3s, fes.x1x2_vector
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
        natural_forces!(PN, E, G, A[i], I2[i], I3[i], J[i], A2s[i], A3s[i], L1, dN, DN)
        local_forces!(LF, PN, L1, aN)
        mul!(elvec, Te, -LF)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elvec, dofnums); 
    end # Loop over elements
    return makevector!(assembler);
end

function restoringforce(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {T<:Number, TI<:Number}
    assembler = SysvecAssembler();
    return restoringforce(self, assembler, geom0, u1, Rfield1, dchi);
end

"""
    distribloads_global(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}, fi) where {T<:Number}
    
Compute the load vector due to distributed loads.

Compute the global load vector corresponding to applied distributed
load. Here it means force per unit length of the beam,
in the configuration u1,Rfield1.

Note: the force intensity must be uniform across the entire element.
Note: the force intensity is given in the global coordinates.
"""
function distribloads_global(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}, fi) where {ASS<:AbstractSysvecAssembler, T<:Number, TI<:Number}
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

function distribloads_global(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}, fi) where {T<:Number, TI<:Number}
    assembler = SysvecAssembler();
    return distribloads_global(self, assembler, geom0, u1, Rfield1, dchi, fi);
end

end # module
