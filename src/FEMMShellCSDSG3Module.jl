module FEMMShellCSDSG3Module

using LinearAlgebra: norm, Transpose, mul!, diag, rank, eigen
using Statistics: mean
using FinEtools
import FinEtools.FESetModule: gradN!, nodesperelem, manifdim
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!, update!, thermalstrain!
using FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, locjac!, add_nnt_ut_only!, add_btsigma!
using ..FESetShellT3Module: FESetShellT3, local_frame!

using Infiltrator

const __nn = 3 # number of nodes
const __ndof = 6 # number of degrees of freedom per node

"""
    FEMMShellCSDSG3{S<:AbstractFESet, F<:Function} <: AbstractFEMM

Class for Cell-Smoothed Discrete Shear Gap shell finite element modeling machine.

Programming developed consistently with the paper
[1] Cui et al, Analysis of plates and shells using an edge-based smoothed finite
element method, Comput Mech (2010) 45:141–156
DOI 10.1007/s00466-009-0429-9

In this reference, the sign next to Ae in equation (44) is wrong.
[2] A superconvergent alpha finite element method (S a FEM) for static and
free vibration analysis of shell structures
Chai et al. (2017).
"""
mutable struct FEMMShellCSDSG3{S<:AbstractFESet, F<:Function, M} <: AbstractFEMM
    integdomain::IntegDomain{S, F} # integration domain data
    material::M # material object
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

function FEMMShellCSDSG3(integdomain::IntegDomain{S, F}, material::M) where {S<:AbstractFESet, F<:Function, M}
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
    return FEMMShellCSDSG3(integdomain, material,
        _loc, _J, _J0,
        _ecoords0, _ecoords1, _edisp1, _evel1, _evel1f, _lecoords0,
        _dofnums, 
        _F0, _Ft, _FtI, _FtJ, _Te,
        _tempelmat1, _tempelmat2, _tempelmat3, _elmat, _elmatTe, _elmato, 
        _elvec, _elvecf, 
        _lloc, _lJ, _lgradN,
        _Bm, _Bb, _Bs, _DpsBmb, _DtBs, _LF, 
        _RI, _RJ, _OS)
end

function make(integdomain, material)
    return FEMMShellCSDSG3(integdomain, material)
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

function _transfmat!(Te, Ft)
    for i in 1:2*__nn
        r = (i-1)*3 .+ (1:3)
        @. Te[r, r] = Ft
    end
    return Te
end

function _triangle_centroid_area(lecoords0)
    lcentroid = (mean(view(lecoords0, :, 1)), mean(view(lecoords0, :, 2)))
    a, b = lecoords0[2, :] .- lecoords0[1, :]
    c, d = lecoords0[3, :] .- lecoords0[1, :]
    A = (a*d - b*c)/2 # area of the triangle
    return lcentroid, A
end

"""
    _Bsmat!(Bs, gradN, N)

Compute the linear transverse shear strain-displacement matrix.
"""
function _Bsmat!(Bs, lcentroid, lecoords, Ae)
    Bs .= 0.0 

    # Area of each sub triangle
    As = Ae/3

    # Subtriangle 1
    s, p, q = 3, 1, 2
    a, b = lecoords[p, :] .- lcentroid # lecoords[s, :]
    c, d = lecoords[q, :] .- lcentroid # lecoords[s, :]
    

    # The first node in the subtriangle (i.e. the centroid)
    m = (1/2/As) * As/Ae * (1/3) # multiplier
    # Node s
    co = (s - 1) * 6 # column offset
    Bs[1, co+3] += m*(b-d);                             Bs[1, co+5] += m*(As) 
    Bs[2, co+3] += m*(c-a); Bs[2, co+4] += m*(-As); 
    co = (p - 1) * 6 # column offset
    Bs[1, co+3] += m*(b-d);                             Bs[1, co+5] += m*(As) 
    Bs[2, co+3] += m*(c-a); Bs[2, co+4] += m*(-As); 
    co = (q - 1) * 6 # column offset
    Bs[1, co+3] += m*(b-d);                             Bs[1, co+5] += m*(As) 
    Bs[2, co+3] += m*(c-a); Bs[2, co+4] += m*(-As); 
    
    # The other two nodes
    m = (1/2/As) * As/Ae # multiplier
    # Node p
    co = (p - 1) * 6 # column offset
    Bs[1, co+3] += m*(d);   Bs[1, co+4] += m*(-b*d/2);  Bs[1, co+5] += m*(a*d/2) 
    Bs[2, co+3] += m*(-c);  Bs[2, co+4] += m*(b*c/2);   Bs[2, co+5] += m*(-a*c/2) 
    # Node q
    co = (q - 1) * 6 # column offset
    Bs[1, co+3] += m*(-b);  Bs[1, co+4] += m*(b*d/2);   Bs[1, co+5] += m*(-b*c/2) 
    Bs[2, co+3] += m*(a);   Bs[2, co+4] += m*(-a*d/2);  Bs[2, co+5] += m*(a*c/2) 
    
    # Subtriangle 2
    s, p, q = 1, 2, 3
    a, b = lecoords[p, :] .- lcentroid # lecoords[s, :]
    c, d = lecoords[q, :] .- lcentroid # lecoords[s, :]
    

    # The first node in the subtriangle (i.e. the centroid)
    m = (1/2/As) * As/Ae * (1/3) # multiplier
    # Node s
    co = (s - 1) * 6 # column offset
    Bs[1, co+3] += m*(b-d);                             Bs[1, co+5] += m*(As) 
    Bs[2, co+3] += m*(c-a); Bs[2, co+4] += m*(-As); 
    co = (p - 1) * 6 # column offset
    Bs[1, co+3] += m*(b-d);                             Bs[1, co+5] += m*(As) 
    Bs[2, co+3] += m*(c-a); Bs[2, co+4] += m*(-As); 
    co = (q - 1) * 6 # column offset
    Bs[1, co+3] += m*(b-d);                             Bs[1, co+5] += m*(As) 
    Bs[2, co+3] += m*(c-a); Bs[2, co+4] += m*(-As); 
    
    # The other two nodes
    m = (1/2/As) * As/Ae # multiplier
    # Node p
    co = (p - 1) * 6 # column offset
    Bs[1, co+3] += m*(d);   Bs[1, co+4] += m*(-b*d/2);  Bs[1, co+5] += m*(a*d/2) 
    Bs[2, co+3] += m*(-c);  Bs[2, co+4] += m*(b*c/2);   Bs[2, co+5] += m*(-a*c/2) 
    # Node q
    co = (q - 1) * 6 # column offset
    Bs[1, co+3] += m*(-b);  Bs[1, co+4] += m*(b*d/2);   Bs[1, co+5] += m*(-b*c/2) 
    Bs[2, co+3] += m*(a);   Bs[2, co+4] += m*(-a*d/2);  Bs[2, co+5] += m*(a*c/2) 


    # Subtriangle 3
    s, p, q = 2, 3, 1
    a, b = lecoords[p, :] .- lcentroid # lecoords[s, :]
    c, d = lecoords[q, :] .- lcentroid # lecoords[s, :]
    

    # The first node in the subtriangle (i.e. the centroid)
    m = (1/2/As) * As/Ae * (1/3) # multiplier
    # Node s
    co = (s - 1) * 6 # column offset
    Bs[1, co+3] += m*(b-d);                             Bs[1, co+5] += m*(As) 
    Bs[2, co+3] += m*(c-a); Bs[2, co+4] += m*(-As); 
    co = (p - 1) * 6 # column offset
    Bs[1, co+3] += m*(b-d);                             Bs[1, co+5] += m*(As) 
    Bs[2, co+3] += m*(c-a); Bs[2, co+4] += m*(-As); 
    co = (q - 1) * 6 # column offset
    Bs[1, co+3] += m*(b-d);                             Bs[1, co+5] += m*(As) 
    Bs[2, co+3] += m*(c-a); Bs[2, co+4] += m*(-As); 
    
    # The other two nodes
    m = (1/2/As) * As/Ae # multiplier
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

function _Bsmat1!(Bs, lcentroid, lecoords, Ae)
    Bs .= 0.0 

    # Orientation 1
    s, p, q = 3, 1, 2
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
function _Bmmat!(Bm, lcentroid, lecoords, Ae)
    # for i in 1:__nn
    #     Bm[1,6*(i-1)+1] = gradN[i,1];
    #                                   Bm[2,6*(i-1)+2] = gradN[i,2];
    #     Bm[3,6*(i-1)+1] = gradN[i,2]; Bm[3,6*(i-1)+2] = gradN[i,1];
    # end

    Bm .= 0.0 

    As = Ae/3

    # Subtriangle 1
    s, p, q = 3, 1, 2
    a, b = lecoords[p, :] .- lcentroid # lecoords[s, :]
    c, d = lecoords[q, :] .- lcentroid # lecoords[s, :]
    
    # The first node in the subtriangle (i.e. the centroid)
    m = (1/2/As) * As/Ae * (1/3) # multiplier
    # Node s
    co = (s - 1) * 6 # column offset
    Bm[1, co+1] += m*(b-d);  
                            Bm[2, co+2] += m*(c-a); 
    Bm[3, co+1] += m*(c-a); Bm[3, co+2] += m*(b-d); 
    co = (p - 1) * 6 # column offset
    Bm[1, co+1] += m*(b-d);  
                            Bm[2, co+2] += m*(c-a); 
    Bm[3, co+1] += m*(c-a); Bm[3, co+2] += m*(b-d); 
    co = (q - 1) * 6 # column offset
    Bm[1, co+1] += m*(b-d);  
                            Bm[2, co+2] += m*(c-a); 
    Bm[3, co+1] += m*(c-a); Bm[3, co+2] += m*(b-d); 
    
    # The other two nodes
    m = (1/2/As) * As/Ae # multiplier
    # Node p
    co = (p - 1) * 6 # column offset
    Bm[1, co+1] += m*(d);  
                            Bm[2, co+2] += m*(-c); 
    Bm[3, co+1] += m*(-c);  Bm[3, co+2] += m*(d); 
    # Node q
    co = (q - 1) * 6 # column offset
    Bm[1, co+1] += m*(-b);  
                            Bm[2, co+2] += m*(a); 
    Bm[3, co+1] += m*(a);  Bm[3, co+2] += m*(-b); 
        
    # Subtriangle 2
    s, p, q = 1, 2, 3
    a, b = lecoords[p, :] .- lcentroid # lecoords[s, :]
    c, d = lecoords[q, :] .- lcentroid # lecoords[s, :]
    
    # The first node in the subtriangle (i.e. the centroid)
    m = (1/2/As) * As/Ae * (1/3) # multiplier
    # Node s
    co = (s - 1) * 6 # column offset
    Bm[1, co+1] += m*(b-d);  
                            Bm[2, co+2] += m*(c-a); 
    Bm[3, co+1] += m*(c-a); Bm[3, co+2] += m*(b-d); 
    co = (p - 1) * 6 # column offset
    Bm[1, co+1] += m*(b-d);  
                            Bm[2, co+2] += m*(c-a); 
    Bm[3, co+1] += m*(c-a); Bm[3, co+2] += m*(b-d); 
    co = (q - 1) * 6 # column offset
    Bm[1, co+1] += m*(b-d);  
                            Bm[2, co+2] += m*(c-a); 
    Bm[3, co+1] += m*(c-a); Bm[3, co+2] += m*(b-d); 
    
    # The other two nodes
    m = (1/2/As) * As/Ae # multiplier
    # Node p
    co = (p - 1) * 6 # column offset
    Bm[1, co+1] += m*(d);  
                            Bm[2, co+2] += m*(-c); 
    Bm[3, co+1] += m*(-c);  Bm[3, co+2] += m*(d); 
    # Node q
    co = (q - 1) * 6 # column offset
    Bm[1, co+1] += m*(-b);  
                            Bm[2, co+2] += m*(a); 
    Bm[3, co+1] += m*(a);  Bm[3, co+2] += m*(-b); 
    
    # Subtriangle 3
    s, p, q = 2, 3, 1
    a, b = lecoords[p, :] .- lcentroid # lecoords[s, :]
    c, d = lecoords[q, :] .- lcentroid # lecoords[s, :]

    # The first node in the subtriangle (i.e. the centroid)
    m = (1/2/As) * As/Ae * (1/3) # multiplier
    # Node s
    co = (s - 1) * 6 # column offset
    Bm[1, co+1] += m*(b-d);  
                            Bm[2, co+2] += m*(c-a); 
    Bm[3, co+1] += m*(c-a); Bm[3, co+2] += m*(b-d); 
    co = (p - 1) * 6 # column offset
    Bm[1, co+1] += m*(b-d);  
                            Bm[2, co+2] += m*(c-a); 
    Bm[3, co+1] += m*(c-a); Bm[3, co+2] += m*(b-d); 
    co = (q - 1) * 6 # column offset
    Bm[1, co+1] += m*(b-d);  
                            Bm[2, co+2] += m*(c-a); 
    Bm[3, co+1] += m*(c-a); Bm[3, co+2] += m*(b-d); 
    
    # The other two nodes
    m = (1/2/As) * As/Ae # multiplier
    # Node p
    co = (p - 1) * 6 # column offset
    Bm[1, co+1] += m*(d);  
                            Bm[2, co+2] += m*(-c); 
    Bm[3, co+1] += m*(-c);  Bm[3, co+2] += m*(d); 
    # Node q
    co = (q - 1) * 6 # column offset
    Bm[1, co+1] += m*(-b);  
                            Bm[2, co+2] += m*(a); 
    Bm[3, co+1] += m*(a);  Bm[3, co+2] += m*(-b); 
        
    return Bm
end

"""
The matrices are identical:
Bm1 = similar(Bm)
_Bmmat1!(Bm1, lcentroid, lecoords0, Ae)
@infiltrate
infil> sum(abs.(Bm - Bm1))
1.1102230246251565e-16  
"""
function _Bmmat1!(Bm, lcentroid, lecoords, Ae)
    # for i in 1:__nn
    #     Bm[1,6*(i-1)+1] = gradN[i,1];
    #                                   Bm[2,6*(i-1)+2] = gradN[i,2];
    #     Bm[3,6*(i-1)+1] = gradN[i,2]; Bm[3,6*(i-1)+2] = gradN[i,1];
    # end

    Bm .= 0.0 
    
    s, p, q = 1, 2, 3
    a, b = lecoords[p, :] .- lecoords[s, :]
    c, d = lecoords[q, :] .- lecoords[s, :]
    
    # The first node in the subtriangle (i.e. the centroid)
    m = (1/2/Ae) # multiplier
    # Node s
    co = (s - 1) * 6 # column offset
    Bm[1, co+1] += m*(b-d);  
                            Bm[2, co+2] += m*(c-a); 
    Bm[3, co+1] += m*(c-a); Bm[3, co+2] += m*(b-d); 
    
    # The other two nodes
    # Node p
    co = (p - 1) * 6 # column offset
    Bm[1, co+1] += m*(d);  
                            Bm[2, co+2] += m*(-c); 
    Bm[3, co+1] += m*(-c);  Bm[3, co+2] += m*(d); 
    # Node q
    co = (q - 1) * 6 # column offset
    Bm[1, co+1] += m*(-b);  
                            Bm[2, co+2] += m*(a); 
    Bm[3, co+1] += m*(a);   Bm[3, co+2] += m*(-b); 
        
    return Bm
end

"""
    _Bbmat!(Bb, gradN)

Compute the linear, displacement independent, curvature-displacement/rotation matrix for a shell quadrilateral element with nfens=3 nodes. Displacements and rotations are in a local coordinate system.
"""
function _Bbmat!(Bb, lcentroid, lecoords, Ae)
    # for i in 1:__nn
    #                                    Bb[1,6*(i-1)+5] = gradN[i,1];
    #     Bb[2,6*(i-1)+4] = -gradN[i,2];
    #     Bb[3,6*(i-1)+4] = -gradN[i,1]; Bb[3,6*(i-1)+5] = gradN[i,2];
    # end
    Bb .= 0.0 

    As = Ae/3

    # Subtriangle 1
    s, p, q = 3, 1, 2
    a, b = lecoords[p, :] .- lcentroid # lecoords[s, :]
    c, d = lecoords[q, :] .- lcentroid # lecoords[s, :]
    
    # The first node in the subtriangle (i.e. the centroid)
    m = (1/2/As) * As/Ae * (1/3) # multiplier
    # Node s
    co = (s - 1) * 6 # column offset
                             Bb[1, co+5] += m*(b-d);  
    Bb[2, co+4] += -m*(c-a); 
    Bb[3, co+4] += -m*(b-d); Bb[3, co+5] += m*(c-a); 
    co = (p - 1) * 6 # column offset
                             Bb[1, co+5] += m*(b-d);  
    Bb[2, co+4] += -m*(c-a); 
    Bb[3, co+4] += -m*(b-d); Bb[3, co+5] += m*(c-a); 
    co = (q - 1) * 6 # column offset
                             Bb[1, co+5] += m*(b-d);  
    Bb[2, co+4] += -m*(c-a); 
    Bb[3, co+4] += -m*(b-d); Bb[3, co+5] += m*(c-a); 
    
    # The other two nodes
    m = (1/2/As) * As/Ae # multiplier
    # Node p
    co = (p - 1) * 6 # column offset
                             Bb[1, co+5] += m*(d);  
    Bb[2, co+4] += -m*(-c); 
    Bb[3, co+4] += -m*(d); Bb[3, co+5] += m*(-c); 
    # Node q
    co = (q - 1) * 6 # column offset
                             Bb[1, co+5] += m*(-b);  
    Bb[2, co+4] += -m*(a); 
    Bb[3, co+4] += -m*(-b); Bb[3, co+5] += m*(a); 
        
    # Subtriangle 2
    s, p, q = 1, 2, 3
    a, b = lecoords[p, :] .- lcentroid # lecoords[s, :]
    c, d = lecoords[q, :] .- lcentroid # lecoords[s, :]
    
    # The first node in the subtriangle (i.e. the centroid)
    m = (1/2/As) * As/Ae * (1/3) # multiplier
    # Node s
    co = (s - 1) * 6 # column offset
                             Bb[1, co+5] += m*(b-d);  
    Bb[2, co+4] += -m*(c-a); 
    Bb[3, co+4] += -m*(b-d); Bb[3, co+5] += m*(c-a); 
    co = (p - 1) * 6 # column offset
                             Bb[1, co+5] += m*(b-d);  
    Bb[2, co+4] += -m*(c-a); 
    Bb[3, co+4] += -m*(b-d); Bb[3, co+5] += m*(c-a); 
    co = (q - 1) * 6 # column offset
                             Bb[1, co+5] += m*(b-d);  
    Bb[2, co+4] += -m*(c-a); 
    Bb[3, co+4] += -m*(b-d); Bb[3, co+5] += m*(c-a); 
    
    # The other two nodes
    m = (1/2/As) * As/Ae # multiplier
    # Node p
    co = (p - 1) * 6 # column offset
                             Bb[1, co+5] += m*(d);  
    Bb[2, co+4] += -m*(-c); 
    Bb[3, co+4] += -m*(d); Bb[3, co+5] += m*(-c); 
    # Node q
    co = (q - 1) * 6 # column offset
                             Bb[1, co+5] += m*(-b);  
    Bb[2, co+4] += -m*(a); 
    Bb[3, co+4] += -m*(-b); Bb[3, co+5] += m*(a); 
    
    # Subtriangle 3
    s, p, q = 2, 3, 1
    a, b = lecoords[p, :] .- lcentroid # lecoords[s, :]
    c, d = lecoords[q, :] .- lcentroid # lecoords[s, :]

    # The first node in the subtriangle (i.e. the centroid)
    m = (1/2/As) * As/Ae * (1/3) # multiplier
    # Node s
    co = (s - 1) * 6 # column offset
                             Bb[1, co+5] += m*(b-d);  
    Bb[2, co+4] += -m*(c-a); 
    Bb[3, co+4] += -m*(b-d); Bb[3, co+5] += m*(c-a); 
    co = (p - 1) * 6 # column offset
                             Bb[1, co+5] += m*(b-d);  
    Bb[2, co+4] += -m*(c-a); 
    Bb[3, co+4] += -m*(b-d); Bb[3, co+5] += m*(c-a); 
    co = (q - 1) * 6 # column offset
                             Bb[1, co+5] += m*(b-d);  
    Bb[2, co+4] += -m*(c-a); 
    Bb[3, co+4] += -m*(b-d); Bb[3, co+5] += m*(c-a); 
    
    # The other two nodes
    m = (1/2/As) * As/Ae # multiplier
    # Node p
    co = (p - 1) * 6 # column offset
                             Bb[1, co+5] += m*(d);  
    Bb[2, co+4] += -m*(-c); 
    Bb[3, co+4] += -m*(d); Bb[3, co+5] += m*(-c); 
    # Node q
    co = (q - 1) * 6 # column offset
                             Bb[1, co+5] += m*(-b);  
    Bb[2, co+4] += -m*(a); 
    Bb[3, co+4] += -m*(-b); Bb[3, co+5] += m*(a); 
    
    return Bb
end

function _Bbmat1!(Bb, lcentroid, lecoords, Ae)
    # for i in 1:__nn
    #                                    Bb[1,6*(i-1)+5] = gradN[i,1];
    #     Bb[2,6*(i-1)+4] = -gradN[i,2];
    #     Bb[3,6*(i-1)+4] = -gradN[i,1]; Bb[3,6*(i-1)+5] = gradN[i,2];
    # end
    Bb .= 0.0 

    # Subtriangle 1
    s, p, q = 1, 2, 3
    a, b = lecoords[p, :] .- lecoords[s, :]
    c, d = lecoords[q, :] .- lecoords[s, :]
    
    # The first node in the subtriangle (i.e. the centroid)
    m = (1/2/Ae)   # multiplier
    # Node s
    co = (s - 1) * 6 # column offset
                             Bb[1, co+5] += m*(b-d);  
    Bb[2, co+4] += -m*(c-a); 
    Bb[3, co+4] += -m*(b-d); Bb[3, co+5] += m*(c-a); 
    
    # The other two nodes
    # Node p
    co = (p - 1) * 6 # column offset
                             Bb[1, co+5] += m*(d);  
    Bb[2, co+4] += -m*(-c); 
    Bb[3, co+4] += -m*(d); Bb[3, co+5] += m*(-c); 
    # Node q
    co = (q - 1) * 6 # column offset
                             Bb[1, co+5] += m*(-b);  
    Bb[2, co+4] += -m*(a); 
    Bb[3, co+4] += -m*(-b); Bb[3, co+5] += m*(a); 
        
    return Bb
end

"""
    stiffness(self::FEMMShellCSDSG3, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the material stiffness matrix.
"""
function stiffness(self::FEMMShellCSDSG3, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {ASS<:AbstractSysmatAssembler, T<:Number, TI<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    @assert npts == 1 "Must use one-point integration rule. "
    loc, J, J0 = self._loc, self._J, self._J0
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    lecoords0 = self._lecoords0
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
    elmat, elmatTe = self._elmat, self._elmatTe
    Bm, Bb, Bs, DpsBmb, DtBs = self._Bm, self._Bb, self._Bs, self._DpsBmb, self._DtBs
    Dps, Dt = _shell_material_stiffness(self.material)
    scf=5/6;  # shear correction factor
    Dt .*= scf
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), dchi.nfreedofs, dchi.nfreedofs);
    for i in 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        _compute_J0!(J0, ecoords0)
        local_frame!(delegateof(fes), Ft, J0)
        fill!(elmat,  0.0); # Initialize element matrix
        locjac!(loc, J, ecoords0, Ns[1], gradNparams[1])
        t = self.integdomain.otherdimension(loc, fes.conn[i], Ns[1])
        mul!(lecoords0, ecoords0, view(Ft, :, 1:2))
        lcentroid, Ae = _triangle_centroid_area(lecoords0)
        _Bmmat1!(Bm, lcentroid, lecoords0, Ae)
        _Bbmat1!(Bb, lcentroid, lecoords0, Ae)
        add_btdb_ut_only!(elmat, Bm, t*Ae, Dps, DpsBmb)
        add_btdb_ut_only!(elmat, Bb, (t^3)/12*Ae, Dps, DpsBmb)
        # Transverse shear
        _Bsmat1!(Bs, lcentroid, lecoords0, Ae)
        # The stabilization expression has a huge effect (at least for the
        # pinched cylinder). What is the recommended multiplier of he^2?
        he = sqrt(Ae)
        add_btdb_ut_only!(elmat, Bs, (t^3/(t^2+0.2*he^2))*Ae, Dt, DtBs)
        # add_btdb_ut_only!(elmat, Bs, t*Ae, Dt, DtBs)
    
        # Apply drilling-rotation artificial stiffness
        kavg4 = mean((elmat[4, 4], elmat[10, 10], elmat[16, 16]))
        kavg5 = mean((elmat[5, 5], elmat[11, 11], elmat[17, 17]))
        kavg = (kavg4 + kavg5) / 1e5
        elmat[6, 6] += kavg 
        elmat[12, 12] += kavg 
        elmat[18, 18] += kavg 
        complete_lt!(elmat)
        # d = [-0.2887, -0.2887, 0.5774, -0.2887, -0.2887, 0.5774]
        # d = fill(0.0, 18)
        # @infiltrate
        # Transformation into global ordinates
        _transfmat!(Te, Ft)
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        # Assembly
        gatherdofnums!(dchi, dofnums, fes.conn[i]); 
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function stiffness(self::FEMMShellCSDSG3, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{TI}) where {T<:Number, TI<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi);
end


"""
    mass(self::FEMMShellCSDSG3,  assembler::A,  geom::NodalField{FFlt}, dchi::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the consistent mass matrix

This is a general routine for the shell FEMM.
"""
function mass(self::FEMMShellCSDSG3,  assembler::A,  geom0::NodalField{FFlt}, dchi::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    loc, J, J0 = self._loc, self._J, self._J0
    ecoords0, ecoords1, edisp1, dofnums = self._ecoords0, self._ecoords1, self._edisp1, self._dofnums
    lecoords0 = self._lecoords0
    F0, Ft, FtI, FtJ, Te = self._F0, self._Ft, self._FtI, self._FtJ, self._Te
    R1I, R1J = self._RI, self._RJ
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

function mass(self::FEMMShellCSDSG3,  geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return mass(self, assembler, geom, u);
end

end # module