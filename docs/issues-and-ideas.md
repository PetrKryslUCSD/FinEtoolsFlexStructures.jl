
Issues and ideas:

-- Documenter:
using FinEtoolsFlexStructures
import Pkg; Pkg.add("DocumenterTools")
using DocumenterTools
Travis.genkeys(user="PetrKryslUCSD", repo="https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl")
Pkg.rm("DocumenterTools")

- Transformation matrix for an element should be calculated just once at the beginning of the integration for each element. The transformation then needs to be applied when the integration loop is finished.
- Use averaging for the transverse sheer strain-displacement matrix.



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
#     transformwith = QTEQTransformer(elmat)
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