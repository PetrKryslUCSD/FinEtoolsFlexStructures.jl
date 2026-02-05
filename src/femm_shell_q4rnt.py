from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal

import numpy as np
import scipy.sparse as sp

from finetools.csys import CSys
from finetools.fields import NodalField, gathervalues_asmat, gathervalues_asvec
from finetools.integ_domain import IntegDomain
from finetools.integ_rules import QuadRule

from .feset_shell_q4rnt import FESetShellQ4RNT
from .femm_shell_t3ff import (
    NodalTriadsE,
    TRANSV_SHEAR_FORMULATION_AVERAGE_B,
    TRANSV_SHEAR_FORMULATION_AVERAGE_K,
    _add_btdb,
    _compute_nodal_normal,
    _complete_lt,
    _e_g,
    _isoparametric,
    _shell_material_stiffness,
)
from .transformer import Layup2ElementAngle

NN = 4  # number of nodes per element
NDOF = 6  # degrees of freedom per node


def _shape_q4(N: np.ndarray, dN_dxi: np.ndarray, dN_deta: np.ndarray, xi: float, eta: float) -> None:
    """
    Bilinear Q4 shape functions on [-1,1]^2 with node ordering:
      1(-1,-1), 2(1,-1), 3(1,1), 4(-1,1).
    """
    xi = float(xi)
    eta = float(eta)
    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta)
    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta)
    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta)
    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta)

    dN_dxi[0] = -0.25 * (1.0 - eta)
    dN_dxi[1] = 0.25 * (1.0 - eta)
    dN_dxi[2] = 0.25 * (1.0 + eta)
    dN_dxi[3] = -0.25 * (1.0 + eta)

    dN_deta[0] = -0.25 * (1.0 - xi)
    dN_deta[1] = -0.25 * (1.0 + xi)
    dN_deta[2] = 0.25 * (1.0 + xi)
    dN_deta[3] = 0.25 * (1.0 - xi)


def _compute_tangents(J0: np.ndarray, ecoords: np.ndarray, dN_dxi: np.ndarray, dN_deta: np.ndarray) -> np.ndarray:
    ecoords = np.asarray(ecoords, dtype=float).reshape(NN, 3)
    dN_dxi = np.asarray(dN_dxi, dtype=float).reshape(NN)
    dN_deta = np.asarray(dN_deta, dtype=float).reshape(NN)
    J0[:, 0] = ecoords.T @ dN_dxi
    J0[:, 1] = ecoords.T @ dN_deta
    return J0


def _map_to_global(XYZ: np.ndarray, ecoords: np.ndarray, N: np.ndarray) -> np.ndarray:
    ecoords = np.asarray(ecoords, dtype=float).reshape(NN, 3)
    N = np.asarray(N, dtype=float).reshape(NN)
    XYZ[:] = N @ ecoords
    return XYZ


def _gradN_e(gradN_e: np.ndarray, J0: np.ndarray, E_G: np.ndarray, dN_dxi: np.ndarray, dN_deta: np.ndarray) -> np.ndarray:
    """
    Gradient of shape functions in the local orthonormal basis (e1,e2) at an integration point.

    Steps:
    - Compute the physical gradient in R^3 on the tangent plane:
        ∇N = J0 * (J0^T J0)^(-1) * [dN/dxi, dN/deta]^T
    - Project onto the orthonormal tangent basis (e1,e2).
    """
    G = J0.T @ J0
    detG = float(G[0, 0] * G[1, 1] - G[0, 1] * G[1, 0])
    if detG == 0.0:
        gradN_e[:, :] = 0.0
        return gradN_e
    invG = np.linalg.inv(G)
    E2 = E_G[:, 0:2]
    tmp = np.zeros(2, dtype=float)
    for a in range(NN):
        tmp[0] = float(dN_dxi[a])
        tmp[1] = float(dN_deta[a])
        g = J0 @ (invG @ tmp)
        gradN_e[a, 0] = float(np.dot(E2[:, 0], g))
        gradN_e[a, 1] = float(np.dot(E2[:, 1], g))
    return gradN_e


def _Bmmat(Bm: np.ndarray, gradN: np.ndarray) -> None:
    Bm[:, :] = 0.0
    for i in range(NN):
        Bm[0, NDOF * i + 0] = gradN[i, 0]
        Bm[1, NDOF * i + 1] = gradN[i, 1]
        Bm[2, NDOF * i + 0] = gradN[i, 1]
        Bm[2, NDOF * i + 1] = gradN[i, 0]


def _Bbmat(Bb: np.ndarray, gradN: np.ndarray) -> None:
    Bb[:, :] = 0.0
    for i in range(NN):
        Bb[0, NDOF * i + 4] = gradN[i, 0]
        Bb[1, NDOF * i + 3] = -gradN[i, 1]
        Bb[2, NDOF * i + 3] = -gradN[i, 0]
        Bb[2, NDOF * i + 4] = gradN[i, 1]


def _Bsmat(Bs: np.ndarray, gradN: np.ndarray, N: np.ndarray) -> None:
    """
    Reissner--Mindlin transverse shear strains:
      γ13 = w,1 + θ2
      γ23 = w,2 - θ1
    in the element local basis (e1,e2,e3).
    """
    Bs[:, :] = 0.0
    for i in range(NN):
        Bs[0, NDOF * i + 2] = gradN[i, 0]
        Bs[0, NDOF * i + 4] = N[i]
        Bs[1, NDOF * i + 2] = gradN[i, 1]
        Bs[1, NDOF * i + 3] = -N[i]


class TransfmatGToA:
    def __init__(self) -> None:
        self.Tblock = np.zeros((3, 3), dtype=float)

    def __call__(self, T: np.ndarray, A_Es: list[np.ndarray], E_G: np.ndarray) -> np.ndarray:
        T[:, :] = 0.0
        for i in range(NN):
            self.Tblock[:, :] = A_Es[i].T @ E_G.T
            offset = i * NDOF
            r = slice(offset + 0, offset + 3)
            T[r, r] = self.Tblock
            r = slice(offset + 3, offset + 6)
            T[r, r] = self.Tblock
        return T


def _transfmat_a_to_e(T: np.ndarray, A_Es: list[np.ndarray], gradN_e: np.ndarray) -> np.ndarray:
    """
    Transformation from nodal 'a' frames to element 'e' frame with drilling decoupling.

    Direct transcription of the T3FF transformation generalized to Q4.
    """
    T[:, :] = 0.0
    for i in range(NN):
        roffst = i * NDOF
        iA_E = A_Es[i]
        iA_E_33 = float(iA_E[2, 2])
        # translation block
        T[roffst : roffst + 3, roffst : roffst + 3] = iA_E
        # rotational 2x2 block with drilling decoupling
        for cl in range(2):
            for rw in range(2):
                T[roffst + 3 + rw, roffst + 3 + cl] = iA_E[rw, cl] - (
                    (1.0 / iA_E_33) * iA_E[rw, 2] * iA_E[cl, 2]
                )
        m1 = (1.0 / iA_E_33) * iA_E[0, 2]
        m2 = (1.0 / iA_E_33) * iA_E[1, 2]
        for j in range(NN):
            coffst = j * NDOF
            for k in range(3):
                a3 = 0.5 * (iA_E[1, k] * gradN_e[j, 0] - iA_E[0, k] * gradN_e[j, 1])
                T[roffst + 3, coffst + k] += m1 * a3
                T[roffst + 4, coffst + k] += m2 * a3
    return T


def _area_from_rule(ecoords: np.ndarray, rule: QuadRule, *, N: np.ndarray, dN_dxi: np.ndarray, dN_deta: np.ndarray) -> float:
    J0 = np.zeros((3, 2), dtype=float)
    area = 0.0
    for q, (xi, eta) in enumerate(rule.points):
        _shape_q4(N, dN_dxi, dN_deta, float(xi), float(eta))
        _compute_tangents(J0, ecoords, dN_dxi, dN_deta)
        g1 = J0[:, 0]
        g2 = J0[:, 1]
        detA = float(np.linalg.norm(np.cross(g1, g2)))
        area += float(rule.weights[q]) * detA
    return float(area)


def _drilling_penalty_kavg(
    elmat: np.ndarray,
    normals: np.ndarray,
    normal_valid: np.ndarray,
    conn: np.ndarray,
    drilling_stiffness_scale: float,
) -> float:
    """
    Compute a robust drilling-penalty stiffness scale.

    The scale should follow the tangential bending stiffness regardless of the
    orientation of the shell normal with respect to the global axes. We estimate
    this from the rotational 3x3 sub-matrices projected onto the tangent plane.
    """
    drilling_stiffness_scale = float(drilling_stiffness_scale)
    if drilling_stiffness_scale == 0.0:
        return 0.0

    I = np.eye(3, dtype=float)
    tangential = []
    for k in range(NN):
        nnode = int(conn[k])
        if not bool(normal_valid[nnode]):
            continue
        nvec = np.asarray(normals[nnode, :], dtype=float).reshape(3)
        nn = float(np.linalg.norm(nvec))
        if nn == 0.0:
            continue
        nvec /= nn
        r = slice(k * NDOF + 3, k * NDOF + 6)
        Krr = np.asarray(elmat[r, r], dtype=float)
        P = I - np.outer(nvec, nvec)
        Kt = P @ Krr @ P
        tangential.append(max(0.0, float(np.trace(Kt)) / 2.0))

    if not tangential:
        return 0.0
    return float(np.mean(tangential) * drilling_stiffness_scale)


@dataclass(slots=True)
class FEMMShellQ4RNT:
    """
    4-node quadrilateral shell element (Q4RNT).

    - Membrane and bending: Gauss 2x2 integration.
    - Transverse shear: reduced integration (1 point at center).
    - Drilling rotation: stabilized with a penalty proportional to the average
      flexural rotational stiffness, applied about the nodal normal.
    - "RNT" refers to the robust nodal triads used for drilling regularization.
    """

    integdomain: IntegDomain
    mcsys: CSys
    material: object

    transv_shear_formulation: int = TRANSV_SHEAR_FORMULATION_AVERAGE_B
    drilling_stiffness_scale: float = 1.0
    threshold_angle: float = 30.0
    mult_el_size: float = 5.0 / 12.0 / 1.5

    _associatedgeometry: bool = False
    _normals: np.ndarray | None = None
    _normal_valid: np.ndarray | None = None

    # Buffers
    _centroid: np.ndarray | None = None
    _J0: np.ndarray | None = None
    _ecoords: np.ndarray | None = None
    _edisp: np.ndarray | None = None
    _E_G: np.ndarray | None = None
    _A_Es: list[np.ndarray] | None = None
    _nvalid: np.ndarray | None = None
    _Tga: np.ndarray | None = None
    _Tae: np.ndarray | None = None
    _T: np.ndarray | None = None
    _gradN_e: np.ndarray | None = None
    _N: np.ndarray | None = None
    _dN_dxi: np.ndarray | None = None
    _dN_deta: np.ndarray | None = None
    _Bm: np.ndarray | None = None
    _Bb: np.ndarray | None = None
    _Bs: np.ndarray | None = None
    _Bhat3: np.ndarray | None = None
    _Bhat2: np.ndarray | None = None
    _DpsB: np.ndarray | None = None
    _DtB: np.ndarray | None = None

    # Scratch
    _XYZ: np.ndarray | None = None

    def __post_init__(self) -> None:
        fes = self.integdomain.fes
        if not isinstance(getattr(fes, "delegate", None), FESetShellQ4RNT):
            raise ValueError("Expected integdomain.fes.delegate to be FESetShellQ4RNT")
        if fes.nodes_per_elem() != NN:
            raise ValueError("Expected quadrilateral connectivity (nel,4)")

        nnmax = int(np.max(fes.conn)) + 1 if fes.count() > 0 else 0
        self._normals = np.zeros((nnmax, 3), dtype=float)
        self._normal_valid = np.ones(nnmax, dtype=bool)

        self._centroid = np.zeros((1, 3), dtype=float)
        self._J0 = np.zeros((3, 2), dtype=float)
        self._ecoords = np.zeros((NN, 3), dtype=float)
        self._edisp = np.zeros(NN * NDOF, dtype=float)
        self._E_G = np.zeros((3, 3), dtype=float)
        self._A_Es = [np.zeros((3, 3), dtype=float) for _ in range(NN)]
        self._nvalid = np.zeros(NN, dtype=bool)
        self._Tga = np.zeros((NN * NDOF, NN * NDOF), dtype=float)
        self._Tae = np.zeros((NN * NDOF, NN * NDOF), dtype=float)
        self._T = np.zeros((NN * NDOF, NN * NDOF), dtype=float)
        self._gradN_e = np.zeros((NN, 2), dtype=float)
        self._N = np.zeros(NN, dtype=float)
        self._dN_dxi = np.zeros(NN, dtype=float)
        self._dN_deta = np.zeros(NN, dtype=float)
        self._Bm = np.zeros((3, NN * NDOF), dtype=float)
        self._Bb = np.zeros((3, NN * NDOF), dtype=float)
        self._Bs = np.zeros((2, NN * NDOF), dtype=float)
        self._Bhat3 = np.zeros((3, NN * NDOF), dtype=float)
        self._Bhat2 = np.zeros((2, NN * NDOF), dtype=float)
        self._DpsB = np.zeros_like(self._Bhat3)
        self._DtB = np.zeros_like(self._Bhat2)
        self._XYZ = np.zeros(3, dtype=float)

    def associategeometry(self, geom: NodalField) -> "FEMMShellQ4RNT":
        threshold_angle = float(self.threshold_angle)
        J0 = self._J0
        normals = self._normals
        normal_valid = self._normal_valid
        ecoords = self._ecoords
        E_G = self._E_G
        N = self._N
        dN_dxi = self._dN_dxi
        dN_deta = self._dN_deta
        assert (
            J0 is not None
            and normals is not None
            and normal_valid is not None
            and ecoords is not None
            and E_G is not None
            and N is not None
            and dN_dxi is not None
            and dN_deta is not None
        )

        normals[:, :] = 0.0
        normal_valid[:] = True

        fes = self.integdomain.fes
        corner_ip = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]], dtype=float)

        nnormal = np.zeros(3, dtype=float)
        for e in range(fes.count()):
            conn = np.asarray(fes.conn[e, :], dtype=int)
            ecoords[:, :] = gathervalues_asmat(geom, conn)
            for loc in range(NN):
                xi, eta = (float(corner_ip[loc, 0]), float(corner_ip[loc, 1]))
                _shape_q4(N, dN_dxi, dN_deta, xi, eta)
                _compute_tangents(J0, ecoords, dN_dxi, dN_deta)
                # Area-weighted contribution: use unnormalized normal from tangents.
                _compute_nodal_normal(nnormal, self.mcsys, geom.values[int(conn[loc]), :], J0, e, loc)
                # Scale by local Jacobian to approximate area weighting.
                detA = float(np.linalg.norm(np.cross(J0[:, 0], J0[:, 1])))
                normals[int(conn[loc]), :] += detA * nnormal

        nn = np.linalg.norm(normals, axis=1)
        mask = nn > 0.0
        normals[mask, :] /= nn[mask][:, None]

        ntolerance = 1.0 - float(np.sqrt(1.0 - np.sin(threshold_angle / 180.0 * np.pi) ** 2))
        for e in range(fes.count()):
            conn = np.asarray(fes.conn[e, :], dtype=int)
            ecoords[:, :] = gathervalues_asmat(geom, conn)
            for loc in range(NN):
                xi, eta = (float(corner_ip[loc, 0]), float(corner_ip[loc, 1]))
                _shape_q4(N, dN_dxi, dN_deta, xi, eta)
                _compute_tangents(J0, ecoords, dN_dxi, dN_deta)
                _e_g(E_G, J0)
                n = int(conn[loc])
                nd = float(np.dot(normals[n, :], E_G[:, 2]))
                if nd < 1.0 - ntolerance:
                    # Invalidate nodal normals at sharp creases to avoid spurious drilling coupling.
                    normal_valid[n] = False

        self._associatedgeometry = True
        return self

    def num_normals(self) -> tuple[int, int]:
        assert self._normal_valid is not None
        total = int(self._normal_valid.shape[0])
        invalid = int(np.sum(~self._normal_valid))
        return total, invalid

    def stiffness(
        self,
        assembler: object,
        geom0: NodalField,
        u1: NodalField,
        Rfield1: NodalField,
        dchi: NodalField,
    ) -> sp.csr_matrix:
        if not self._associatedgeometry:
            raise RuntimeError("Call associategeometry() first")

        fes = self.integdomain.fes
        normals = self._normals
        normal_valid = self._normal_valid
        centroid = self._centroid
        J0 = self._J0
        ecoords = self._ecoords
        E_G = self._E_G
        A_Es = self._A_Es
        nvalid = self._nvalid
        Tga = self._Tga
        Tae = self._Tae
        T = self._T
        gradN_e = self._gradN_e
        N = self._N
        dN_dxi = self._dN_dxi
        dN_deta = self._dN_deta
        Bm = self._Bm
        Bb = self._Bb
        Bs = self._Bs
        Bhat3 = self._Bhat3
        Bhat2 = self._Bhat2
        DpsB = self._DpsB
        DtB = self._DtB
        assert (
            normals is not None
            and normal_valid is not None
            and centroid is not None
            and J0 is not None
            and ecoords is not None
            and E_G is not None
            and A_Es is not None
            and nvalid is not None
            and Tga is not None
            and Tae is not None
            and T is not None
            and gradN_e is not None
            and N is not None
            and dN_dxi is not None
            and dN_deta is not None
            and Bm is not None
            and Bb is not None
            and Bs is not None
            and Bhat3 is not None
            and Bhat2 is not None
            and DpsB is not None
            and DtB is not None
        )

        nodal_triads_e = NodalTriadsE()
        transfmat_g_to_a = TransfmatGToA()

        Dps, Dt = _shell_material_stiffness(self.material)
        Dt = Dt.copy()
        Dt *= 5.0 / 6.0

        rule_mb = self.integdomain.rule
        if not hasattr(rule_mb, "points") or not hasattr(rule_mb, "weights"):
            raise ValueError("Expected a quadrilateral rule with points/weights in integdomain.rule")

        # Design choice for Q4RNT: full integration for membrane/bending, reduced for shear.
        shear_rule = QuadRule(1)
        mult_el_size = float(self.mult_el_size)
        transv_shear_formulation = int(self.transv_shear_formulation)
        drilling_stiffness_scale = float(self.drilling_stiffness_scale)

        ndofs = dchi.nfreedofs()
        rows: list[int] = []
        cols: list[int] = []
        data: list[float] = []

        elmat = np.zeros((NN * NDOF, NN * NDOF), dtype=float)

        for e in range(fes.count()):
            conn = np.asarray(fes.conn[e, :], dtype=int)
            ecoords[:, :] = gathervalues_asmat(geom0, conn)
            centroid[0, :] = np.mean(ecoords, axis=0)

            area_total = _area_from_rule(ecoords, QuadRule(2), N=N, dN_dxi=dN_dxi, dN_deta=dN_deta)

            elmat[:, :] = 0.0

            # Membrane + bending (Gauss 2x2).
            for q, (xi, eta) in enumerate(np.asarray(rule_mb.points, dtype=float).reshape(-1, 2)):
                wq = float(np.asarray(rule_mb.weights, dtype=float).reshape(-1)[q])
                _shape_q4(N, dN_dxi, dN_deta, float(xi), float(eta))
                _compute_tangents(J0, ecoords, dN_dxi, dN_deta)
                detA = float(np.linalg.norm(np.cross(J0[:, 0], J0[:, 1])))
                if detA == 0.0:
                    continue
                _e_g(E_G, J0)
                _gradN_e(gradN_e, J0, E_G, dN_dxi, dN_deta)

                t = float(self.integdomain.otherdimension(centroid, conn, np.array([xi, eta], dtype=float)))

                nodal_triads_e(A_Es, nvalid, E_G, normals, normal_valid, conn)
                transfmat_g_to_a(Tga, A_Es, E_G)
                _transfmat_a_to_e(Tae, A_Es, gradN_e)
                T[:, :] = Tae @ Tga

                _Bmmat(Bm, gradN_e)
                Bhat3[:, :] = Bm @ T
                _add_btdb(elmat, Bhat3, (t * detA) * wq, Dps, DpsB)

                _Bbmat(Bb, gradN_e)
                Bhat3[:, :] = Bb @ T
                _add_btdb(elmat, Bhat3, ((t**3) / 12.0 * detA) * wq, Dps, DpsB)

            # Shear (reduced integration at center).
            xi, eta = 0.0, 0.0
            _shape_q4(N, dN_dxi, dN_deta, xi, eta)
            _compute_tangents(J0, ecoords, dN_dxi, dN_deta)
            detA = float(np.linalg.norm(np.cross(J0[:, 0], J0[:, 1])))
            if detA != 0.0:
                _e_g(E_G, J0)
                _gradN_e(gradN_e, J0, E_G, dN_dxi, dN_deta)
                t = float(self.integdomain.otherdimension(centroid, conn, np.array([xi, eta], dtype=float)))

                nodal_triads_e(A_Es, nvalid, E_G, normals, normal_valid, conn)
                transfmat_g_to_a(Tga, A_Es, E_G)
                _transfmat_a_to_e(Tae, A_Es, gradN_e)
                T[:, :] = Tae @ Tga

                shear_fac = (t**3) / (t**2 + mult_el_size * 2.0 * area_total) if area_total > 0.0 else t

                if transv_shear_formulation == TRANSV_SHEAR_FORMULATION_AVERAGE_K:
                    # Use the same energy-averaging idea as T3FF by sampling shear
                    # at the 4 corners and averaging stiffness contributions.
                    corner_ip = np.array(
                        [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]], dtype=float
                    )
                    for (xk, ek) in corner_ip:
                        _shape_q4(N, dN_dxi, dN_deta, float(xk), float(ek))
                        _compute_tangents(J0, ecoords, dN_dxi, dN_deta)
                        detAk = float(np.linalg.norm(np.cross(J0[:, 0], J0[:, 1])))
                        if detAk == 0.0:
                            continue
                        _e_g(E_G, J0)
                        _gradN_e(gradN_e, J0, E_G, dN_dxi, dN_deta)
                        nodal_triads_e(A_Es, nvalid, E_G, normals, normal_valid, conn)
                        transfmat_g_to_a(Tga, A_Es, E_G)
                        _transfmat_a_to_e(Tae, A_Es, gradN_e)
                        T[:, :] = Tae @ Tga
                        _Bsmat(Bs, gradN_e, N)
                        Bhat2[:, :] = Bs @ T
                        _add_btdb(elmat, Bhat2, (shear_fac * detAk) * (float(shear_rule.weights[0]) / 4.0), Dt, DtB)
                else:
                    _Bsmat(Bs, gradN_e, N)
                    Bhat2[:, :] = Bs @ T
                    _add_btdb(elmat, Bhat2, (shear_fac * detA) * float(shear_rule.weights[0]), Dt, DtB)

            # Drilling stabilization (in global rotations about the nodal normal).
            if drilling_stiffness_scale != 0.0:
                kavg = _drilling_penalty_kavg(elmat, normals, normal_valid, conn, drilling_stiffness_scale)
                if kavg != 0.0:
                    for k in range(NN):
                        nnode = int(conn[k])
                        if not bool(normal_valid[nnode]):
                            continue
                        nvec = np.asarray(normals[nnode, :], dtype=float).reshape(3)
                        if float(np.linalg.norm(nvec)) == 0.0:
                            continue
                        # Coordinate-free drilling penalty aligned with the nodal normal.
                        Ke = kavg * np.outer(nvec, nvec)
                        r = slice(k * NDOF + 3, k * NDOF + 6)
                        elmat[r, r] += Ke

            _complete_lt(elmat)

            gdofs = dchi.gatherdofnums(conn)
            free = np.nonzero(gdofs >= 0)[0]
            if free.size == 0:
                continue
            g = gdofs[free]
            ke = elmat[np.ix_(free, free)]
            rr, cc = np.meshgrid(g, g, indexing="ij")
            rows.extend(rr.reshape(-1).tolist())
            cols.extend(cc.reshape(-1).tolist())
            data.extend(ke.reshape(-1).tolist())

        K = sp.coo_matrix((data, (rows, cols)), shape=(ndofs, ndofs)).tocsr()
        return K

    def mass(self, assembler: object, geom0: NodalField, dchi: NodalField) -> sp.csr_matrix:
        if not self._associatedgeometry:
            raise RuntimeError("Call associategeometry() first")

        fes = self.integdomain.fes
        centroid = self._centroid
        J0 = self._J0
        ecoords = self._ecoords
        N = self._N
        dN_dxi = self._dN_dxi
        dN_deta = self._dN_deta
        assert centroid is not None and J0 is not None and ecoords is not None
        assert N is not None and dN_dxi is not None and dN_deta is not None

        rho = float(self.material.massdensity())
        rule = QuadRule(2)
        ndofs = dchi.nfreedofs()
        rows: list[int] = []
        cols: list[int] = []
        data: list[float] = []

        for e in range(fes.count()):
            conn = np.asarray(fes.conn[e, :], dtype=int)
            ecoords[:, :] = gathervalues_asmat(geom0, conn)
            centroid[0, :] = np.mean(ecoords, axis=0)

            area_total = _area_from_rule(ecoords, rule, N=N, dN_dxi=dN_dxi, dN_deta=dN_deta)
            if area_total == 0.0:
                continue
            t = float(self.integdomain.otherdimension(centroid, conn, np.array([0.0, 0.0], dtype=float)))

            tmass = rho * (t * area_total) / float(NN)
            rmass = rho * ((t**3) / 12.0 * area_total) / float(NN)

            gdofs = dchi.gatherdofnums(conn)
            for k in range(NN):
                for d in range(3):
                    g = int(gdofs[k * NDOF + d])
                    if g >= 0:
                        rows.append(g)
                        cols.append(g)
                        data.append(tmass)
                for d in range(3, 6):
                    g = int(gdofs[k * NDOF + d])
                    if g >= 0:
                        rows.append(g)
                        cols.append(g)
                        data.append(rmass)

        M = sp.coo_matrix((data, (rows, cols)), shape=(ndofs, ndofs)).tocsr()
        return M

    def inspectintegpoints(
        self,
        geom0: NodalField,
        u: NodalField,
        felist: np.ndarray,
        inspector: Callable,
        idat,
        quantity: Literal["moment", "membrane", "shear", "bending_moment", "transverse_shear", "membrane_force"] = "moment",
        *,
        outputcsys: CSys | None = None,
    ):
        if not self._associatedgeometry:
            raise RuntimeError("Call associategeometry() first")

        fes = self.integdomain.fes
        normals = self._normals
        normal_valid = self._normal_valid
        centroid = self._centroid
        J0 = self._J0
        ecoords = self._ecoords
        edisp = self._edisp
        E_G = self._E_G
        A_Es = self._A_Es
        nvalid = self._nvalid
        Tga = self._Tga
        Tae = self._Tae
        gradN_e = self._gradN_e
        N = self._N
        dN_dxi = self._dN_dxi
        dN_deta = self._dN_deta
        Bm = self._Bm
        Bb = self._Bb
        Bs = self._Bs
        XYZ = self._XYZ
        assert (
            normals is not None
            and normal_valid is not None
            and centroid is not None
            and J0 is not None
            and ecoords is not None
            and edisp is not None
            and E_G is not None
            and A_Es is not None
            and nvalid is not None
            and Tga is not None
            and Tae is not None
            and gradN_e is not None
            and N is not None
            and dN_dxi is not None
            and dN_deta is not None
            and Bm is not None
            and Bb is not None
            and Bs is not None
            and XYZ is not None
        )

        nodal_triads_e = NodalTriadsE()
        transfmat_g_to_a = TransfmatGToA()
        lla = Layup2ElementAngle()

        Dps, Dt = _shell_material_stiffness(self.material)
        Dt = Dt.copy()
        Dt *= 5.0 / 6.0

        if outputcsys is None:
            outputcsys = self.mcsys

        quant = "moment"
        if quantity in ("bending", "moment", "bending_moment"):
            quant = "moment"
            rule = QuadRule(2)
        elif quantity in ("membrane", "membrane_force"):
            quant = "membrane"
            rule = QuadRule(2)
        elif quantity in ("shear", "transverse", "transverse_shear"):
            quant = "shear"
            rule = QuadRule(1)
        else:
            raise ValueError(f"Unsupported quantity '{quantity}'")

        mult_el_size = float(self.mult_el_size)

        edisp_a = np.zeros_like(edisp)
        edisp_e = np.zeros_like(edisp)
        out = np.zeros(3, dtype=float)
        o2_e = np.zeros((2, 2), dtype=float)

        for ilist in range(len(felist)):
            e = int(felist[ilist])
            conn = np.asarray(fes.conn[e, :], dtype=int)
            ecoords[:, :] = gathervalues_asmat(geom0, conn)
            centroid[0, :] = np.mean(ecoords, axis=0)
            edisp[:] = gathervalues_asvec(u, conn)

            area_total = _area_from_rule(ecoords, QuadRule(2), N=N, dN_dxi=dN_dxi, dN_deta=dN_deta)

            for q, (xi, eta) in enumerate(np.asarray(rule.points, dtype=float).reshape(-1, 2)):
                _shape_q4(N, dN_dxi, dN_deta, float(xi), float(eta))
                _map_to_global(XYZ, ecoords, N)
                _compute_tangents(J0, ecoords, dN_dxi, dN_deta)
                _e_g(E_G, J0)
                _gradN_e(gradN_e, J0, E_G, dN_dxi, dN_deta)
                t = float(self.integdomain.otherdimension(centroid, conn, np.array([xi, eta], dtype=float)))

                nodal_triads_e(A_Es, nvalid, E_G, normals, normal_valid, conn)
                transfmat_g_to_a(Tga, A_Es, E_G)
                edisp_a[:] = Tga @ edisp
                _transfmat_a_to_e(Tae, A_Es, gradN_e)
                edisp_e[:] = Tae @ edisp_a

                outputcsys.updatecsmat(XYZ.reshape(1, 3), J0, e, q)
                ocsm, ocsn = lla(E_G, outputcsys.csmat)
                o2_e[0, 0] = o2_e[1, 1] = ocsm
                o2_e[0, 1] = ocsn
                o2_e[1, 0] = -ocsn

                if quant == "moment":
                    _Bbmat(Bb, gradN_e)
                    kurv = Bb @ edisp_e
                    mom = (t**3 / 12.0) * (Dps @ kurv)
                    m = np.array([[mom[0], mom[2]], [mom[2], mom[1]]], dtype=float)
                    mo = o2_e.T @ m @ o2_e
                    out[:] = [mo[0, 0], mo[1, 1], mo[0, 1]]
                elif quant == "membrane":
                    _Bmmat(Bm, gradN_e)
                    strn = Bm @ edisp_e
                    frc = t * (Dps @ strn)
                    f = np.array([[frc[0], frc[2]], [frc[2], frc[1]]], dtype=float)
                    fo = o2_e.T @ f @ o2_e
                    out[:] = [fo[0, 0], fo[1, 1], fo[0, 1]]
                elif quant == "shear":
                    _Bsmat(Bs, gradN_e, N)
                    shr = Bs @ edisp_e
                    shear_fac = (t**3) / (t**2 + mult_el_size * 2.0 * area_total) if area_total > 0.0 else t
                    frc = shear_fac * (Dt @ shr)
                    fo = o2_e.T @ frc
                    out[0:2] = [fo[0], fo[1]]
                    out[2] = 0.0
                else:  # pragma: no cover
                    raise ValueError(f"Unsupported quantity '{quantity}'")

                idat = inspector(idat, e, conn, ecoords, out.copy(), XYZ.copy())

        return idat


def make_q4rnt(integdomain: IntegDomain, material, mcsys: CSys | None = None) -> FEMMShellQ4RNT:
    if mcsys is None:
        mcsys = CSys(3, 3, _isoparametric)
    return FEMMShellQ4RNT(integdomain=integdomain, mcsys=mcsys, material=material)
