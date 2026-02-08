"""
    FEMMShellQ4RNTModule

4-node quadrilateral shell element (Q4RNT) implementation for Reissner-Mindlin plate theory.

- Membrane and bending: Gauss 2×2 full integration
- Transverse shear: 1-point reduced integration at center
- Drilling rotation: stabilized with penalty proportional to average flexural rotational stiffness
- RNT: Robust nodal triads for drilling regularization
"""
module FEMMShellQ4RNTModule

using LinearAlgebra
using SparseArrays
using Statistics

const NN = 4  # nodes per element
const NDOF = 6  # DOF per node (3 translation + 3 rotation)

# Transverse shear formulation flags
const TRANSV_SHEAR_FORMULATION_AVERAGE_B = 1
const TRANSV_SHEAR_FORMULATION_AVERAGE_K = 2

export FEMMShellQ4RNT, make_q4rnt
export TRANSV_SHEAR_FORMULATION_AVERAGE_B, TRANSV_SHEAR_FORMULATION_AVERAGE_K

"""
    _shape_q4!(N, dN_dxi, dN_deta, xi, eta)

Bilinear Q4 shape functions on [-1,1]² reference element.
Node ordering: 1(-1,-1), 2(1,-1), 3(1,1), 4(-1,1)
"""
function _shape_q4!(N::Vector, dN_dxi::Vector, dN_deta::Vector, xi::Real, eta::Real)
    xi = Float64(xi)
    eta = Float64(eta)
    
    N[1] = 0.25 * (1.0 - xi) * (1.0 - eta)
    N[2] = 0.25 * (1.0 + xi) * (1.0 - eta)
    N[3] = 0.25 * (1.0 + xi) * (1.0 + eta)
    N[4] = 0.25 * (1.0 - xi) * (1.0 + eta)
    
    dN_dxi[1] = -0.25 * (1.0 - eta)
    dN_dxi[2] = 0.25 * (1.0 - eta)
    dN_dxi[3] = 0.25 * (1.0 + eta)
    dN_dxi[4] = -0.25 * (1.0 + eta)
    
    dN_deta[1] = -0.25 * (1.0 - xi)
    dN_deta[2] = -0.25 * (1.0 + xi)
    dN_deta[3] = 0.25 * (1.0 + xi)
    dN_deta[4] = 0.25 * (1.0 - xi)
    
    return nothing
end

"""
    _compute_tangents!(J0, ecoords, dN_dxi, dN_deta)

Compute tangent vectors to surface at an integration point.
J0 is 3×2: columns are dX/dξ and dX/dη tangent vectors.
"""
function _compute_tangents!(J0::Matrix, ecoords::Matrix, dN_dxi::Vector, dN_deta::Vector)
    J0[:, 1] = ecoords' * dN_dxi
    J0[:, 2] = ecoords' * dN_deta
    return nothing
end

"""
    _map_to_global!(XYZ, ecoords, N)

Map reference element point to global coordinates.
"""
function _map_to_global!(XYZ::Vector, ecoords::Matrix, N::Vector)
    XYZ[:] = ecoords' * N
    return nothing
end

"""
    _gradN_e!(gradN_e, J0, E_G, dN_dxi, dN_deta)

Gradient of shape functions in the local orthonormal element basis (e1, e2).

Computes: ∇N = J0 * (J0^T J0)^(-1) * [dN/dξ, dN/dη]^T
Then projects onto the orthonormal tangent basis (e1, e2).
"""
function _gradN_e!(gradN_e::Matrix, J0::Matrix, E_G::Matrix, dN_dxi::Vector, dN_deta::Vector)
    G = J0' * J0
    detG = G[1, 1] * G[2, 2] - G[1, 2] * G[2, 1]
    
    if abs(detG) < 1.0e-14
        fill!(gradN_e, 0.0)
        return nothing
    end
    
    invG = inv(G)
    E2 = @view E_G[:, 1:2]
    
    for a = 1:NN
        tmp = [dN_dxi[a], dN_deta[a]]
        g = J0 * (invG * tmp)
        gradN_e[a, 1] = dot(E2[:, 1], g)
        gradN_e[a, 2] = dot(E2[:, 2], g)
    end
    
    return nothing
end

"""
    _e_g!(E_G, J0)

Compute orthonormal basis (e1, e2, e3) from tangent vectors.
E_G is 3×3 where columns are the orthonormal basis vectors.
"""
function _e_g!(E_G::Matrix, J0::Matrix)
    g1 = @view J0[:, 1]
    g2 = @view J0[:, 2]
    
    # Local e3: normal to the surface
    e3 = cross(g1, g2)
    len_e3 = norm(e3)
    if len_e3 > 0.0
        e3 = e3 / len_e3
    end
    
    # Local e1: normalize g1
    e1 = g1 / norm(g1)
    
    # Local e2: orthogonal to both e1 and e3
    e2 = cross(e3, e1)
    
    E_G[:, 1] = e1
    E_G[:, 2] = e2
    E_G[:, 3] = e3
    
    return nothing
end

"""
    _Bmmat!(Bm, gradN)

Membrane strain-displacement matrix.
Relates 2D membrane strain to nodal DOF: [ε11, ε22, 2ε12]^T = Bm * u
"""
function _Bmmat!(Bm::Matrix, gradN::Matrix)
    fill!(Bm, 0.0)
    for i = 1:NN
        Bm[1, NDOF * (i - 1) + 1] = gradN[i, 1]
        Bm[2, NDOF * (i - 1) + 2] = gradN[i, 2]
        Bm[3, NDOF * (i - 1) + 1] = gradN[i, 2]
        Bm[3, NDOF * (i - 1) + 2] = gradN[i, 1]
    end
    return nothing
end

"""
    _Bbmat!(Bb, gradN)

Bending strain-displacement matrix.
Relates curvature to nodal DOF: [κ11, κ22, 2κ12]^T = Bb * u
"""
function _Bbmat!(Bb::Matrix, gradN::Matrix)
    fill!(Bb, 0.0)
    for i = 1:NN
        Bb[1, NDOF * (i - 1) + 5] = gradN[i, 1]
        Bb[2, NDOF * (i - 1) + 4] = -gradN[i, 2]
        Bb[3, NDOF * (i - 1) + 4] = -gradN[i, 1]
        Bb[3, NDOF * (i - 1) + 5] = gradN[i, 2]
    end
    return nothing
end

"""
    _Bsmat!(Bs, gradN, N)

Transverse shear strain-displacement matrix (Reissner-Mindlin).
Relates transverse shear strains to nodal DOF:
  γ13 = w,1 + θ2
  γ23 = w,2 - θ1
"""
function _Bsmat!(Bs::Matrix, gradN::Matrix, N::Vector)
    fill!(Bs, 0.0)
    for i = 1:NN
        Bs[1, NDOF * (i - 1) + 3] = gradN[i, 1]
        Bs[1, NDOF * (i - 1) + 5] = N[i]
        Bs[2, NDOF * (i - 1) + 3] = gradN[i, 2]
        Bs[2, NDOF * (i - 1) + 4] = -N[i]
    end
    return nothing
end

"""
    TransfmatGToA

Transformation from global 'g' to nodal 'a' frames.
"""
mutable struct TransfmatGToA
    Tblock::Matrix{Float64}
    
    function TransfmatGToA()
        return new(zeros(3, 3))
    end
end

function (tfm::TransfmatGToA)(T::Matrix, A_Es::Vector, E_G::Matrix)
    fill!(T, 0.0)
    for i = 1:NN
        tfm.Tblock[:, :] = A_Es[i]' * E_G'
        offset = (i - 1) * NDOF
        r1 = (offset + 1):(offset + 3)
        r2 = (offset + 4):(offset + 6)
        T[r1, r1] = tfm.Tblock
        T[r2, r2] = tfm.Tblock
    end
    return T
end

"""
    _transfmat_a_to_e!(T, A_Es, gradN_e)

Transformation matrix from nodal 'a' frames to element 'e' frame with drilling decoupling.
Direct adaptation of T3FF transformation to Q4.
"""
function _transfmat_a_to_e!(T::Matrix, A_Es::Vector, gradN_e::Matrix)
    fill!(T, 0.0)
    
    for i = 1:NN
        roffst = (i - 1) * NDOF
        iA_E = A_Es[i]
        iA_E_33 = iA_E[3, 3]
        
        # Translation block: 3×3
        T[(roffst + 1):(roffst + 3), (roffst + 1):(roffst + 3)] = iA_E
        
        # Rotational 2×2 block with drilling decoupling
        for cl = 1:2, rw = 1:2
            T[roffst + 3 + rw, roffst + 3 + cl] = iA_E[rw, cl] - (iA_E[rw, 3] * iA_E[cl, 3] / iA_E_33)
        end
        
        m1 = iA_E[1, 3] / iA_E_33
        m2 = iA_E[2, 3] / iA_E_33
        
        for j = 1:NN
            coffst = (j - 1) * NDOF
            for k = 1:3
                a3 = 0.5 * (iA_E[2, k] * gradN_e[j, 1] - iA_E[1, k] * gradN_e[j, 2])
                T[roffst + 4, coffst + k] += m1 * a3
                T[roffst + 5, coffst + k] += m2 * a3
            end
        end
    end
    
    return T
end

"""
    _add_btdb!(elmat, B, weight, D, scratch)

Add contribution B^T D B to element matrix with given weight.
scratch is a work vector.
"""
function _add_btdb!(elmat::Matrix, B::Matrix, weight::Real, D::Matrix, scratch::Matrix)
    mul!(scratch, D, B)
    mul!(elmat, B', scratch, weight, 1.0)
    return nothing
end

"""
    _complete_lt!(A)

Complete a lower-triangular symmetric matrix by copying to upper triangle.
"""
function _complete_lt!(A::Matrix)
    n = size(A, 1)
    for i = 1:n, j = (i + 1):n
        A[i, j] = A[j, i]
    end
    return nothing
end

"""
    _area_from_rule(ecoords, rule_points, rule_weights, N, dN_dxi, dN_deta)

Compute element area by numerical integration over the reference element.
"""
function _area_from_rule(ecoords::Matrix, rule_points::Matrix, rule_weights::Vector,
                         N::Vector, dN_dxi::Vector, dN_deta::Vector)
    J0 = zeros(3, 2)
    area = 0.0
    
    for q = 1:size(rule_points, 1)
        xi, eta = rule_points[q, 1], rule_points[q, 2]
        _shape_q4!(N, dN_dxi, dN_deta, xi, eta)
        _compute_tangents!(J0, ecoords, dN_dxi, dN_deta)
        
        g1 = @view J0[:, 1]
        g2 = @view J0[:, 2]
        detA = norm(cross(g1, g2))
        area += rule_weights[q] * detA
    end
    
    return area
end

"""
    _drilling_penalty_kavg(elmat, normals, normal_valid, conn, drilling_stiffness_scale)

Compute drilling penalty stiffness scale based on average tangential bending stiffness.
"""
function _drilling_penalty_kavg(elmat::Matrix, normals::Matrix, normal_valid::Vector,
                                conn::Vector, drilling_stiffness_scale::Real)
    drilling_stiffness_scale = Float64(drilling_stiffness_scale)
    if drilling_stiffness_scale == 0.0
        return 0.0
    end
    
    I3 = Matrix{Float64}(I, 3, 3)
    tangential = Float64[]
    
    for k = 1:NN
        nnode = Int(conn[k])
        if !normal_valid[nnode]
            continue
        end
        
        nvec = @view normals[nnode, :]
        nn = norm(nvec)
        if nn == 0.0
            continue
        end
        nvec_normalized = nvec / nn
        
        r = ((k - 1) * NDOF + 4):((k - 1) * NDOF + 6)
        Krr = @view elmat[r, r]
        P = I3 - nvec_normalized * nvec_normalized'
        Kt = P * Krr * P
        
        push!(tangential, max(0.0, tr(Kt) / 2.0))
    end
    
    if isempty(tangential)
        return 0.0
    end
    
    return mean(tangential) * drilling_stiffness_scale
end

"""
    FEMMShellQ4RNT

4-node quadrilateral shell finite element with Reissner-Mindlin plate theory.

Attributes:
- `integdomain`: Integration domain with quadrature rule
- `mcsys`: Material coordinate system
- `material`: Material with stiffness and density properties
- `transv_shear_formulation`: Shear formulation flag (AVERAGE_B or AVERAGE_K)
- `drilling_stiffness_scale`: Drilling penalty parameter (0 to disable)
- `threshold_angle`: Angle threshold for crease detection (degrees)
- `mult_el_size`: Multiplier for shear penalty coupling
"""
mutable struct FEMMShellQ4RNT
    integdomain::Any
    mcsys::Any
    material::Any
    
    transv_shear_formulation::Int = TRANSV_SHEAR_FORMULATION_AVERAGE_B
    drilling_stiffness_scale::Float64 = 1.0
    threshold_angle::Float64 = 30.0
    mult_el_size::Float64 = 5.0 / 12.0 / 1.5
    
    _associatedgeometry::Bool = false
    _normals::Matrix{Float64} = Matrix{Float64}(undef, 0, 3)
    _normal_valid::Vector{Bool} = Vector{Bool}(undef, 0)
    
    # Geometry buffers
    _centroid::Matrix{Float64} = zeros(1, 3)
    _J0::Matrix{Float64} = zeros(3, 2)
    _ecoords::Matrix{Float64} = zeros(NN, 3)
    _edisp::Vector{Float64} = zeros(NN * NDOF)
    _E_G::Matrix{Float64} = zeros(3, 3)
    _A_Es::Vector{Matrix{Float64}} = [zeros(3, 3) for _ = 1:NN]
    _nvalid::Vector{Bool} = zeros(Bool, NN)
    _Tga::Matrix{Float64} = zeros(NN * NDOF, NN * NDOF)
    _Tae::Matrix{Float64} = zeros(NN * NDOF, NN * NDOF)
    _T::Matrix{Float64} = zeros(NN * NDOF, NN * NDOF)
    _gradN_e::Matrix{Float64} = zeros(NN, 2)
    _N::Vector{Float64} = zeros(NN)
    _dN_dxi::Vector{Float64} = zeros(NN)
    _dN_deta::Vector{Float64} = zeros(NN)
    _Bm::Matrix{Float64} = zeros(3, NN * NDOF)
    _Bb::Matrix{Float64} = zeros(3, NN * NDOF)
    _Bs::Matrix{Float64} = zeros(2, NN * NDOF)
    _Bhat3::Matrix{Float64} = zeros(3, NN * NDOF)
    _Bhat2::Matrix{Float64} = zeros(2, NN * NDOF)
    _DpsB::Matrix{Float64} = zeros(3, NN * NDOF)
    _DtB::Matrix{Float64} = zeros(2, NN * NDOF)
    _XYZ::Vector{Float64} = zeros(3)
end

function FEMMShellQ4RNT(integdomain, mcsys, material; kw...)
    obj = FEMMShellQ4RNT(integdomain, mcsys, material)
    
    fes = integdomain.fes
    nnmax = length(fes) > 0 ? maximum(fes.conn[:]) : 0
    obj._normals = zeros(nnmax, 3)
    obj._normal_valid = ones(Bool, nnmax)
    
    for (k, v) in pairs(kw)
        setproperty!(obj, k, v)
    end
    
    return obj
end

"""
    associategeometry!(femm, geom)

Associate geometry field and compute nodal normals with crease detection.
"""
function associategeometry!(femm::FEMMShellQ4RNT, geom)
    threshold_angle = femm.threshold_angle
    fes = femm.integdomain.fes
    
    # First pass: accumulate area-weighted normals
    fill!(femm._normals, 0.0)
    fill!(femm._normal_valid, true)
    
    corner_ip = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]]
    
    for e = 1:length(fes)
        conn = fes.conn[e, :]
        for i = 1:NN
            femm._ecoords[i, :] = geom.values[conn[i], :]
        end
        
        for (loc, (xi, eta)) in enumerate(corner_ip)
            _shape_q4!(femm._N, femm._dN_dxi, femm._dN_deta, xi, eta)
            _compute_tangents!(femm._J0, femm._ecoords, femm._dN_dxi, femm._dN_deta)
            
            # Nodal normal (unnormalized)
            g1 = @view femm._J0[:, 1]
            g2 = @view femm._J0[:, 2]
            nnormal = cross(g1, g2)
            detA = norm(nnormal)
            
            nnode = conn[loc]
            femm._normals[nnode, :] .+= detA .* nnormal
        end
    end
    
    # Normalize nodal normals
    for i = 1:size(femm._normals, 1)
        nn = norm(@view femm._normals[i, :])
        if nn > 0.0
            femm._normals[i, :] ./= nn
        end
    end
    
    # Second pass: detect creases
    ntolerance = 1.0 - sqrt(1.0 - sin(threshold_angle / 180.0 * π)^2)
    
    for e = 1:length(fes)
        conn = fes.conn[e, :]
        for i = 1:NN
            femm._ecoords[i, :] = geom.values[conn[i], :]
        end
        
        for (loc, (xi, eta)) in enumerate(corner_ip)
            _shape_q4!(femm._N, femm._dN_dxi, femm._dN_deta, xi, eta)
            _compute_tangents!(femm._J0, femm._ecoords, femm._dN_dxi, femm._dN_deta)
            _e_g!(femm._E_G, femm._J0)
            
            nnode = conn[loc]
            nd = dot(femm._normals[nnode, :], femm._E_G[:, 3])
            
            if nd < 1.0 - ntolerance
                femm._normal_valid[nnode] = false
            end
        end
    end
    
    femm._associatedgeometry = true
    return femm
end

"""
    num_normals(femm)

Return (total_normals, invalid_normals) count.
"""
function num_normals(femm::FEMMShellQ4RNT)
    total = size(femm._normals, 1)
    invalid = count(.!femm._normal_valid)
    return (total, invalid)
end

"""
    stiffness(femm, assembler, geom0, u1, Rfield1, dchi)

Compute element stiffness matrix and assemble global matrix.
"""
function stiffness(femm::FEMMShellQ4RNT, assembler, geom0, u1, Rfield1, dchi)
    if !femm._associatedgeometry
        error("Call associategeometry!() first")
    end
    
    fes = femm.integdomain.fes
    normals = femm._normals
    normal_valid = femm._normal_valid
    centroid = femm._centroid
    J0 = femm._J0
    ecoords = femm._ecoords
    E_G = femm._E_G
    A_Es = femm._A_Es
    nvalid = femm._nvalid
    Tga = femm._Tga
    Tae = femm._Tae
    T = femm._T
    gradN_e = femm._gradN_e
    N = femm._N
    dN_dxi = femm._dN_dxi
    dN_deta = femm._dN_deta
    Bm = femm._Bm
    Bb = femm._Bb
    Bs = femm._Bs
    Bhat3 = femm._Bhat3
    Bhat2 = femm._Bhat2
    DpsB = femm._DpsB
    DtB = femm._DtB
    
    # Get material stiffness matrices (placeholder - should be implemented)
    Dps = zeros(3, 3)  # Membrane/bending stiffness
    Dt = zeros(2, 2)   # Shear stiffness
    Dt .*= 5.0 / 6.0
    
    # Quadrature rules
    rule_mb_points = [-1/sqrt(3), 1/sqrt(3)]
    rule_mb = [(xi, eta) for xi in rule_mb_points for eta in rule_mb_points]
    rule_mb_weights = ones(length(rule_mb))
    
    # 1-point rule for shear
    shear_rule_weights = [4.0]
    
    mult_el_size = Float64(femm.mult_el_size)
    transv_shear_formulation = Int(femm.transv_shear_formulation)
    drilling_stiffness_scale = Float64(femm.drilling_stiffness_scale)
    
    ndofs = maximum(dchi) > 0 ? maximum(dchi) : 0
    rows = Int[]
    cols = Int[]
    data = Float64[]
    
    elmat = zeros(NN * NDOF, NN * NDOF)
    
    nodal_triads_e = nothing  # Placeholder - requires NodalTriadsE implementation
    transfmat_g_to_a = TransfmatGToA()
    
    for e = 1:length(fes)
        conn = fes.conn[e, :]
        for i = 1:NN
            ecoords[i, :] = geom0.values[conn[i], :]
        end
        centroid[1, :] = vec(mean(ecoords; dims=1))
        
        # Calculate total area
        rule_points_2x2 = [(xi, eta) for xi in rule_mb_points for eta in rule_mb_points]
        rule_weights_2x2 = ones(length(rule_points_2x2))
        area_total = _area_from_rule(ecoords, hcat(rule_points_2x2...)', rule_weights_2x2, N, dN_dxi, dN_deta)
        
        fill!(elmat, 0.0)
        
        # Membrane + bending (2×2 Gauss)
        for (q, (xi, eta)) in enumerate(rule_points_2x2)
            wq = 1.0
            _shape_q4!(N, dN_dxi, dN_deta, Float64(xi), Float64(eta))
            _compute_tangents!(J0, ecoords, dN_dxi, dN_deta)
            g1 = @view J0[:, 1]
            g2 = @view J0[:, 2]
            detA = norm(cross(g1, g2))
            if detA == 0.0
                continue
            end
            _e_g!(E_G, J0)
            _gradN_e!(gradN_e, J0, E_G, dN_dxi, dN_deta)
            
            t = 1.0  # Placeholder for thickness
            
            # Compute nodal triads (placeholder)
            for i = 1:NN
                A_Es[i] = Matrix{Float64}(I, 3, 3)
                nvalid[i] = true
            end
            
            transfmat_g_to_a(Tga, A_Es, E_G)
            _transfmat_a_to_e!(Tae, A_Es, gradN_e)
            mul!(T, Tae, Tga)
            
            _Bmmat!(Bm, gradN_e)
            mul!(Bhat3, Bm, T)
            _add_btdb!(elmat, Bhat3, (t * detA) * wq, Dps, DpsB)
            
            _Bbmat!(Bb, gradN_e)
            mul!(Bhat3, Bb, T)
            _add_btdb!(elmat, Bhat3, ((t^3) / 12.0 * detA) * wq, Dps, DpsB)
        end
        
        # Shear (1-point reduced integration at center)
        xi, eta = 0.0, 0.0
        _shape_q4!(N, dN_dxi, dN_deta, xi, eta)
        _compute_tangents!(J0, ecoords, dN_dxi, dN_deta)
        g1 = @view J0[:, 1]
        g2 = @view J0[:, 2]
        detA = norm(cross(g1, g2))
        if detA != 0.0
            _e_g!(E_G, J0)
            _gradN_e!(gradN_e, J0, E_G, dN_dxi, dN_deta)
            t = 1.0  # Placeholder for thickness
            
            for i = 1:NN
                A_Es[i] = Matrix{Float64}(I, 3, 3)
                nvalid[i] = true
            end
            
            transfmat_g_to_a(Tga, A_Es, E_G)
            _transfmat_a_to_e!(Tae, A_Es, gradN_e)
            mul!(T, Tae, Tga)
            
            shear_fac = area_total > 0.0 ? (t^3) / (t^2 + mult_el_size * 2.0 * area_total) : t
            
            if transv_shear_formulation == TRANSV_SHEAR_FORMULATION_AVERAGE_K
                corner_ip = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]]
                for (xk, ek) in corner_ip
                    _shape_q4!(N, dN_dxi, dN_deta, Float64(xk), Float64(ek))
                    _compute_tangents!(J0, ecoords, dN_dxi, dN_deta)
                    g1_k = @view J0[:, 1]
                    g2_k = @view J0[:, 2]
                    detAk = norm(cross(g1_k, g2_k))
                    if detAk == 0.0
                        continue
                    end
                    _e_g!(E_G, J0)
                    _gradN_e!(gradN_e, J0, E_G, dN_dxi, dN_deta)
                    
                    for i = 1:NN
                        A_Es[i] = Matrix{Float64}(I, 3, 3)
                        nvalid[i] = true
                    end
                    
                    transfmat_g_to_a(Tga, A_Es, E_G)
                    _transfmat_a_to_e!(Tae, A_Es, gradN_e)
                    mul!(T, Tae, Tga)
                    
                    _Bsmat!(Bs, gradN_e, N)
                    mul!(Bhat2, Bs, T)
                    _add_btdb!(elmat, Bhat2, (shear_fac * detAk) * (shear_rule_weights[1] / 4.0), Dt, DtB)
                end
            else
                _Bsmat!(Bs, gradN_e, N)
                mul!(Bhat2, Bs, T)
                _add_btdb!(elmat, Bhat2, (shear_fac * detA) * shear_rule_weights[1], Dt, DtB)
            end
        end
        
        # Drilling stabilization
        if drilling_stiffness_scale != 0.0
            kavg = _drilling_penalty_kavg(elmat, normals, normal_valid, conn, drilling_stiffness_scale)
            if kavg != 0.0
                for k = 1:NN
                    nnode = Int(conn[k])
                    if !normal_valid[nnode]
                        continue
                    end
                    nvec = @view normals[nnode, :]
                    if norm(nvec) == 0.0
                        continue
                    end
                    Ke = kavg .* (nvec * nvec')
                    r = ((k - 1) * NDOF + 4):((k - 1) * NDOF + 6)
                    elmat[r, r] .+= Ke
                end
            end
        end
        
        _complete_lt!(elmat)
        
        # Extract free DOFs and assemble
        gdofs = vec(dchi[conn, :])
        free_mask = gdofs .>= 0
        if any(free_mask)
            free_idx = findall(free_mask)
            g = gdofs[free_idx]
            ke = elmat[free_idx, free_idx]
            
            for i in eachindex(free_idx), j in eachindex(free_idx)
                push!(rows, g[i])
                push!(cols, g[j])
                push!(data, ke[i, j])
            end
        end
    end
    
    K = sparse(rows, cols, data, ndofs, ndofs)
    return K
end

"""
    mass(femm, assembler, geom0, dchi)

Compute element mass matrix and assemble global matrix.
"""
function mass(femm::FEMMShellQ4RNT, assembler, geom0, dchi)
    if !femm._associatedgeometry
        error("Call associategeometry!() first")
    end
    
    fes = femm.integdomain.fes
    centroid = femm._centroid
    J0 = femm._J0
    ecoords = femm._ecoords
    N = femm._N
    dN_dxi = femm._dN_dxi
    dN_deta = femm._dN_deta
    
    rho = 1.0  # Placeholder for material density
    
    # 2×2 Gauss rule
    rule_points_2x2 = [(-1/sqrt(3), -1/sqrt(3)), (1/sqrt(3), -1/sqrt(3)), 
                       (1/sqrt(3), 1/sqrt(3)), (-1/sqrt(3), 1/sqrt(3))]
    rule_weights_2x2 = ones(length(rule_points_2x2))
    
    ndofs = maximum(dchi) > 0 ? maximum(dchi) : 0
    rows = Int[]
    cols = Int[]
    data = Float64[]
    
    for e = 1:length(fes)
        conn = fes.conn[e, :]
        for i = 1:NN
            ecoords[i, :] = geom0.values[conn[i], :]
        end
        centroid[1, :] = vec(mean(ecoords; dims=1))
        
        area_total = _area_from_rule(ecoords, hcat(rule_points_2x2...)', rule_weights_2x2, N, dN_dxi, dN_deta)
        if area_total == 0.0
            continue
        end
        
        t = 1.0  # Placeholder for thickness
        
        tmass = rho * (t * area_total) / Float64(NN)
        rmass = rho * ((t^3) / 12.0 * area_total) / Float64(NN)
        
        gdofs = vec(dchi[conn, :])
        
        for k = 1:NN
            for d = 1:3
                g = Int(gdofs[(k - 1) * NDOF + d])
                if g >= 0
                    push!(rows, g)
                    push!(cols, g)
                    push!(data, tmass)
                end
            end
            for d = 4:6
                g = Int(gdofs[(k - 1) * NDOF + d])
                if g >= 0
                    push!(rows, g)
                    push!(cols, g)
                    push!(data, rmass)
                end
            end
        end
    end
    
    M = sparse(rows, cols, data, ndofs, ndofs)
    return M
end

"""
    inspectintegpoints(femm, geom0, u, felist, inspector, idat; quantity="moment", outputcsys=nothing)

Inspect integration point quantities (moment, membrane force, transverse shear).
"""
function inspectintegpoints(femm::FEMMShellQ4RNT, geom0, u, felist, inspector, idat;
                            quantity::String="moment", outputcsys=nothing)
    if !femm._associatedgeometry
        error("Call associategeometry!() first")
    end
    
    fes = femm.integdomain.fes
    normals = femm._normals
    normal_valid = femm._normal_valid
    centroid = femm._centroid
    J0 = femm._J0
    ecoords = femm._ecoords
    edisp = femm._edisp
    E_G = femm._E_G
    A_Es = femm._A_Es
    nvalid = femm._nvalid
    Tga = femm._Tga
    Tae = femm._Tae
    gradN_e = femm._gradN_e
    N = femm._N
    dN_dxi = femm._dN_dxi
    dN_deta = femm._dN_deta
    Bm = femm._Bm
    Bb = femm._Bb
    Bs = femm._Bs
    XYZ = femm._XYZ
    
    # Get material stiffness (placeholder)
    Dps = zeros(3, 3)
    Dt = zeros(2, 2)
    Dt .*= 5.0 / 6.0
    
    if outputcsys === nothing
        outputcsys = femm.mcsys
    end
    
    quant = "moment"
    if quantity in ("bending", "moment", "bending_moment")
        quant = "moment"
        rule_points = [(-1/sqrt(3), -1/sqrt(3)), (1/sqrt(3), -1/sqrt(3)), 
                       (1/sqrt(3), 1/sqrt(3)), (-1/sqrt(3), 1/sqrt(3))]
    elseif quantity in ("membrane", "membrane_force")
        quant = "membrane"
        rule_points = [(-1/sqrt(3), -1/sqrt(3)), (1/sqrt(3), -1/sqrt(3)), 
                       (1/sqrt(3), 1/sqrt(3)), (-1/sqrt(3), 1/sqrt(3))]
    elseif quantity in ("shear", "transverse", "transverse_shear")
        quant = "shear"
        rule_points = [(0.0, 0.0)]
    else
        error("Unsupported quantity '$quantity'")
    end
    
    mult_el_size = Float64(femm.mult_el_size)
    
    edisp_a = zeros(NN * NDOF)
    edisp_e = zeros(NN * NDOF)
    out = zeros(3)
    o2_e = zeros(2, 2)
    
    rule_points_2x2 = [(-1/sqrt(3), -1/sqrt(3)), (1/sqrt(3), -1/sqrt(3)), 
                       (1/sqrt(3), 1/sqrt(3)), (-1/sqrt(3), 1/sqrt(3))]
    rule_weights_2x2 = ones(length(rule_points_2x2))
    
    for ilist in eachindex(felist)
        e = Int(felist[ilist])
        conn = fes.conn[e, :]
        for i = 1:NN
            ecoords[i, :] = geom0.values[conn[i], :]
            edisp[(i-1)*NDOF+1:i*NDOF] = vec(u.values[conn[i], :])
        end
        centroid[1, :] = vec(mean(ecoords; dims=1))
        
        area_total = _area_from_rule(ecoords, hcat(rule_points_2x2...)', rule_weights_2x2, N, dN_dxi, dN_deta)
        
        for (q, (xi, eta)) in enumerate(rule_points)
            _shape_q4!(N, dN_dxi, dN_deta, Float64(xi), Float64(eta))
            _map_to_global!(XYZ, ecoords, N)
            _compute_tangents!(J0, ecoords, dN_dxi, dN_deta)
            _e_g!(E_G, J0)
            _gradN_e!(gradN_e, J0, E_G, dN_dxi, dN_deta)
            
            t = 1.0  # Placeholder for thickness
            
            for i = 1:NN
                A_Es[i] = Matrix{Float64}(I, 3, 3)
                nvalid[i] = true
            end
            
            transfmat_g_to_a = TransfmatGToA()
            transfmat_g_to_a(Tga, A_Es, E_G)
            mul!(edisp_a, Tga, edisp)
            
            _transfmat_a_to_e!(Tae, A_Es, gradN_e)
            mul!(edisp_e, Tae, edisp_a)
            
            # Output coordinate system transform (placeholder)
            ocsm = 1.0
            ocsn = 0.0
            o2_e[1, 1] = ocsm
            o2_e[2, 2] = ocsm
            o2_e[1, 2] = ocsn
            o2_e[2, 1] = -ocsn
            
            if quant == "moment"
                _Bbmat!(Bb, gradN_e)
                kurv = Bb * edisp_e
                mom = ((t^3 / 12.0) .* (Dps * kurv))
                m = [mom[1] mom[3]; mom[3] mom[2]]
                mo = o2_e' * m * o2_e
                out[:] = [mo[1, 1], mo[2, 2], mo[1, 2]]
            elseif quant == "membrane"
                _Bmmat!(Bm, gradN_e)
                strn = Bm * edisp_e
                frc = t .* (Dps * strn)
                f = [frc[1] frc[3]; frc[3] frc[2]]
                fo = o2_e' * f * o2_e
                out[:] = [fo[1, 1], fo[2, 2], fo[1, 2]]
            elseif quant == "shear"
                _Bsmat!(Bs, gradN_e, N)
                shr = Bs * edisp_e
                shear_fac = area_total > 0.0 ? (t^3) / (t^2 + mult_el_size * 2.0 * area_total) : t
                frc = shear_fac .* (Dt * shr)
                fo = o2_e' * frc
                out[1:2] = [fo[1], fo[2]]
                out[3] = 0.0
            end
            
            idat = inspector(idat, e, conn, ecoords, copy(out), copy(XYZ))
        end
    end
    
    return idat
end

"""
    make_q4rnt(integdomain, material; mcsys=nothing)

Factory function to create a Q4RNT shell element.
"""
function make_q4rnt(integdomain, material; mcsys=nothing)
    if mcsys === nothing
        mcsys = (zeros(3, 3), zeros(3, 3), zeros(3, 3))
    end
    return FEMMShellQ4RNT(integdomain, mcsys, material)
end

end # module FEMMShellQ4RNTModule
