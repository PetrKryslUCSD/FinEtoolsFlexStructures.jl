using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield

# Unsupported 4x4 patch of the unit sphere, interior nodes moved off the grid
function sphere_patch(n = 4, a = 0.5)
    h = 2a / n
    id(i, j) = j * (n + 1) + i + 1
    xyz = zeros((n + 1)^2, 3)
    for j in 0:n, i in 0:n
        x, y = -a + i * h, -a + j * h
        if 0 < i < n && 0 < j < n
            x += 0.15h * sin(2.1i + 1.3j)
            y += 0.15h * cos(1.7i - 2.3j)
        end
        xyz[id(i, j), :] .= (x, y, sqrt(1 - x^2 - y^2))
    end
    conn = zeros(Int, 2n^2, 3)
    e = 0
    for j in 0:n-1, i in 0:n-1
        conn[e += 1, :] .= (id(i, j), id(i + 1, j), id(i + 1, j + 1))
        conn[e += 1, :] .= (id(i, j), id(i + 1, j + 1), id(i, j + 1))
    end
    return FENodeSet(xyz), FESetT3(conn)
end

# max|K v| / (max|K| max|v|) for the six rigid-body motions of the free patch
function rigid_body_residuals()
    fens, fes = sphere_patch()
    accepttodelegate(fes, FESetShellT3())
    femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), 0.1),
                                    MatDeforElastIso(DeforModelRed3D, 1.0e7, 0.3))
    femm.drilling_stiffness_scale = 0.0  # the drilling springs act on the nodal drilling rotation by design
    geom0 = NodalField(fens.xyz)
    dchi = NodalField(zeros(count(fens), 6))
    numberdofs!(dchi)
    FEMMShellT3FFModule.associategeometry!(femm, geom0)
    K = FEMMShellT3FFModule.stiffness(femm, SysmatAssemblerFFBlock(nfreedofs(dchi)), geom0,
                                      NodalField(zeros(count(fens), 3)), initial_Rfield(fens), dchi)
    x0 = vec(sum(fens.xyz, dims = 1)) ./ count(fens)
    res = zeros(6)
    for m in 1:6
        v = zeros(nfreedofs(dchi))
        for k in 1:count(fens)
            d = dchi.dofnums[k, :]
            if m <= 3                       # unit translation along x, y, z
                v[d[m]] = 1.0
            else                            # unit rotation about x, y, z through x0
                w = zeros(3); w[m - 3] = 1.0
                v[d[1:3]] .= cross(w, fens.xyz[k, :] .- x0)
                v[d[4:6]] .= w
            end
        end
        res[m] = maximum(abs, K * v) / (maximum(abs, K) * maximum(abs, v))
    end
    return res
end

println("as released: ", round.(rigid_body_residuals(), sigdigits = 2))

# The same function with the two terms corrected
@eval FEMMShellT3FFModule function _transfmat_a_to_e!(T, A_Es, gradN_e)
    T .= zero(eltype(T))
    for i = 1:__NN
        roffst = (i - 1) * __NDOF
        iA_E = A_Es[i]
        iA_E_33 = iA_E[3, 3]
        for cl = 1:3, rw = 1:3
            T[roffst+rw, roffst+cl] = iA_E[rw, cl]
        end
        for cl = 1:2, rw = 1:2
            T[roffst+3+rw, roffst+3+cl] =
                iA_E[rw, cl] - (1 / iA_E_33) * iA_E[rw, 3] * iA_E[3, cl]   # was iA_E[cl, 3]
        end
        m1 = (1 / iA_E_33) * iA_E[1, 3]
        m2 = (1 / iA_E_33) * iA_E[2, 3]
        for j = 1:__NN
            coffst = (j - 1) * __NDOF
            jA_E = A_Es[j]                                                 # was iA_E
            for k = 1:3
                a3 = 1 / 2 * (jA_E[2, k] * gradN_e[j, 1] - jA_E[1, k] * gradN_e[j, 2])
                T[roffst+4, coffst+k] += m1 * a3
                T[roffst+5, coffst+k] += m2 * a3
            end
        end
    end
    return T
end

println("corrected:   ", round.(Base.invokelatest(rigid_body_residuals), sigdigits = 2))
