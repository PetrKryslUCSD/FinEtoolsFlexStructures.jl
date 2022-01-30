module test_loops

using LinearAlgebra
using SparseArrays
using SparseMatricesCSR
using ThreadedSparseCSR
using BenchmarkTools

struct ThreadBuff{IT,JT,VT}
    # mrange::UnitRange{Int64}
    I::IT
    J::JT
    V::VT
    result::Vector{Float64}
end

function coomul!(A, IB, JB, VB, C)
    nz = length(IB)
    @assert length(JB) == nz
    @assert length(VB) == nz
    A .= zero(eltype(A))
    @inbounds @simd for i in 1:nz
        ir = IB[i]
        jc = JB[i]
        A[ir] += VB[i] * C[jc]
    end
end

function parcoomul!(R, IB, JB, VB, tbs, U)
    function innercoomul!(A, IB, JB, VB, C)
        nz = length(IB)
        @assert length(JB) == nz
        @assert length(VB) == nz
        A .= zero(eltype(A))
        @inbounds @simd for i in 1:nz
            ir = IB[i]
            jc = JB[i]
            A[ir] += VB[i] * C[jc]
        end
    end
    tasks = []
    for th in 1:length(tbs)
        push!(tasks, Threads.@spawn begin
            innercoomul!(tbs[th].result, tbs[th].I, tbs[th].J, tbs[th].V, U)
        end)
    end
    Threads.wait(tasks[1])
    R .= tbs[1].result
    for th in 2:length(tbs)
        Threads.wait(tasks[th])
        R .+= tbs[th].result
    end
    R
end

function parcsrmul!(y::AbstractVector, A::SparseMatrixCSR, v::AbstractVector)
  A.n == size(v, 1) || throw(DimensionMismatch())
  A.m == size(y, 1) || throw(DimensionMismatch())
  fill!(y, zero(eltype(y)))
  o = getoffset(A)
  Threads.@threads for row = 1:size(y, 1)
    @inbounds for nz in nzrange(A,row)
      col = A.colval[nz]+o
      y[row] += A.nzval[nz]*v[col]
    end
  end
  return y
end

function _execute()
    N = 2000000
    K = sprand(N, N, 0.00004)
    K = K + K'
    # @show K
    # K = sparse([5, 5, 10, 4, 3, 6, 9, 10, 1, 2, 8, 4, 5, 10, 4, 2, 4, 8], [1, 2, 2, 3, 4, 4, 4, 4, 5, 5, 5, 6, 8, 8, 9, 10, 10, 10], [2.27568e-01, 1.10446e-01, 5.29148e-01, 3.24832e-01, 3.24832e-01, 6.16751e-03, 8.60792e-01, 1.01654e+00, 2.27568e-01, 1.10446e-01, 7.70209e-01, 6.16751e-03, 7.70209e-01, 9.17418e-01, 8.60792e-01, 5.29148e-01, 1.01654e+00, 9.17418e-01], 10, 10)
    U = rand(N)
    R = zeros(N)

    @show nth = Base.Threads.nthreads()
    nz = nnz(K)
    @show chunk = Int(floor(nz / nth))

    IK, JK, VK = findnz(K)

    threadbuffs = ThreadBuff[]
    for th in 1:nth
        mrange = th < nth ? (chunk*(th-1)+1:chunk*(th)+1-1) : (chunk*(th-1)+1:nz)
        push!(threadbuffs, ThreadBuff(view(IK, mrange), view(JK, mrange), view(VK, mrange), deepcopy(R)))
    end

    Kr = SparseMatricesCSR.sparsecsr(IK, JK, VK, N, N)

    oU = deepcopy(U)
    oU .= rand()

    U .= oU
    mul!(R, K, U)
    R1 = deepcopy(R)

    U .= oU
    coomul!(R, IK, JK, VK, U)
    @show norm(R - R1) / norm(R1)

    U .= oU
    parcoomul!(R, IK, JK, VK, threadbuffs, U)
    @show norm(R - R1) / norm(R1)

    U .= oU
    parcsrmul!(R, Kr, U)
    @show norm(R - R1) / norm(R1)

    @btime ThreadedSparseCSR.tmul!($R, $Kr, $U)
    @btime ThreadedSparseCSR.bmul!($R, $Kr, $U)
    @btime mul!($R, $K, $U)
    @btime $coomul!($R, $IK, $JK, $VK, $U)
    @btime $parcoomul!($R, $IK, $JK, $VK, $threadbuffs, $U)
    @btime SparseMatricesCSR.mul!($R, $Kr, $U)
    @btime $parcsrmul!($R, $Kr, $U)

    nothing
end

function allrun()
    _execute()
    return true
end # function allrun

end # module
nothing

 using .Main.test_loops
 test_loops.allrun()                                                 







