module example

using LinearAlgebra
using SparseArrays

struct ThreadBuffer{KVT}
    kcolumns::KVT
    uv::SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true} 
    result::Vector{Float64}
end

function parmul!(R, tbs)
    tasks = [];
    for th in 1:length(tbs)
        tb = tbs[th]
        push!(tasks, Threads.@spawn begin 
            mul!(tb.result, tb.kcolumns, tb.uv)
        end);
    end
    Threads.wait(tasks[1]);
    R .= tbs[1].result
    for th in 2:length(tbs)
        Threads.wait(tasks[th]);
        R .+= tbs[th].result
    end
    R
end

function parloop!(N)
    K = sprand(N, N, 0.1)
    U1 = rand(N)
    R = fill(0.0, N)
    nth = Base.Threads.nthreads()
    @info "$nth threads used"
    chunk = Int(floor(N / nth))
    threadbuffs = ThreadBuffer[];
    for th in 1:nth
        colrange = th < nth ? (chunk*(th-1)+1:chunk*(th)+1-1) : (chunk*(th-1)+1:length(U1))
        push!(threadbuffs, ThreadBuffer(view(K, :, colrange), view(U1, colrange), deepcopy(U1)));
        # push!(threadbuffs, ThreadBuffer(K[:, colrange], view(U1, colrange), deepcopy(U1)));
    end

    R .= parmul!(R, threadbuffs)
    @show norm(R - K*U1)
    true 
end

end # module
nothing

using .example; example.parloop!(10000)




