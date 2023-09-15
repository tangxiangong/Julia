function sum_single(x)
    s = 0.0
    for i in x
        s += i
    end
    return s
end

function sum_multi_bad(x::Vector{Float64})
    s = 0.0
    Threads.@threads for i in x
        s += i
    end
    return s
end

function sum_multi_good(x)
    chunks = Iterators.partition(x, length(x) ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn sum_single(chunk)
    end
    chunk_sums = fetch.(tasks)
    return sum(chunk_sums)
end


x = rand(1_000_000)
sum_single(x)
sum_multi_bad(x)
sum_multi_good(x)
# using BenchmarkTools

# @benchmark 