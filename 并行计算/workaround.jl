using .Threads
using LoopVectorization

function threaded_workaround(n::Integer)
    partial = zeros(Int64, nthreads())
    @threads for _ in 1:n
        x = rand()
        y = rand()
        if x^2 + y^2 <= 1
            @inbounds partial[threadid()] += 1
        end
    end
    return 4*sum(partial)/n
end

function serial(n::Integer)
    count = 0
    for _ in 1:n
        x = rand()
        y = rand()
        if x^2 + y^2 <= 1
            count += 1
        end
    end
    return 4*count/n
end

function sqrt_sum(A)
    s = zero(eltype(A))
    for i in eachindex(A)
        @inbounds s += sqrt(A[i])
    end
    return s
end

function threaded_sqrt_sum_workaround(A)
    partial = zeros(eltype(A), nthreads())
    @threads for i in eachindex(A)
        @inbounds partial[threadid()] += sqrt(A[i])
    end
    s = zero(eltype(A))
    for i in eachindex(partial)
        s += partial[i]
    end
    return s
end

function threaded_sqrt_sum_atomic(A)
    s = Atomic{eltype(A)}(zero(eltype(A)))
    @threads for i in eachindex(A)
        @inbounds atomic_add!(s, sqrt(A[i]))
    end
    return s[]
end