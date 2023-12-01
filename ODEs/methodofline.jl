using LinearAlgebra, SparseArrays, BenchmarkTools, IterativeSolvers

function Euler(T, u₀, f, τ, M)
    t = collect(0:τ:T)
    # u = zeros(length(u₀), length(t))
    # @inbounds u[:, 1] .= u₀
    A = I - τ * M
    b = similar(u₀)
    x = similar(b)
    x .= u₀
    @inbounds for n in firstindex(t)+1:lastindex(t)
        b .= τ * f(t[n]) + x
        cg!(x, A, b)
    end
    t, x
end 

h = 1/128
x = h:h:1-h
u₀ = sinpi.(x)
f(t) = (π^2-1)*exp(-t)*sinpi.(x)
M = spdiagm(0=>-2*ones(length(x)), 1=>ones(length(x)-1), -1=>ones(length(x)-1))/h^2
t, u = Euler(1, u₀, f, 1/64, M)
uₑ(t) = exp(-t)*sinpi.(x)
norm(u[:, end]-uₑ(t[end]))

@btime Euler(1, u₀, f, 1/256, M);  # 2.1 ms
