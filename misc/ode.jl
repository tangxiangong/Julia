function solver(T::Float64, u₀::Float64, τ::Float64=0.01)
    t = collect(0.0:τ:T)
    n = length(t)
    u = Vector{Float64}(undef, n)
    u[1] = u₀
    @inbounds for k in 2:n
        u[k] = (u[k-1] + t[k] * τ)/(1+τ)
    end
    return u
end
