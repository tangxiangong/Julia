import SpecialFunctions: gamma


function mesh(time_length::T, τ::Float64, β::Float64) where T <: Real
    @assert time_length > 0 && τ > 0 && 1/2 < β < 1

    t = Vector{Float64}()
    sizehint!(t, round(Int, time_length/β))
    append!(t, 0.0)
    append!(t, time_length*(τ/time_length)^(1/(1-β)))
    last_time = t[end]
    while last_time < time_length
        time_step = τ * (last_time/time_length)^β
        last_time += time_step
        append!(t, last_time)
    end
    t[end] = time_length
    t
end

function ∂ᵅ(α, func, t, τ, β=0.6)
    ts = mesh(t, τ, β)
    τₙ = t[end] - t[end-1]
    c = 1/gamma(2-α)
    result = 0.0
    @inbounds for k in firstindex(ts)+1:lastindex(ts)-1
        ω = (t-t[k+1])^(1-α) - 2(t-t[k])^(1-α) + (t-t[k-1])^(1-α)
        result += ω * func(t[k])
    end
    result += τₙ^(1-α)*func(t) + ((t-t[begin+1])^(1-α)-t^(1-α))*func(0)
    result *= c
    result
end

@doc raw"""
solve equation
```math
\partial_t^{\alpha} u(t) - \lambda u(t) = f(t)
```
"""
function solver(α, T, λ, u₀, f, params)
    β, τ = params
    t = mesh(T, τ, β)
    u = zeros(length(t))
    u[1] = u₀
    constant = 1/gamma(2-α)
    for n in firstindex(t)+1:lastindex(t)
        tₙ = t[n]
        τₙ = tₙ - t[n-1]
        b = tₙ^(1-α) - (tₙ-t[begin+1])^(1-α)
        for k in 2:n-1
            c = (tₙ-t[k+1])^(1-α) - 2(tₙ-t[k])^(1-α) + (tₙ - t[k-1])^(1-α)
            b -= c * u[k]
        end
        b *= constant
        b += f(tₙ)
        u[n] = (τₙ^(1-α)*constant - λ)\b
    end
    t, u
end

f(t) = t
t, u = solver(0.7, 1, 1, 0, f, (0.7, 1/32));

using Plots
plot(t, u)