using Distributions, Plots, LoopVectorization

function skewed_stable_rand(alpha)
    c = π/2
    w = rand(Exponential())
    γ = -c + 2 * c * rand()
    s₁ = sin(alpha * (γ + c))
    s₂ = (cos(γ))^(1/alpha)
    s₃ = cos(γ - alpha * (γ + c))
    s₄ = (s₃ / w)^((1-alpha)/alpha)
    rnd = s₄ * s₁ / s₂
    return rnd
end

function subordinator(length_t, param, τ=0.01)
    t = collect(0.0:τ:length_t)
    n = length(t)
    noise = zeros(Float64, n-1)
    @inbounds @simd for i in eachindex(noise)
        noise[i] = τ^(1/param) * skewed_stable_rand(param)
    end
    x = cumsum([0.0; noise])
    return t, x, noise
end
t, s, noise = subordinator(100, 0.7)

# plot(t, s)

# function inv_subordinator(length_t, param, τ=0.01)
#     t, s = subordinator(length_t, param)
#     count = 1
#     while s[end] < length_t
#         t, s = subordinator(2^count*length_t, param)
#         count += 1
#         @show s[end]
#     end
    
# end