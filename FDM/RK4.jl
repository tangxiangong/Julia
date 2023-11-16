function solve(T, f, y₀, τ)
    t = collect(0:τ:T)
    y = zeros(length(t))
    y[1] = y₀
    @inbounds for n in firstindex(t)+1:lastindex(t)
        k₁ = τ * f(t[n-1], y[n-1], T)
        k₂ = τ * f(t[n-1]+τ/2, y[n-1]+k₁/2, T)
        k₃ = τ * f(t[n-1]+τ/2, y[n-1]+k₂/2, T)
        k₄ = τ * f(t[n], y[n-1]+k₃, T)
        y[n] = y[n-1] + (k₁ + 2k₂ + 2k₃ + k₄)/6
    end
    t, y
end

function solve!(t, y, T, f, y₀, τ)
    t[1] = 0.0
    y[1] = y₀
    @inbounds for n in firstindex(t)+1:lastindex(t)
        t[n] = t[n-1] + τ
        k₁ = τ * f(t[n-1], y[n-1], T)
        k₂ = τ * f(t[n-1]+τ/2, y[n-1]+k₁/2, T)
        k₃ = τ * f(t[n-1]+τ/2, y[n-1]+k₂/2, T)
        k₄ = τ * f(t[n], y[n-1]+k₃, T)
        y[n] = y[n-1] + (k₁ + 2k₂ + 2k₃ + k₄)/6
    end
    nothing
end


T = 10.0
y₀ = 1.0
τ = 1/256
f(t, x, T) = x + exp(-t/T) * t * sin(t)
t = zeros(2561)
y = zeros(2561)
@btime solve(T, f, y₀, τ);
@benchmark solve!(t, y, T, f, y₀, τ)