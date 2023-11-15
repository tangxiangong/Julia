struct ODE
    T::Float64
    f::Function
    y₀::Float64
end

abstract type Solver end

abstract type RKSolver <: Solver end

struct RK4 <: RKSolver
    τ::Float64
end

function solve(ode::ODE, solver::RK4)
    T, f, y₀ = ode.T, ode.f, ode.y₀
    τ = solver.τ
    t = collect(0:τ:T)
    y = zeros(length(t))
    y[1] = y₀
    @fastmath @inbounds for n in firstindex(t)+1:lastindex(t)
        k₁ = τ * f(t[n-1], y[n-1])
        k₂ = τ * f(t[n-1]+τ/2, y[n-1]+k₁/2)
        k₃ = τ * f(t[n-1]+τ/2, y[n-1]+k₂/2)
        k₄ = τ * f(t[n], y[n-1]+k₃)
        y[n] = y[n-1] + (k₁ + 2k₂ + 2k₃ + k₄)/6
    end
    t, y
end

function solve(T, f, y₀, τ)
    t = collect(0:τ:T)
    y = zeros(length(t))
    y[1] = y₀
    @fastmath @inbounds for n in firstindex(t)+1:lastindex(t)
        k₁ = τ * f(t[n-1], y[n-1])
        k₂ = τ * f(t[n-1]+τ/2, y[n-1]+k₁/2)
        k₃ = τ * f(t[n-1]+τ/2, y[n-1]+k₂/2)
        k₄ = τ * f(t[n], y[n-1]+k₃)
        y[n] = y[n-1] + (k₁ + 2k₂ + 2k₃ + k₄)/6
    end
    t, y
end

function solve!(t, y, f, y₀, τ)
    t[1] = 0.0
    y[1] = y₀
    @fastmath @inbounds for n in firstindex(t)+1:lastindex(t)
        t[n] = t[n-1] + τ
        k₁ = τ * f(t[n-1], y[n-1])
        k₂ = τ * f(t[n-1]+τ/2, y[n-1]+k₁/2)
        k₃ = τ * f(t[n-1]+τ/2, y[n-1]+k₂/2)
        k₄ = τ * f(t[n], y[n-1]+k₃)
        y[n] = y[n-1] + (k₁ + 2k₂ + 2k₃ + k₄)/6
    end
    nothing
end



T = 10.0
y₀ = 1.0
τ = 1/256

rk4solver = RK4(τ)
f_(t, x, T) = x + exp(-t/T) * t * sin(t)
ode = ODE(T, (t, x)->f_(t, x, T), y₀)
# @btime solve(ode, rk4solver);   # 2ms
@btime solve(T, (t, x)->f_(t, x, T), y₀, τ);

t = zeros(2561)
y = zeros(2561)
@btime solve!(t, y, (t, x)->f_(t, x, T), y₀, τ);