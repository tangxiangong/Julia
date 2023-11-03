include("meshes.jl")
include("solver.jl")

T = 1
domain = (0, 1, 0, 1)
h = 1/256
τ = 1/32
α = 0.6
u₀(x) = sinpi(x[1]) * sinpi(x[2])
f(x, t) = (2 * π^2 - 1) * sinpi(x[1]) * sinpi(x[2]) * exp(-t)

mesh = spatial_mesh(domain, h)
FEMs = generateFEMs()

uₕ₁ = solver(T, τ, α, FEMs, f, u₀)
uₕ₂ = solver(T, τ/2, α, FEMs, f, u₀)

uₑ(x, t) = exp(-t) * sinpi(x[1]) * sinpi(x[2])
u = x -> uₑ(x, T)

~, ~, dΩ = FEMs

e₁ = ∫((u - uₕ₁)*(u - uₕ₁))dΩ |> sum |> sqrt
e₂ = ∫((u - uₕ₂)*(u - uₕ₂))dΩ |> sum |> sqrt

log(e₁/e₂)/log(2)

# u₂ = solver(T, τ/2, α, FEMs, f, u₀)
# u₃ = solver(T, τ/4, α, FEMs, f, u₀)
# ū₁ = interpolate_everywhere(Interpolable(u₁), U₂)
# e₁ = ∫((u₂ - ū₁) * (u₂ - ū₁))dΩ₂ |> sum |> sqrt

# ū₂ = interpolate_everywhere(Interpolable(u₂), U₃)
# e₂ = ∫((u₃ - ū₂) * (u₃ - ū₂))dΩ₃ |> sum |> sqrt