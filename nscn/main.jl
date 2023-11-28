include("./src/error.jl")
include("./src/grids.jl")
include("./src/solver.jl")
include("./src/spaces.jl")


domain = (0, 1, 0, 1)
h = 1/128
mesh = spatialmesh(domain, h; name="square_128")
fes = FESpaces(mesh)
(X, Y), _, (Ω, dΩ) = fes
initials = (x->VectorValue(x[1], x[2]), x->sinpi(x[1]) * sinpi(x[2]))
j⃗ = VectorValue(1.0, 2.0)
v₁, u₁ = solver((1.0, 1/32, 0.6), fes, j⃗, initials)

# v₂, u₂ = solver((1.0, 1/64, 0.6), fes, j⃗, initials)


# errorᵥ = norm_error(v₁, v₂, dΩ, "vector") 
# errorᵤ = norm_error(u₁, u₂, dΩ, "scalar")

