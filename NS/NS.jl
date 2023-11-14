using Gridap, GridapGmsh, IterativeSolvers
import GridapGmsh: gmsh


temporalmesh(T, τ) = 0.0:τ:T

function spatialmesh(domain::NTuple{4, <:Real}, h::Float64; path::String = "./", name::String = "model")
    a, b, c, d = domain
    @assert a < b && c < d
    @assert h > 0   
    gmsh.initialize()
    gmsh.model.add(name)
    p₁ = gmsh.model.geo.addPoint(a, c, 0, h)
    p₂ = gmsh.model.geo.addPoint(b, c, 0, h)
    p₃ = gmsh.model.geo.addPoint(b, d, 0, h)
    p₄ = gmsh.model.geo.addPoint(a, d, 0, h)
    ℓ₁ = gmsh.model.geo.addLine(p₁, p₂)
    ℓ₂ = gmsh.model.geo.addLine(p₂, p₃)
    ℓ₃ = gmsh.model.geo.addLine(p₃, p₄)
    ℓ₄ = gmsh.model.geo.addLine(p₄, p₁)
    curve_loop = gmsh.model.geo.addCurveLoop([ℓ₁, ℓ₂, ℓ₃, ℓ₄])
    plane_surface = gmsh.model.geo.addPlaneSurface([curve_loop])
    gmsh.model.geo.synchronize()
    DirichletBC_group = gmsh.model.addPhysicalGroup(1, [ℓ₁, ℓ₂, ℓ₃, ℓ₄])
    gmsh.model.addPhysicalGroup(2, [plane_surface], 1)
    gmsh.model.setPhysicalName(1, DirichletBC_group, "boundary")
    gmsh.model.mesh.generate(2)
    file_path = string(path, name, ".msh")
    gmsh.write(file_path)   
    gmsh.finalize()
    file_path
end

Gridap.:(⊗)(U::FESpace, V::FESpace) = MultiFieldFESpace([U, V])

function FESpaces(; mesh="./model.msh", order::Int=2)
    model = GmshDiscreteModel(mesh)
    # Reference Finite Elements
    reffeᵥ = ReferenceFE(lagrangian, VectorValue{2, Float64}, order; space=:P)
    reffeₚ = ReferenceFE(lagrangian, Float64, order-1; space=:P)
    # Test Fuction Spaces
    Vₕ = TestFESpace(model, reffeᵥ, conformity=:H1, dirichlet_tags="boundary")
    Qₕ = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean)
    # Trial Function Spaces
    Uₕ = TrialFESpace(Vₕ, VectorValue(0.0, 0.0))
    Pₕ = TrialFESpace(Qₕ)
    # MultiFiled Trial/Test Function Spaces
    X = Uₕ ⊗ Pₕ
    Y = Vₕ ⊗ Qₕ
    # Triangulation and Measurable Computional Domain, Quadrature.
    Ωₕ = Triangulation(model)
    dΩₕ = Measure(Ωₕ, order+1)
    (X, Y), (Uₕ, Pₕ), (Ωₕ, dΩₕ)
end

function assemble_mass_stiff_matrix(X, Y, dΩₕ)
    mass((u⃗, p), (v⃗, q)) = ∫(u⃗ ⋅ v⃗)dΩₕ
    stiff((u⃗, p), (v⃗, q)) = ∫(∇(u⃗) ⊙ ∇(v⃗) - p * (∇ ⋅ v⃗) + (∇ ⋅ u⃗) * q)dΩₕ
    M = assemble_matrix(mass, X, Y)
    S = assemble_matrix(stiff, X, Y)
    M, S
end

function solver(temporalgrids, fes, u⃗⁰, f⃗)
    t = temporalmesh(temporalgrids...)
    τ = diff(t)
    (X, Y), (_, __), (___, dΩₕ) = fes
    M, S = assemble_mass_stiff_matrix(X, Y, dΩₕ)
    A = similar(M)
    C = similar(M)
    num_dofs = num_free_dofs(X)
    b = zeros(num_dofs)
    solution = zeros(num_dofs)
    u⃗ⁿ = u⃗⁰
    @inbounds for n in eachindex(τ)
        conv((u⃗, p), (v⃗, q)) = ∫(v⃗ ⊙ ((∇(u⃗))' ⋅ u⃗ⁿ))dΩₕ
        C .= assemble_matrix(conv, X, Y)
        A .= M + τ[n] * (S + C)
        f⃗ⁿ(x) = τ[n] * f⃗(x, t[n+1])
        ℓ((v⃗, q)) = ∫(f⃗ⁿ ⋅ v⃗ + u⃗ⁿ ⋅ v⃗)dΩₕ
        b .= assemble_vector(ℓ, Y)
        cg!(solution, A, b)
        u⃗ⁿ, ____ = FEFunction(X, solution)
    end
    u⃗ⁿ
end


u⃗⁰(x) = VectorValue(-π*sinpi(x[1])*cospi(x[2]), π*cospi(x[1])*sinpi(x[2]))
f⃗(x, t) = VectorValue(π^3*exp(-2t)*sinpi(x[1])*cospi(x[1])-2π^3*exp(-t)*sinpi(x[1])*cospi(x[2]), -2π*exp(-t)*cospi(x[1])*sinpi(x[2])+π^3*exp(-2t)*sinpi(x[2])*cospi(x[2])+2π^3*exp(-t)*cospi(x[1])*sinpi(x[2]))
u⃗ₑ(x, t) = VectorValue(-π*exp(-t)*sinpi(x[1])*cospi(x[2]), π*exp(-t)*cospi(x[1])*sinpi(x[2]))


T = 1.0
τ = 1/32
domain = (0, 1, 0, 1)
h = 1/64

spatialmesh(domain, h)

fes = FESpaces()
@time u⃗ = solver((T, τ), fes, u⃗⁰, f⃗)

(X, Y), (Uₕ, Pₕ), (Ωₕ, dΩₕ) = fes
num_free_dofs(X)