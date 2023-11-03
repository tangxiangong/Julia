using Gridap, GridapGmsh
using IterativeSolvers

include("meshes.jl")

function generateFEMs(mesh="./model.msh", order::Int=1)
    model = GmshDiscreteModel(mesh)
    reffe = ReferenceFE(lagrangian, Float64, order)
    V = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags="boundary")
    U = TrialFESpace(V, x->0)
    Ω = Triangulation(model)
    dΩ = Measure(Ω, order+1)
    return U, V, dΩ
end


function solver(T, τ, α, FEMs, f::Function, u₀::Function)
    U, V, dΩ = FEMs 
    mass_bilinear(u, v) = ∫(u * v)dΩ
    stiff_bilinear(u, v) = ∫(∇(u) ⊙ ∇(v))dΩ
    mass = assemble_matrix(mass_bilinear, U, V)
    stiff = assemble_matrix(stiff_bilinear, U, V)
    num_dofs = num_free_dofs(U)
    t = temporal_mesh(T, τ, α)
    τs = diff(t)
    num_time_stepping = length(t) - 1
    A = similar(mass)
    rhs = zeros(num_dofs)
    uₕ = zeros(num_dofs)
    # Uₕ = zeros(num_dofs, num_time_stepping)
    last_uₕ = u₀
    for n in 1:num_time_stepping
        A .= mass + τs[n] * stiff
        ft = x -> τs[n] * f(x, t[n+1])
        ℓ(v) = ∫(ft * v)dΩ + ∫(last_uₕ * v)dΩ
        rhs .= assemble_vector(ℓ, V)
        cg!(uₕ, A, rhs)
        last_uₕ = FEFunction(U, uₕ)
        # last_u = FEFunction(U, @view Uₕ[:, n])
    end
    last_uₕ
end