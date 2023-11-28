using Gridap, IterativeSolvers  
include("grids.jl")


function Vector(n::Int, U::TrialFESpace)
    num_dofs = num_free_dofs(U)
    Fefunction = FEFunction(U, zeros(num_dofs))
    [Fefunction for _ in 1:n]
end


function solver(temporalgrids, fes, j⃗, initials)
    t = temporalmesh(temporalgrids...)
    τ = diff(t)
    (X, Y), (V, U), (~, dΩ) = fes
    mass((v, p, u), (w, q, φ)) = ∫(v⋅w + u*φ)dΩ
    stiff((v, p, u), (w, q, φ)) = ∫(∇(v)⊙∇(w) - (∇⋅w)*p - (j⃗⋅w)*u + (∇⋅v)*q + (v⋅∇(u))*φ + ∇(u)⋅∇(φ))dΩ
    M = assemble_matrix(mass, X, Y)
    S = assemble_matrix(stiff, X, Y)
    v₀, u₀ = initials
    A = similar(M)
    C = similar(M)
    num_dofs = num_free_dofs(X)
    b = zeros(num_dofs)
    solution = zeros(num_dofs)
    vs = Vector(length(t), V)
    us = Vector(length(t), U)
    vs[1] = interpolate_everywhere(v₀, V)
    us[1] = interpolate_everywhere(u₀, U)
    @inbounds for n in eachindex(τ)
        conv((v, p, u), (w, q, φ)) = ∫(w⊙((∇(v))'⋅vs[n]))dΩ
        C .= assemble_matrix(conv, X, Y)
        A .= M + τ[n] * (S + C)
        ℓ((w, q, φ)) = ∫(vs[n] ⋅ w + us[n] * φ)dΩ
        b .= assemble_vector(ℓ, Y)
        cg!(solution, A, b)
        vs[n+1], _, us[n+1] = FEFunction(X, solution)
    end
    vs, us
end