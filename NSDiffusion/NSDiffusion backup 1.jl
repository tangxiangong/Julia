using LinearAlgebra, MKL, Gridap, GridapGmsh, IterativeSolvers
import GridapGmsh: gmsh

#=
    网格剖分
=#
function temporalmesh(time_length::T, τ::Float64, α::Float64) where T <: Real
    @assert time_length > 0 && τ > 0 && 1/2 < α < 1

    t = Vector{Float64}()
    sizehint!(t, round(Int, time_length/α))
    append!(t, 0.0)
    append!(t, time_length*(τ/time_length)^(1/(1-α)))
    last_time = t[end]
    while last_time < time_length
        time_step = τ * (last_time/time_length)^α
        last_time += time_step
        append!(t, last_time)
    end
    t[end] = time_length
    t
end

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

#= 
    有限元空间
=#
function FESpaces(mesh="./model.msh", orders::NTuple{2, Int}=(2, 1))
    model = GmshDiscreteModel(mesh)
    # Reference Finite Elements
    reffeᵥ = ReferenceFE(lagrangian, VectorValue{2, Float64}, orders[1]; space=:P)
    reffeₚ = ReferenceFE(lagrangian, Float64, orders[1]-1; space=:P)
    reffeᵤ = ReferenceFE(lagrangian, Float64, orders[2])
    # Test Fuction Spaces
    W = TestFESpace(model, reffeᵥ, conformity=:H1, dirichlet_tags="boundary")
    Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean)
    Z = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags="boundary")
    # Trial Function Spaces
    V = TrialFESpace(W, VectorValue(0.0, 0.0))
    P = TrialFESpace(Q)
    U = TrialFESpace(Z, 0.0)
    # MultiFiled Trial/Test Function Spaces
    X = MultiFieldFESpace([V, P, U])
    Y = MultiFieldFESpace([W, Q, Z])
    # Triangulation and Measurable Computional Domain
    Ω = Triangulation(model)
    dΩ = Measure(Ω, maximum(orders)+1)
    (X, Y), (V, U), (Ω, dΩ)
end

#=
    辅助函数
=#
function Vector(n::Int, U::TrialFESpace)
    num_dofs = num_free_dofs(U)
    Fefunction = FEFunction(U, zeros(num_dofs))
    [Fefunction for _ in 1:n]
end

#=
    质量、刚度矩阵组装
=#
function assemble_mass_stiff_matrix(j⃗, X, Y, dΩ)
    mass((v, p, u), (w, q, φ)) = ∫(v⋅w + u*φ)dΩ
    stiff((v, p, u), (w, q, φ)) = ∫(∇(v)⊙∇(u) - (∇⋅w)*p -(j⃗⋅w)*u + (∇⋅v)*q + (v⋅∇(u))*φ + ∇(u)⋅∇(φ))dΩ
    M = assemble_matrix(mass, X, Y)
    S = assemble_matrix(stiff, X, Y)
    M, S
end

#=
    右端项组装及时间步进求解
=#
"""
    solver(temporalgrids, fes, j⃗, initials, sources)

# 参数
- `temporalgrids`: (T, τ, α)
- `fes`: ((X, Y), (V, U), (Ω, dΩ))
- `initials`: (v₀, u₀)
- `sources`: (f, g)
"""
function solver(temporalgrids, fes, j⃗, initials, sources)
    t = temporalmesh(temporalgrids...)
    τ = diff(t)
    (X, Y), (V, U), (~, dΩ) = fes
    M, S = assemble_mass_stiff_matrix(j⃗, X, Y, dΩ)
    v₀, u₀ = initials
    f, g = sources
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
        # 对流项
        conv((v, p, u), (w, q, φ)) = ∫(w⊙((∇(v))'⋅vs[n]))dΩ
        C .= assemble_matrix(conv, X, Y)
        # 总系数矩阵
        A .= M + τ[n] * (S + C)
        fⁿ(x) = τ[n] * f(x, t[n+1])
        gⁿ(x) = τ[n] * g(x, t[n+1])
        # 右端泛函
        ℓ((w, q, φ)) = ∫(fⁿ⋅w + gⁿ*φ + vs[n]⋅w + us[n]*φ)dΩ
        # 总右端向量
        b .= assemble_vector(ℓ, Y)
        # 求解
        cg!(solution, A, b)
        # 有限元函数
        vs[n+1], ~, us[n+1] = FEFunction(X, solution)
    end
    vs, us
end

""""
    计算 L2 和 H1 误差

若两个函数均为有限元空间中的函数，则对应的基函数必须一致， 否则，需要进行插值处理。
```julia
# 若 u₁ u₂ 分别对应于剖分 Ω₁ （疏） Ω₂ (密)， 和空间 U₁ U₂。
# 将函数 u₁ 插值到空间 U₂ 中:
̄u₁ = interpolate(u₁, U₂)
# 计算误差
norm_error(̄u₁, u₂, dΩ)
```
"""
function norm_error(u₁, u₂, dΩ, scarlar_or_vector="scalar")
    @assert scarlar_or_vector in ["scalar", "vector"]

    e = u₁ - u₂
    if scarlar_or_vector == "scalar"
        L²_error = ∫(e * e)dΩ |> sum |> sqrt
        H¹_error = ∫(e * e + ∇(e) ⋅ ∇(e))dΩ |> sum |> sqrt
    else
        L²_error = ∫(e ⋅ e)dΩ |> sum |> sqrt
        H¹_error = ∫(e ⋅ e + ∇(e) ⊙ ∇(e))dΩ |> sum |> sqrt
    end
    L²_error, H¹_error
end



#= 示例 =#

domain = (0, 1, 0, 1)
h = 1/128
mesh = spatialmesh(domain, h)
fes = FESpaces(mesh)
initals = (x->VectorValue(x[1], x[2]), x->sinpi(x[1]) * sinpi(x[2]))
sources = ((x, t)->VectorValue(exp(-t)*x[1], exp(-t)*x[2]), (x, t)->exp(-t)*norm(x))
j⃗ = VectorValue(1.0, 2.0)
v₁, u₁ = solver((1.0, 1/32, 0.6), fes, j⃗, initials, sources)
v₂, u₂ = solver((1.0, 1/64, 0.6), fes, j⃗, initials, sources)
# 误差
~, ~, (Ω, dΩ) = fes
errorᵥ = norm_error(v₁, v₂, dΩ, "vector") 
errorᵤ = norm_error(u₁, u₂, dΩ, "scalar")
