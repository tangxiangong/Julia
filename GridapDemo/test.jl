using Gridap, GridapGmsh
import GridapGmsh:gmsh


function spatial_mesh(height::Real=1, width::Real=1, lc::Float64=0.01, path::String="./", name::String="model")
    gmsh.initialize()
    gmsh.model.add(name)
    
    p₁ = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc)
    p₂ = gmsh.model.geo.addPoint(width, 0.0, 0.0, lc)
    p₃ = gmsh.model.geo.addPoint(width , height, 0.0, lc)
    p₄ = gmsh.model.geo.addPoint(0.0, height, 0.0, lc)

    l₁ = gmsh.model.geo.addLine(p₁, p₂)
    l₂ = gmsh.model.geo.addLine(p₂, p₃)
    l₃ = gmsh.model.geo.addLine(p₃, p₄)
    l₄ = gmsh.model.geo.addLine(p₄, p₁)

    curve = gmsh.model.geo.addCurveLoop([l₁, l₂, l₃, l₄])
    surface = gmsh.model.geo.addPlaneSurface([curve])
    
    gmsh.model.geo.synchronize()

    dirichlet = gmsh.model.addPhysicalGroup(1, [l₁, l₂, l₃, l₄])
    gmsh.model.addPhysicalGroup(2, [surface], 1)
    gmsh.model.setPhysicalName(1, dirichlet, "boundary")

    gmsh.model.mesh.generate(2)

    file_path = string(path, name, ".msh")
    gmsh.write(file_path)
    gmsh.finalize()

    return file_path
end

function temporal_mesh(T, τ)
    t = collect(0:τ:T)
    τs = diff(t)
    return t, τs
end

function solver(T::Real, height::Real, width::Real, τ::Float64, h::Float64, f::Function, u₀::Function; mesh::Union{String, Nothing}=nothing, order::Int=1, degree::Union{Int, Nothing}=nothing)
    if isnothing(degree)
        degree = order + 2
    end
    if isnothing(mesh)
        nx, ny = round(Int, width/h), round(Int, height/h)
        model = CartesianDiscreteModel((0, width, 0, height), (nx, ny))
    else
        model = GmshDiscreteModel(mesh)
    end

    reffe = ReferenceFE(lagrangian, Float64, order)
    V = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags="boundary")
    U = TrialFESpace(V, x->0)

    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    mass_a(u, v) = ∫(u * v) * dΩ
    stiff_a(u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ

    # 质量矩阵
    mass = assemble_matrix(mass_a, U, V)
    # 刚度矩阵
    stiff = assemble_matrix(stiff_a, U, V)
    # 自由度
    num_dofs = num_free_dofs(U)
    t, τs = temporal_mesh(T, τ)
    num_time_stepping = length(t)-1

    A = similar(mass)
    rhs = zeros(num_dofs)
    U = zeros(num_dofs, num_time_stepping)
    u_pre = interpolate_everywhere(u₀, V)

    for n in 1:num_time_stepping
        A .= mass + τs[n] * stiff
        ft = x -> τs[n] * f(x, t[n+1])
        L(v) = ∫(ft * v) * dΩ + ∫(u_pre * v) * dΩ
        rhs .= assemble_vector(L, V)
        U[:, n] .= A\rhs
        u_pre = FEFunction(V, @view U[:, n])
    end
    U
end

# domain = (0, 1, 0, 1)
# partition = (10, 10)
# model = CartesianDiscreteModel(domain, partition)
# reffe = ReferenceFE(lagrangian, Float64, 1)
# V = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags="boundary")
# U = TrialFESpace(V, x->0)
# Ω = Triangulation(model)
# dΩ = Measure(Ω, 3)    

# f(x) = 2 * π^2 * sin(π*x[1]) * sin(π*x[2])
# a(u, v) = ∫(∇(v) ⊙ ∇(u)) * dΩ
# L(v) = ∫(v * f) * dΩ

# # op = AffineFEOperator(a, L, U, V)
# fs = assemble_vector(L, V)
# # num_free_dofs(U)
# # assemble_matrix(a, U, V) == get_matrix(op)
# # assemble_vector(L, V) == get_vector(op)
# # r = rand(81);
# fh = FEFunction(U, fs)
# # get_free_dof_values(rf) == r
# fh(Point(1/2,1/2))
f(x, t) = (2 * π^2 * sinpi(x[1]) * sinpi(x[2])) * t
u₀(x) = sinpi(x[1]) * sinpi(x[2])
@btime solver(1, 1, 1, 0.01, 0.001, f, u₀);