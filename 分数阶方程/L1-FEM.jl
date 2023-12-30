using Gridap, GridapGmsh, IterativeSolvers
import GridapGmsh: gmsh
import SpecialFunctions: gamma


function spatialmesh(domain::NTuple{4, <:Real}, h::Float64; path::String = "./model/", name::String = "model")
    a, b, c, d = domain
    @assert a < b && c < d
    @assert h > 0   
    !isdir(path) && (mkdir(path))
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

function FESpaces(mesh; order::Int=1)
    model = GmshDiscreteModel(mesh)
    reffe = ReferenceFE(lagrangian, Float64, order)
    V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags="boundary")
    U = TrialFESpace(V, 0.0)
    Ω = Triangulation(model)
    dΩ = Measure(Ω, order+1)
    (U, V), (Ω, dΩ)
end

function solver(tempargs, fes, initial, source)
    T, α, τ = tempargs
    constant = 1/(τ^α * gamma(2-α))
    
    t = collect(0:τ:T)
    
    (U, V), (~, dΩ) = fes
    mass(u, v) = ∫(u*v)dΩ
    stiff(u, v) = ∫(∇(u)⋅∇(v))dΩ
    M = assemble_matrix(mass, U, V)
    S = assemble_matrix(stiff, U, V)
    A = constant * M + S

    u₀ = initial
    init()

    num_dofs = num_free_dofs(U)
    b = zeros(num_dofs)
    solution = zeros(num_dofs)
    @inbounds for n in firstindex(t)+1:lastindex(t)
        ℓ((w, q, φ)) = ∫(vs[n] ⋅ w + us[n] * φ)dΩ
        b .= assemble_vector(ℓ, Y)
        cg!(solution, A, b)
        vs[n+1], _, us[n+1] = FEFunction(X, solution)
    end
    vs, us
end