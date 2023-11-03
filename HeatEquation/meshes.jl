import GridapGmsh:gmsh


function temporal_mesh(time_length::T, τ::Float64, α::Float64) where T <: Real
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


function spatial_mesh(domain::NTuple{4, <:Real}, h::Float64, path::String = "./", name::String = "model")
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