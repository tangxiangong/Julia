# using Gridap

# domain = (0.0, 1.0, 0.0, 1.0)
# partition = (10, 10)
# model = CartesianDiscreteModel(domain, partition)

# writevtk(model, "test")

# using Gridap, GridapMakie, GLMakie, GridapGmsh
# # using FileIO    

# # # mkdir("models")
# # # mkdir("images")

# # domain = (1,2,0,1)
# # partition = (10,10)

# model = GmshDiscreteModel("square_regular.msh") |> simplexify
# # model = CartesianDiscreteModel(domain,partition) |> simplexify
# Ωₕ = Triangulation(model)
# writevtk(model, "models/model")

# fig = plot(Ωₕ)
# wireframe!(Ωₕ, color=:black, linewidth=2) 
# scatter!(Ωₕ, marker=:star8, markersize=20, color=:blue)
# save("images/model.png", fig)

