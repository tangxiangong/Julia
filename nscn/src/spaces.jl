using Gridap, GridapGmsh


function FESpaces(mesh; orders::NTuple{2, Int}=(2, 1))
    model = GmshDiscreteModel(mesh)
    reffeᵥ = ReferenceFE(lagrangian, VectorValue{2, Float64}, orders[1]; space=:P)
    reffeₚ = ReferenceFE(lagrangian, Float64, orders[1]-1; space=:P)
    reffeᵤ = ReferenceFE(lagrangian, Float64, orders[2])
    W = TestFESpace(model, reffeᵥ, conformity=:H1, dirichlet_tags="boundary")
    Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean)
    Z = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags="boundary")
    V = TrialFESpace(W, VectorValue(0.0, 0.0))
    P = TrialFESpace(Q)
    U = TrialFESpace(Z, 0.0)
    X = MultiFieldFESpace([V, P, U])
    Y = MultiFieldFESpace([W, Q, Z])
    Ω = Triangulation(model)
    dΩ = Measure(Ω, maximum(orders))
    (X, Y), (V, U), (Ω, dΩ)
end