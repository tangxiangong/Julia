using Gridap

n = 100
domain = (0, 1, 0, 1)
partition = (n, n)
model = CartesianDiscreteModel(domain, partition)

dim = 2
order = 2
reffeᵤ = ReferenceFE(lagrangian, VectorValue{dim, Float64}, order)
V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags="boundary")

reffeₚ = ReferenceFE(lagrangian, Float64, order-1; space=:P)
Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean)

U = TrialFESpace(V, VectorValue(0, 0))
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

degree = order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ, degree)