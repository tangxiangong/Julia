using Gridap, IterativeSolvers

model = CartesianDiscreteModel((0, 1, 0, 1), (100, 100))

reffe = ReferenceFE(lagrangian, Float64, 1)
V = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags="boundary")
U = TrialFESpace(V, x->0)

Ω = Triangulation(model)
dΩ = Measure(Ω, 2)


f(x) = sinpi(x[1]) * sinpi(x[2])
a(u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ

L(v) = ∫(f * v) * dΩ 

num_dofs = num_free_dofs(U)

op = AffineFEOperator(a, L, U, V)

using BenchmarkTools
@btime solve($op);

A = get_matrix(op);
b = get_vector(op);

@btime $A\$b;
x = zeros(num_dofs);
@btime gmres!($x, $A, $b);
@btime cg!($x, $A, $b);