using OrdinaryDiffEq, LinearAlgebra, SparseArrays, BenchmarkTools

const N = 32
const xyd_brusselator = range(0,stop=1,length=N)
brusselator_f(x, y, t) = (((x-0.3)^2 + (y-0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.
limit(a, N) = a == N+1 ? 1 : a == 0 ? N : a
function brusselator_2d_loop(du, u, p, t)
  A, B, alpha, dx = p
  alpha = alpha/dx^2
  @inbounds for I in CartesianIndices((N, N))
    i, j = Tuple(I)
    x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
    ip1, im1, jp1, jm1 = limit(i+1, N), limit(i-1, N), limit(j+1, N), limit(j-1, N)
    du[i,j,1] = alpha*(u[im1,j,1] + u[ip1,j,1] + u[i,jp1,1] + u[i,jm1,1] - 4u[i,j,1]) +
                B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1] + brusselator_f(x, y, t)
    du[i,j,2] = alpha*(u[im1,j,2] + u[ip1,j,2] + u[i,jp1,2] + u[i,jm1,2] - 4u[i,j,2]) +
                A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]
    end
end
p = (3.4, 1., 10., step(xyd_brusselator))

function init_brusselator_2d(xyd)
  N = length(xyd)
  u = zeros(N, N, 2)
  for I in CartesianIndices((N, N))
    x = xyd[I[1]]
    y = xyd[I[2]]
    u[I,1] = 22*(y*(1-y))^(3/2)
    u[I,2] = 27*(x*(1-x))^(3/2)
  end
  u
end
u0 = init_brusselator_2d(xyd_brusselator)
prob_ode_brusselator_2d = ODEProblem(brusselator_2d_loop,u0,(0.,11.5),p)

du = similar(u0)
brusselator_2d_loop(du, u0, p, 0.0)
du[34] # 802.9807693762164
du[1058] # 985.3120721709204
du[2000] # -403.5817880634729
du[end] # 1431.1460373522068
du[521] # -323.1677459142322

du2 = similar(u0)
brusselator_2d_loop(du2, u0, p, 1.3)
du2[34] # 802.9807693762164
du2[1058] # 985.3120721709204
du2[2000] # -403.5817880634729
du2[end] # 1431.1460373522068
du2[521] # -318.1677459142322

using Symbolics
du0 = copy(u0)
jac_sparsity = Symbolics.jacobian_sparsity((du,u)->brusselator_2d_loop(du,u,p,0.0),du0,u0)

f = ODEFunction(brusselator_2d_loop;jac_prototype=float.(jac_sparsity))
prob_ode_brusselator_2d_sparse = ODEProblem(f,u0,(0.,11.5),p,tstops=[1.1])

using IncompleteLU
function incompletelu(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    Pl = ilu(convert(AbstractMatrix,W), τ = 150.0)
  else
    Pl = Plprev
  end
  Pl,nothing
end

# Required due to a bug in Krylov.jl: https://github.com/JuliaSmoothOptimizers/Krylov.jl/pull/477
Base.eltype(::IncompleteLU.ILUFactorization{Tv,Ti}) where {Tv,Ti} = Tv

sol = solve(prob_ode_brusselator_2d_sparse,KenCarp47(linsolve=KrylovJL_GMRES(),precs=incompletelu,concrete_jac=true),abstol=1e-12,reltol=1e-12)
sol.t[104] # 1.1
sol[104][1] # 0.5305455584661771
sol[end][1] # 3.2723204220222195

@btime solve(prob_ode_brusselator_2d_sparse,KenCarp47(linsolve=KrylovJL_GMRES(),precs=incompletelu,concrete_jac=true),save_everystep=false,abstol=1e-8,reltol=1e-8);
# 1.422 s (304441 allocations: 140.77 MiB)

using DifferentialEquations
@btime solve(prob_ode_brusselator_2d,alg_hints=[:stiff],save_everystep=false,abstol=1e-8,reltol=1e-8);
# 2.076 s (551306 allocations: 23.87 MiB)

using Sundials, ModelingToolkit
prob_ode_brusselator_2d_mtk = ODEProblem(modelingtoolkitize(prob_ode_brusselator_2d_sparse),[],(0.0,11.5),jac=true,sparse=true);

using LinearAlgebra
u0 = prob_ode_brusselator_2d_mtk.u0
p  = prob_ode_brusselator_2d_mtk.p
const jaccache = prob_ode_brusselator_2d_mtk.f.jac(u0,p,0.0)
const WW = I - 1.0*jaccache

prectmp = ilu(WW, τ = 50.0)
const preccache = Ref(prectmp)

function psetupilu(p, t, u, du, jok, jcurPtr, gamma)
  if jok
    prob_ode_brusselator_2d_mtk.f.jac(jaccache,u,p,t)
    jcurPtr[] = true

    # W = I - gamma*J
    @. WW = -gamma*jaccache
    idxs = diagind(WW)
    @. @view(WW[idxs]) = @view(WW[idxs]) + 1

    # Build preconditioner on W
    preccache[] = ilu(WW, τ = 5.0)
  end
end

function precilu(z,r,p,t,y,fy,gamma,delta,lr)
  ldiv!(z,preccache[],r)
end

@btime solve(prob_ode_brusselator_2d_sparse, CVODE_BDF(linear_solver=:GMRES,prec=precilu,psetup=psetupilu,prec_side=1),save_everystep=false);
# 72.949 ms (15412 allocations: 55.75 MiB)
