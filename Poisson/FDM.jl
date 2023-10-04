using LinearAlgebra
# function poisson1d(T, a, b, f, u₀, τ, h)
# 	t = collect(0.0:τ:T)
# 	x = collect(a+h:h:b-h)
# 	n, m = length(t), length(x)
# 	ev = -2*ones(m)
# 	tv = ones(m-1)
# 	M = SymTridiagonal(ev, tv)/h^2
# 	U = zeros(Float64, m, n)  # julia 矩阵按列存储
# 	@. U[:, 1] = u₀(x)
# 	A = I-τ*M
# 	@inbounds for k in 2:n
# 		U[:, k] = A\(τ * f.(x, t[k]) + @view(U[:, k-1]))
# 	end
# 	return (t, x, U)
# end
function poisson1dgpt(T, a, b, f, u₀, τ, h)
	t = collect(0.0:τ:T)
	x = collect(a+h:h:b-h)
	n, m = length(t), length(x)
	ev = -2*ones(m)
	tv = ones(m-1)
	M = SymTridiagonal(ev, tv)/h^2
	U = zeros(Float64, m, n)  # julia 矩阵按列存储
	@. U[:, 1] = u₀(x)
	A = I-τ*M
    # fx_tk = zeros(m)
	tmp = zeros(m)
	@inbounds for k in 2:n
		fx_tk = f.(x, t[k])
		@. tmp = τ * fx_tk + @view(U[:, k-1])
		ldiv!(A, tmp)
		U[:, k] .= tmp
	end
	return (t, x, U)
end
T = 1
a = 0
b = 1
f(x, t) = (π^2-1/T)*exp(-t/T)*sin(π*x)
u₀(x) = sin(π*x)
h = 0.001
τ = 0.001