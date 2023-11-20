function bm(length_t, τ=0.01)
    t = collect(0.0:τ:length_t)
    n = length(t)
    noise = zeros(Float64, n-1)
    x = zeros(Float64, n)
    @inbounds for i in 2:n
        noise[i-1] = sqrt(τ) * randn()
        x[i] = x[i-1] + noise[i-1]
    end
    return t, x
end

function msd(Ts, N, τ=0.01) 
	T_max = Ts[end]
	index = @. Int(Ts / τ) + 1
	s = zeros(Float64, length(Ts))
	Threads.@threads for _ in 1:N
		~, x = bm(T_max)
		@. s += x[index]^2
	end
	s ./= N
	return s
end
using LoopVectorization
function mynorm(a::Vector{Float64}, b::Vector{Float64})   
    @assert length(a) == length(b)
    s = 0.0
    @turbo for i in eachindex(a)
        s += (a[i] - b[i])^2
    end
    s = sqrt(s)
    return s
end

Ts = [10, 20, 30, 50, 100]
N = 1000
@btime msd(Ts, N)