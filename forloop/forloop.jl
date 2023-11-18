using BenchmarkTools

function seriel(f::Function, x::Vector)
    y = similar(x)
    @inbounds for k in eachindex(x)
        y[k] = f(x[k])
    end
    y
end

function seriel!(y::Vector, f::Function, x::Vector)::Nothing
    @inbounds for k in eachindex(x)
        y[k] = f(x[k])
    end
    nothing
end


n = 100_000
x = rand(n)
y = Vector{Float64}(undef, n)
println("Julia for 循环")
@btime seriel(sin, $x);
println("Julia inplace for 循环")
@btime seriel!($y, sin, $x);
println("Julia 向量化")
@btime sin.($x);