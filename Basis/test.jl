x₀ = 1.0
x = Float64[x₀]
sizehint!(x, 1000)
current = x₀
next = current/2 + 2/current
push!(x, next) 
while abs(next - current) > 1e-5
    global current, next = next, current/2 + 2/current
    push!(x, next)
end