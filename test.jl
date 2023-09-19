function f(x::Vector{Float64})
    a = Float64[]
    for i in x
        push!(a, i)
    end
    return nothing 
end

function g(x::Vector{Float64})
    a = Float64[]
    sizehint!(a, 1000000)
    for i in x
        push!(a, i)
    end
    return nothing
end
