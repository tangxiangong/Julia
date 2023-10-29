struct Point{T, N}
    ndim::Int
    data::Vector{T}
end
Point(data::Vector) = Point{eltype(data), length(data)}(length(data), data)