struct Point{N, T}
    v::NTuple{N, T}
    function Point(v...)
        if isempty(v)
            new{0, Int}(v)
        else
            new{length(v), typeof(v)}(v)
        end
    end
end

Base.isempty(p::Point) = isempty(p.v)
function Base.show(io::IO, p::Point)
    if isempty(p)
        print(io, "Empty Point.")
    else
        print(io, "Point{$(p.N), $(p.T)}$(p.v)")
    end
end


p = Point(1, 2, 3)

function add(x, y, z...)
    isempty(z) && return x + y
    s = x + y
    for k in z
        s += k
    end
    return s
end
