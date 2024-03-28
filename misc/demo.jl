#=
多
行
注
释
示
例
=#

# 单行注释

import Base:+, -, show

struct Point{N, T}
    data::NTuple{N, T}
end
Point(data...) = Point{length(data), eltype(data)}(data)

show(io::IO, p::Point) = print(io, "$(p.data)")

+(p₁::Point, p₂::Point) = Point(p₁.data .+ p₂.data)

Point(1, 2) + Point(2, 2)