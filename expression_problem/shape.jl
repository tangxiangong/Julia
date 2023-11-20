#=
假设已有一个圆类和一个矩形类，都实现了计算其周长的方法。
如何在不修改原有代码的前提条件下，新增一个三角形类、实现其计算周长的算法，并且对这三个类型新增计算面积的方法。
=#

# 已有代码
abstract type AbstractShape end
perimeter(::AbstractShape) = println("算法未实现")

struct Circle{T<:Real} <: AbstractShape
    radius::T
    function Circle(radius::T) where T <: Real
        @assert radius > 0
        new{T}(radius)
    end
end
perimeter(c::Circle) = 2 * π * c.radius

struct Rectangle{T₁<:Real, T₂<:Real} <: AbstractShape
    height::T₁
    width::T₂
    function Rectangle(height::T₁, width::T₂) where {T₁ <: Real, T₂ <: Real}
        @assert height > 0 && width > 0
        new{T₁, T₂}(height, width)
    end
end
perimeter(r::Rectangle) = 2(r.height + r.width)

# 新增代码
struct Triangle{T₁<:Real, T₂<:Real, T₃<:Real} <: AbstractShape
    a::T₁
    b::T₂
    c::T₃
    function Triangle(a::T₁, b::T₂, c::T₃) where {T₁ <: Real, T₂ <: Real, T₃ <: Real}
        @assert a > 0 && b > 0 && c > 0
        @assert a + b > c && a + c > b && b + c > a
        new{T₁, T₂, T₃}(a, b, c)
    end
end
perimeter(t::Triangle) = t.a + t.b + t.c

area(::AbstractShape) = println("算法未实现")

area(c::Circle) = π * c.radius^2

area(r::Rectangle) = r.height * r.width

function area(t::Triangle)
    s = perimeter(t)/2
    sqrt(s * (s - t.a) * (s - t.b) * (s - t.c))
end

