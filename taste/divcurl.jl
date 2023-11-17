import Base: show


function print_aux(a::NTuple{4, <:Real})
    str = ""

    a[1] != 0 && abs(a[1]) != 1 && (str *= "$(a[1])x")
    a[1] == 1 && (str *= "x")
    a[1] == -1 && (str *= "-x")
    
    a[2] > 0 && a[2] != 1 && (str *= " + $(a[2])y")
    a[2] == 1 && (str *= " + y")
    a[2] < 0 && a[2] != -1 && (str *= " - $(abs(a[2]))y")
    a[2] == -1 && (str *= " - y")

    a[3] > 0 && a[3] != 1 && (str *= " + $(a[3])z")
    a[3] == 1 && (str *= " + z")
    a[3] < 0 && a[3] != -1 && (str *= " - $(abs(a[3]))z")
    a[3] == -1 && (str *= " - z")

    isempty(str) && return "$(a[4])"
    a[4] > 0 && (str *= " + $(a[4])")
    a[4] < 0 && (str *= " - $(abs(a[4]))")

    str
end

"""
    三元线性标量函数

```math
    p(x, y, z) = a_1 x + a_2 y + a_3 z + a_4
```
"""
struct ScalarFunction
    a::NTuple{4, <:Real}
end

(p::ScalarFunction)(x::Real, y::Real, z::Real) = p.a[1] * x + p.a[2] * y + p.a[3] * z + p.a[4]

function show(io::IO, p::ScalarFunction)
    str = print_aux(p.a)
    print(io, str)
end


"""
    三元线性向量函数

```math
p(x, y, z) = (a_1 x + a_2 y + a_3 z + a_4, b_1 x + b_2 y + b_3 z + b_4, c_1 x + c_2 y + c_3 z + c_4)
```
"""
struct VectorFunction
    a::NTuple{4, <:Real}
    b::NTuple{4, <:Real}
    c::NTuple{4, <:Real}
end

function (p::VectorFunction)(x::Real, y::Real, z::Real)
    p₁ = p.a[1] * x + p.a[2] * y + p.a[3] * z + p.a[4]
    p₂ = p.b[1] * x + p.b[2] * y + p.b[3] * z + p.b[4]
    p₃ = p.c[1] * x + p.c[2] * y + p.c[3] * z + p.c[4]
    [p₁ p₂ p₃]
end

function show(io::IO, p::VectorFunction)
    s₁ = print_aux(p.a)
    s₂ = print_aux(p.b)
    s₃ = print_aux(p.c)
    str = "($(s₁), $(s₂), $(s₃))"
    print(io, str)
end 

function ∇(p::ScalarFunction)
    T = eltype(p.a)
    O = zero(T)
    a = (O, O, O, p.a[1])
    b = (O, O, O, p.a[2])
    c = (O, O, O, p.a[3])
    VectorFunction(a, b, c)
end

function ⋅(::typeof(∇), p::VectorFunction)
    T = promote_type(eltype(p.a), eltype(p.b), eltype(p.c))
    div = p.a[1] + p.b[2] + p.c[3]
    O = zero(T)
    ScalarFunction((O, O, O, div))
end

function ×(::typeof(∇), p::VectorFunction)
    T = promote_type(eltype(p.a), eltype(p.b), eltype(p.c)) 
    O = zero(T)
    a = (O, O, O, p.c[2] - p.b[3])
    b = (O, O, O, p.a[3] - p.c[1])
    c = (O, O, O, p.b[1] - p.a[2])
    VectorFunction(a, b, c)
end
p = VectorFunction((1, 2, 3, 4), (4, 3, 2, 1), (5, 6, 7, 8))
∇ ⋅ p
∇ × p