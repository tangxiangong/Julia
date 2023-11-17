∇(f::Function) = f

⋅(::typeof(∇), f::Function) = x -> 2 * f(x)

×(::typeof(∇), f::Function) = x -> 3 * f(x)

f = sinpi

Δ(f::Function) = ∇ ⋅ (∇(f))

f(1/6)

(∇ ⋅ f)(1/6)

(∇ × f)(1/6)

(Δ(f))(1) == 2 * (∇(f))(1) == 2 * f(1)