function bisection(f::Function, a::Float64, b::Float64, tol::Float64=1e-16)
    if f(a) * f(b) > 0
        error("Function must have opposite signs at the endpoints")
    end
    c = (a + b) / 2
    while abs(f(c)) > tol
        if f(a) * f(c) < 0
            b = c
        else
            a = c
        end
        c = (a + b) / 2
    end
    return c
end
