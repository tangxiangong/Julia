struct Fraction <: Real
    num::Int    # 分子
    den::Int    # 分母
    function Fraction(num::Int, den::Int)
        @assert !iszero(den)
        iszero(num) && return new(0, 1)
        gcd_value = gcd(num, den)
        if sign(num) == sign(den) == -1
            new(abs(num) ÷ gcd_value , abs(den) ÷ gcd_value)
        else
            new(num ÷ gcd_value , den ÷ gcd_value)
        end
    end
end

Base.iszero(f::Fraction) = iszero(f.num)

Base.zero(::Type{Fraction}) = Fraction(0, 1)

Base.similar(::Fraction) = zero(Fraction)

isnegative(f::Fraction) = sign(f.den * f.num) == -1

ispositive(f::Fraction) = sign(f.den * f.num) == 1

function Base.show(io::IO, ::MIME"plain/txt", f::Fraction)
    if iszero(f)
        str = "0" 
    elseif ispositive(f)
        str = "$(f.num)/$(f.den)"
    else
        str = "-$(abs(f.num))/$(abs(f.den))"
    end
    print(io, str)
end

Base.:(==)(lhs::Fraction, rhs::Fraction) = lhs.den == rhs.den && lhs.num == rhs.num

Base.:-(f::Fraction) = Fraction(-f.num, f.den)

function Base.:+(lhs::Fraction, rhs::Fraction)
    iszero(lhs) && return Fraction(rhs.num, rhs.den)
    iszero(rhs) && return Fraction(lhs.num, lhs.den)
    num = lhs.num * rhs.den + rhs.num * lhs.den
    den = lhs.den * rhs.den
    Fraction(num, den)
end

Base.:+(lhs::Fraction, rhs::Int) = nothing

Base.:-(lhs::Real, rhs::Real) = lhs + (-rhs)
