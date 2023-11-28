import Base:show, +, -, *, ==, zero, one, oneunit, rand

struct ℤ{N}
    buffer::Int
    N::Int
    value::Int
    function ℤ{N}(buffer::Int) where N
        value = buffer % N
        new{N}(buffer, N, value)
    end
end

show(io::IO, ::MIME"text/plain", n::ℤ{N}) where N = print(io, "(n.N))")
zero(::Type{ℤ{N}}) where N = ℤ{N}(0)
one(::Type{ℤ{N}}) where N = ℤ{N}(1)
zero(::ℤ{N}) where N = zero(ℤ{N})
oneunit(T::Type{ℤ{N}}) where N = one(T)
-(num::ℤ{N}) where N = ℤ{N}(N-num.value)
+(lhs::ℤ{N}, rhs::ℤ{N}) where N = ℤ{N}(lhs.value + rhs.value)
-(lhs::ℤ{N}, rhs::ℤ{N}) where N = lhs + (-rhs)
*(lhs::ℤ{N}, rhs::ℤ{N}) where N = ℤ{N}(lhs.value * rhs.value)
function rand(::Type{ℤ{N}}, dims::Tuple{Vararg{Int64, M}}) where {N, M}
    nums = [ℤ{N}(n) for n in 0:N-1]
    rand(nums, dims)
end