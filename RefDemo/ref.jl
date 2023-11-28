function swap!(a::Ref{T}, b::Ref{T}) where T
    c::T = a[]
    a[] = b[]
    b[] = c
    nothing
end

a = Ref(1.2)
b = Ref(2.3)
swap!(a, b)
a[]
b[]

function square(x::AbstractVector, n::Int)
    s = 0.0
    for k in x
        s += 1/k^n
    end
    s
end
square(x::Ref, n::Int) = square(x[], n)

x = collect(1:10)
# 不能直接广播
square.(x, [1, 2, 3])
# 使用 Ref 可以广播
square.(Ref(x), [1, 2, 3]) == map(n->square(x, n), [1, 2, 3])