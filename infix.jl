±(a, b) = (a+b, a-b)
∓(a, b) = (a-b, a+b)
# issubset == ⊆
⊂(A::AbstractSet, B::AbstractSet) = issubset(A, B)
# Kronecker product
⊗ = kron

1 ± 2
1 ∓ 2

Set([1, 2]) ⊂ Set([1, 2, 3])
Set([1, 2]) ⊂ Set([1])

A = [1 0; 0 1]
B = [1 1; 1 1]
A ⊗ B