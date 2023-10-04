function strassen(A::Matrix{<:Number}, B::Matrix{<:Number})
    size(A, 2) == size(B, 1) || throw(DomainError("矩阵维数不一致!"))
    m, l, n = size(A, 1), size(A, 2), size(B, 2)
    C = Matrix{promote_type(eltype(A), eltype(B))}(undef, m, n)
    
    return C
end