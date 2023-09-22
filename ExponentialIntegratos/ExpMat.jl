@doc raw"""
    ExpMatrix{T} where T <: Number
为了后续计算指数积分子与向量的乘积 $\varphi_n(A)v$ 的乘法运算符重载而创建的(参数化)数据类型, 
存储了要计算指数的矩阵 `basis` (即 A). 内置的构造函数要求矩阵 `basis` 必须为方阵.
这里 $\varphi_n$ ($n$ 对应 `oeder`) 为指数型 BDF$n$ 方法所对应的生成函数:
```math
\varphi_0(z) = \mathrm{e}^{z}, \quad \varphi_1(z) = \frac{\mathrm{e}^{z}-1}{z}, \quad \cdots. 
```
"""
struct ExpMatrix{T <: Number}
    basis :: Matrix{T}
    order :: Int
    function ExpMatrix{T}(mat::Matrix{T}; order::Int) where T <: Number
        @assert size(mat, 1) == size(mat, 2) "必须为方阵"
        @assert order >= 0 "必须为非负整数"
        new{T}(mat, order)
    end
end
"""外置构造函数, 将类型参数隐去."""
ExpMatrix(mat::Matrix{T}; order::Int = 0) where T <: Number = ExpMatrix{T}(mat, order=order)
"""构造 \$\varphi_n(tA)\$ 型指数积分子 """
ExpMatrix(t::Number, mat::Matrix; order::Int = 0) = ExpMatrix(t .* mat; order=order)
import Base.*
"""
    *(expmat::ExpMatrix, b::AbstractVector; method::String="Contour-Krylov")
多重派发实现乘法运算符 `*` 的重载, 利用 `method` 方法计算 \$\varphi_n(tA)b\$.
# 参数
- `expmat` : 未计算的 \$\varphi_n(tA)\$ 或者 \$\varphi_n(A)\$
- `b` : 向量
- `method` : 实现算法, 默认为 `Contour-Krylov`, 具体含义为
 - `Krylov` : 利用 Krylov 子空间迭代计算乘积
 - `Contour` : 利用周线数值积分方法计算乘积
 - `Coutour-Krylov` : 结合周线数值积分方法和 Krylov 子空间迭代计算乘积
"""
function *(expmat::ExpMatrix, b::AbstractVector; method::String="Contour-Krylov")
    dim = size(expmat.basis, 1)
    @assert dim == length(b) "矩阵和向量的维度不一致"
    methods = ["Krylov", "Contour", "Contour-Krylov"]
    @assert method in methods "输入的可选关键字参数有误, 请选择 `Krylov`, `Contour` 或者 `Contour-Krylov`"
    A = expmat.basis
    return 1
end

A = ExpMatrix(1, rand(3, 3))