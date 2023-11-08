abstract type NDArray end
struct MyVector <: NDArray name::String end
struct MyMatrix <: NDArray name::String end

⋅(a::NDArray, b::NDArray) = "不确定"
⋅(a::MyVector, b::MyVector) = "向量点乘"
⋅(a::MyMatrix, b::MyVector) = "矩阵乘向量"
⋅(a::MyVector, b::MyMatrix) = "维数不一致"
⋅(a::MyMatrix, b::MyMatrix) = "矩阵乘积"

result(a::NDArray, b::NDArray) =  println("$(a.name)和$(b.name)相乘的结果为$(a ⋅ b)")

v = MyVector("向量")
m = MyMatrix("矩阵")

result(v, v); result(m, v); result(v, m); result(m, m)
