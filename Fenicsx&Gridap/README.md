# FEniCSx 与 Gridap.jl 的用法区别

## 网格剖分
以求解区域 $[0, 1]\times [0, 1]$ 为例.

### FEniCSx
`MPI` 多进程剖分网格, 通过命令 `mpi` 的参数 `-n` 指定所使用的核心数:  `mpi -n 2 python xxx.py`. 

```python
import numpy as np
from mpi4py import MPI
from dolfinx import mesh

nx, ny = 8, 8  # 剖分单元数量
model = mesh.create_unit_square(MPI.COMM_WORLD, nx, ny, mesh.CellType.quadrilateral)
```
`mesh.CellType.quadrilateral` 表示矩形单元, 三角单元可用 `mesh.CellType.triangle`. 
此外, 创建矩形求解区域 $[a, b]\times [c, d]$ 可使用 `mesh.create_rectangle`:
```python
mesh.create_rectangle(MPI.COMM_WORLD, [np.array([a, c]), np.array([b, d])], [nx, ny], mesh.CellType.triangle)
```

导入 `Gmsh` 的网格请看 [Interfacing with GMSH in DOLFINx](https://jsdokken.com/dolfinx-tutorial/chapter1/membrane_code.html).

### Gridap.jl
`Gridap.jl` 可以使用 `CartesianDiscreteModel` 创建矩形区域的笛卡尔网格.
```julia
using Gridap
nx, ny = 8, 8
domain = (0, 1, 0, 1)
partition = (nx, ny)
model = CartesianDiscreteModel(domain, partition)
```
还可使用 `GridapGmsh` 导入 `Gmsh` 生成的 `.msh` 文件.
```julia
using GridapGmsh
model = GmshDiscreteModel("xxx.msh")
```

## 定义有限元函数空间
以 Dirichlet 边界条件的 Poisson 方程为例, 边界条件
$$
u\\_d(x, y) = 1 + \sin{x} + \cos{y}.
$$

### FEniCSx
```python
from dolfinx import fem
from dolfinx.fem import FunctionSpace
# P^{degree} 有限元空间
degree = 1
V = FunctionSpace(model, ("Lagrange", degree))  
# 指定边界条件
tdim = model.topology.dim
fdim = tdim - 1
model.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(model.topology)
boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
ud = fem.Function(V)
ud.interpolate(lambda x: 1 + np.sin(x[0]) + np.cos(x[1]))
bc = fem.dirichletbc(uD, boundary_dofs)
```

### Gridap.jl
```julia
degree = 1
# 参考单元
reffe = ReferenceFE(lagrangian, Float64, degree)
# 测试函数空间
V = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags="boundary")
# 试用函数空间
U = TrialFESpace(V, x->1+sin(x[1])+cos(x[2]))
```

## 定义双线性型与右端泛函
源项 
$$
f(x, y) = \sin{x} + \cos{y}.
$$

### FEniCSx
```python
import ufl
from ulf import dot, dx, grad
# 源项
f = fem.Function(V)
f.interpolate(lambda x: np.sin(x[0]) + np.cos(x[1]))
# 弱形式
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
a = dot(grad(u), grad(v)) * dx
L = f * v * dx
```

### Gridap.jl
```julia
# 数值积分
order = degree + 1
Ω = Triangulation(model)
dΩ = Measure(Ω, order)
# 源项
f(x) = sin(x[1]) + cos(x[2])
# 弱形式
a(u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
L(v) = ∫(f * v)dΩ
```

## 组装线性系统与求解

### FEniCSx
```python
from dolfinx.fem.petsc import LinearProblem
problem = LinearProblem(a, L, bcs=[bc])
uh = problem.solve()
```
手动组装
```python
bilinear_form = fem.form(a)
linear_form = fem.form(L)
# 组装矩形/向量
A = assemble_matrix(bilinear_form, bcs=[bc])
A.assemble()
b = assemble_vector(linear_form, bcs=[bc])
b.assemble()
```

### Gridap.jl
```julia
op = AffineFEOperator(a, L, U, V)
uh = solve(op)
```
手动组装
```julia
A = assemble_matrix(a, U, V)
b = assemble_vector(L, V)
```