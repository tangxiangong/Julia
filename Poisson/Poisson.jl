using Gridap, GridapGmsh
import GridapGmsh:gmsh

"""
    generate_mesh(height::Float64=1.0, width::Float64=1.0, lc::Float64=0.05, 
                dirichlet_tags::String="boundary", path::String="./", name::String="model") :: String 
调用 Gmsh API 生成第一象限中原点为左下顶点、高度为 `height`、宽度为 `width` 的矩形区域的三角网格, 单元尺寸为 `lc`, 
Dirichlet 边界的标签为 `dirichlet_tags`, 并且将网格数据保存在路径 `path` 下的名为 `name` 且后缀名为 `.msh` 的文件中.
"""
function generate_mesh(height::Float64=1.0, width::Float64=1.0, lc=0.01, 
    dirichlet_tags::String="boundary", path::String="./", name::String="model") :: String 
    # 初始化
    gmsh.initialize()
    # 新建名为 name 的模型
    gmsh.model.add(name)
    # 添加节点
    p₁ = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc)
    p₂ = gmsh.model.geo.addPoint(width, 0.0, 0.0, lc)
    p₃ = gmsh.model.geo.addPoint(width , height, 0.0, lc)
    p₄ = gmsh.model.geo.addPoint(0.0, height, 0.0, lc)
    # 按逆时针将节点连接成线
    l₁ = gmsh.model.geo.addLine(p₁, p₂)
    l₂ = gmsh.model.geo.addLine(p₂, p₃)
    l₃ = gmsh.model.geo.addLine(p₃, p₄)
    l₄ = gmsh.model.geo.addLine(p₄, p₁)
    # 将线按逆时针组合成闭合曲线形成一个面
    curve = gmsh.model.geo.addCurveLoop([l₁, l₂, l₃, l₄])
    surface = gmsh.model.geo.addPlaneSurface([curve])
    # 同步化, 生成数据结构来存储模型中点线面的信息
    gmsh.model.geo.synchronize()
    # 指定物理边界
    dirichlet = gmsh.model.addPhysicalGroup(1, [l₁, l₂, l₃, l₄])
    gmsh.model.addPhysicalGroup(2, [surface], 1)
    gmsh.model.setPhysicalName(1, dirichlet, dirichlet_tags)
    # 生成(2d)网格
    gmsh.model.mesh.generate(2)
    # 存储网格
    file_path = string(path, name, ".msh")
    gmsh.write(file_path)
    gmsh.finalize()
    return file_path
end

"""
    PoissonSolver(mesh::String, f::Function, dirichlet_tags::String="boundary")
网格剖分为 `mesh`, 右端项为 `f` 且 Dirichlet 边界为 `dirichlet_tags` 的齐次 Dirichlet 边界的 Poisson 方程的求解器.
"""
function PoissonSolver(mesh::String, f::Function, dirichlet_tags::String="boundary") :: FEFunction
    # 导入网格
    model = GmshDiscreteModel(mesh)
    # 基函数的次数
    order = 1
    # 参考单元上的基函数类型
    reffe = ReferenceFE(lagrangian, Float64, order)
    # 测试函数空间
    V = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags=dirichlet_tags)
    # 试用函数空间
    U = TrialFESpace(V, x->0)
    # 数值积分的精度
    degree = 2
    # 将离散模型三角化、可测化
    Ω = Triangulation(model)
    # 积分单元的测度
    dΩ = Measure(Ω, degree)
    # 双线性型
    a(u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ 
    # 右端项
    b(v) = ∫(v*f) * dΩ  
    # 组装得到线性系统
    op = AffineFEOperator(a,b,U,V)
    # 求解器
    solver = LinearFESolver(LUSolver())
    # 数值解
    uₕ = solve(solver, op)
    writevtk(Ω, "result", cellfields=["solution"=>uₕ])
    return uₕ
end

f(x) = 2 * π^2 * sin(π*x[1]) * sin(π*x[2])
mesh = generate_mesh()
uₕ = PoissonSolver(mesh, f)
# u(x) = sin(π*x[1])*sin(π*x[2])
# e = u - uₕ 
# el2 = sqrt(sum(∫(e * e) * dΩ))
# println("L2 error = $el2")
