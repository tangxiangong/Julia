### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 209f1308-7ee8-11ee-37a4-efac0b62e389
using LinearAlgebra, MKL, Gridap

# ╔═╡ 32bfbf9c-fa5b-4f99-ad83-c612c9b0efe2
md"""
# Demo
## 模型 && 数值格式
```math
\begin{cases}
\partial_t \vec{v} + (\vec{v}\cdot∇)\vec{v} - \Delta\vec{v} + \nabla p = \vec{f} + \vec{j}u, & \\
\nabla\cdot\vec{v} = 0, & \\
\partial_t u + \vec{v}\cdot\nabla u - \Delta u = g,
\end{cases}
```
求解区域 ``(x,t)\in\Omega\times (0,T]``. 
记向前差商算子为 ``\bar{\partial}_{\tau}``, 则欧拉--有限元格式为 
```math
\begin{cases}
(\bar{\partial}_{\tau}\vec{v}_h^n, \vec{w}_h) + ((\vec{v}_h^{n-1}\cdot\nabla)\vec{v}_h^n, \vec{w}_h)+(\nabla\vec{v}_h^n,\nabla\vec{w}_h) - (p,\nabla\cdot\vec{w}_h) - (\vec{j}u_h^n,\vec{w}_h) = (\vec{f}^n, \vec{w}_h), & \\
(\nabla\cdot\vec{v}_h^n, q_h) = 0, & \\
(\bar{\partial}_{\tau}u_h^n,\varphi_h) + (\vec{v}_h^n\cdot\nabla u_h^n, \varphi_h) + (\nabla u_h^n, \nabla\varphi_h) = (g^n, \varphi_h).
\end{cases}
```
"""

# ╔═╡ 81d7ad67-3362-4154-be4e-25f1a3faa35a
md"""
质量矩阵对应的双线性型为
```math
\text{mass}((\vec{v},p,u),(\vec{w},q,\varphi)) = (\vec{v}, \vec{w}) + (u, \varphi).
```
"""

# ╔═╡ cb8e8bc1-c2c5-4c62-8124-1d0fc520f181
md"""
刚度矩阵对应的双线性型为
```math
\begin{aligned}
\text{stiff}((\vec{v}, p, u), (\vec{w}, q, φ)) &= (\nabla\vec{v},\nabla\vec{w})-(p,\nabla\cdot\vec{w})-(\vec{j}u,\vec{w}) \\
&+(\nabla\cdot\vec{v}, q) + (\vec{v}\cdot\nabla u, \varphi) + (\nabla u,\nabla\varphi).
\end{aligned}
```
"""

# ╔═╡ a177b3b1-e4a3-446b-94b5-aff5634e6571
md"""
对流项对应的双线性型为 
```math
\text{conv}((\vec{v},p,u),(\vec{w},q,\varphi)) = ((\vec{a}\cdot\nabla)\vec{v},\vec{w}).
```
"""

# ╔═╡ ec129a2a-0072-4a4a-b83c-21bda90176b3
md"## 程序实现"

# ╔═╡ 1e6c7b24-524e-423c-9992-1d028abd0309
md"""
### 1. 网格剖分
设 ``\Omega=(0, 1)^2``, ``T=1``, 剖分尺寸 ``h=1/128``, ``\tau=1/32``.
"""

# ╔═╡ c635e775-77d2-4f4a-8472-29aa7b4acb81
domain = (0, 1, 0, 1)    # 正方形求解区域

# ╔═╡ 034cdf88-ca69-4938-96b4-904487530cc2
n = 128     # 垂直/水平方向的单元数

# ╔═╡ d8c681df-f4b8-482a-b1d5-49b0f9f02be8
model = CartesianDiscreteModel(domain, (n, n))    # 网格剖分

# ╔═╡ d9b279ea-1f36-45d4-9661-3cd87582426e
md"""
### 2. 有限元空间
"""

# ╔═╡ 64f868d6-e5c5-4b82-83a8-901c83ee1371
order = 2

# ╔═╡ 34dbe9ce-a135-4880-8fdc-37e98c784f27
begin
	reffeᵥ = ReferenceFE(lagrangian, VectorValue{2, Float64}, order)  # P2
	reffeₚ = ReferenceFE(lagrangian, Float64, order-1; space=:P)      # P1
	reffeᵤ = ReferenceFE(lagrangian, Float64, 1)   # 线性元
	# Test function space
	W = TestFESpace(model, reffeᵥ, conformity=:H1, dirichlet_tags="boundary")
	Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean)
	Z = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags="boundary")
	# Trial function space
	V = TrialFESpace(W, VectorValue(0.0, 0.0))
	P = TrialFESpace(Q)
	U = TrialFESpace(Z, 0.0)
	X = MultiFieldFESpace([V, P, U])
	Y = MultiFieldFESpace([W, Q, Z])
end

# ╔═╡ ca063757-098a-4051-9f7b-75b4ca136fca
num_dofs = num_free_dofs(X)   # 总自由度

# ╔═╡ a73ec73e-be63-450c-bb0f-ec43cca65621
md"### 3. 双线性型"

# ╔═╡ 6dfd5b56-3162-4e9f-82fd-b53c88ec9a48
degree = order + 1; Ω = Triangulation(model); dΩ = Measure(Ω, degree)

# ╔═╡ 2b4de08b-10a9-4d1b-8c92-b515a91f7702
mass((v, p, u), (w,q,φ)) = ∫(v⋅w + u*φ)dΩ

# ╔═╡ ded709b0-f743-4ba0-8eb7-9ecbbe28b2c6
M = assemble_matrix(mass, X, Y)    # 质量矩阵

# ╔═╡ fc85ccd9-0969-4688-80e1-3b9d072a318b
md"""
### 4. 右端项与求解
"""

# ╔═╡ 6f78b279-879c-46ae-a7ca-ba06f3a847b0
md"""
初值 ``\vec{v}(0) = (-\pi\sin(\pi x)\cos(\pi y), \pi\cos(\pi x)\sin(\pi y))``, ``u(0)=\sin(\pi x)\sin(\pi y)``.

右端项 ``\vec{f}(x, y,t) = (\exp(-t)x, \exp(-t)y)``, ``g(x,y,t)=\exp(-t)*x*y``.

参数 ``\vec{j}=(1, 1)``.
"""

# ╔═╡ 5029fc67-4660-4e38-b4ba-87cffe858e75
begin
	v₀(x) = VectorValue(-π*sin(π*x[1])*cos(π*x[2]), π*cos(π*x[1])*sin(π*x[2]))
	u₀(x) = sin(π*x[1]) * sin(π * x[2])
	f(x, t) = VectorValue(exp(-t) * x[1], exp(-t)*x[2])
	g(x, t) = exp(-t)*x[1]*x[2]
	j⃗ = VectorValue(1., 1.)
	nothing
end

# ╔═╡ bafefb66-34a5-48e8-914c-9391953a310d
stiff((v, p, u), (w, q, φ)) = ∫(∇(v)⊙∇(w) - (∇⋅w)*p - (j⃗⋅w)*u + (∇⋅v)*q + (v⋅∇(u))*φ + ∇(u)⋅∇(φ))dΩ

# ╔═╡ b19e69cc-0fd2-4afb-8284-f9a18e2fc765
S = assemble_matrix(stiff, X, Y)

# ╔═╡ 9a3844f0-f958-4732-93c6-be434c19149f
function Vector{FEFunction}(n::Int, U::TrialFESpace)
	num_dofs = num_free_dofs(U)
	fefunction = FEFunction(U, zeros(num_dofs))
	[fefunction for _ in 1:n]
end

# ╔═╡ c56564a6-8d92-4b73-825f-a245d7950592
# 时间网格剖分
function temporalmesh(T, τ, α)
	@assert 1/2 < α < 1 && T > 0 && τ > 0
	t = Vector{Float64}()
	sizehint!(t, round(Int, T/τ))
	append!(t, 0.0)
	append!(t, T*(τ/T)^(1/(1-α)))
	last_time = t[end]
	while last_time < T
		τn = τ * (last_time/T)^α
		last_time += τn
		append!(t, last_time)
	end
	t[end] = T
	t
end

# ╔═╡ f9c935e6-aae7-4eb3-90fd-2657b1a4f86e
T = 1.0; α = 0.6; τ = 1/32; t = temporalmesh(T, τ, α); length(t)

# ╔═╡ ce06b4cd-28b8-4617-96b4-0d0eb0e4794b
begin 
	τs = diff(t)
	previous_v = v₀
	previous_u = u₀
	A = similar(M)
	temp = similar(M)
	b = zeros(num_dofs)
	solution = zeros(num_dofs)
	vs = Vector{FEFunction}(length(t), V)
	us = Vector{FEFunction}(length(t), U)
	vs[1] = interpolate(v₀, V)
	us[1] = interpolate(u₀, U)
	@inbounds for k in eachindex(τs)
		stiff_conv((v, p, u), (w, q, φ)) = ∫(w⊙((∇(v))'⋅vs[k]))dΩ
		temp .= assemble_matrix(stiff_conv, X, Y)
		A .= M + τs[k] * (S + temp)
		fn(x) = τs[k]*f(x, t[k+1])
		gn(x) = τs[k]*g(x, t[k+1])
		ℓ((w, q, φ)) = ∫(fn⋅w + gn*φ + vs[k]⋅w + us[k]*φ)dΩ
		b .= assemble_vector(ℓ, Y)
		solution .= A\b
		# global previous_v, ~, previous_u = FEFunction(X, solution)
		vs[k+1], ~, us[k+1] = FEFunction(X, solution)
	end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Gridap = "56d4f2e9-7ea1-5844-9cf6-b9c51ca7ce8e"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MKL = "33e6dc65-8f57-5167-99aa-e5a354878fb2"

[compat]
Gridap = "~0.17.20"
MKL = "~0.6.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0-rc1"
manifest_format = "2.0"
project_hash = "89402dc3b71e9067b7b8ca9cfc0026a2c5c3a578"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "02f731463748db57cc2ebfbd9fbc9ce8280d3433"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "16267cf279190ca7c1b30d020758ced95db89cd0"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.5.1"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "af43df5704827c8618afd36eb56fcab20d3041ee"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.4.3"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AutoHashEquals]]
deps = ["Pkg"]
git-tree-sha1 = "daaeb6f7f77b88c072a83a2451801818acb5c63b"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "2.1.0"

[[deps.BSON]]
git-tree-sha1 = "2208958832d6e1b59e49f53697483a84ca8d664e"
uuid = "fbb218c0-5317-5bc6-957e-2ee96dd4b1f0"
version = "0.3.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BlockArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "54cd829dd26330c42e1cf9df68470dd4df602c61"
uuid = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
version = "0.16.38"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "cd67fc487743b2f0fd4380d4cbd3a24660d0eec8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.3"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "5225c965635d8c21168e32a12954675e7bea1151"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.10"

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

    [deps.Distances.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "b4fbdd20c889804969571cc589900803edda16b7"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.7.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "58d83dd5a78a36205bdfddb82b1bb67682e64487"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "0.4.9"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "35f0c0f345bff2c6d636f95fdb136323b5a796ef"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.7.0"
weakdeps = ["SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "c6e4a1fbe73b31a3dea94b1da449503b8830c306"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.21.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Gridap]]
deps = ["AbstractTrees", "BSON", "BlockArrays", "Combinatorics", "DataStructures", "DocStringExtensions", "FastGaussQuadrature", "FileIO", "FillArrays", "ForwardDiff", "JLD2", "JSON", "LineSearches", "LinearAlgebra", "NLsolve", "NearestNeighbors", "PolynomialBases", "QuadGK", "Random", "SparseArrays", "SparseMatricesCSR", "StaticArrays", "Test", "WriteVTK"]
git-tree-sha1 = "f504a1dc2a1474784e6d9854692d320be171f18f"
uuid = "56d4f2e9-7ea1-5844-9cf6-b9c51ca7ce8e"
version = "0.17.20"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ad37c091f7d7daf900963171600d7c1c5c3ede32"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "PrecompileTools", "Printf", "Reexport", "Requires", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "9bbb5130d3b4fa52846546bca4791ecbdfb52730"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.38"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "e129d9391168c677cd4800f5c0abb1ed8cb3794f"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL]]
deps = ["Artifacts", "Libdl", "LinearAlgebra", "MKL_jll"]
git-tree-sha1 = "100521a1d2181cb39036ee1a6955d6b9686bb363"
uuid = "33e6dc65-8f57-5167-99aa-e5a354878fb2"
version = "0.6.1"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "eb006abbd7041c28e0d16260e50a24f8f9104913"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2023.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "2c3726ceb3388917602169bed973dbc97f1b51a8"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.13"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PolynomialBases]]
deps = ["ArgCheck", "AutoHashEquals", "FFTW", "FastGaussQuadrature", "LinearAlgebra", "Requires", "SimpleUnPack", "SpecialFunctions"]
git-tree-sha1 = "aa1877430a7e8b0c7a35ea095c415d462af0870f"
uuid = "c74db56a-226d-5e98-8bb0-a6049094aeea"
version = "0.4.21"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9ebcd48c498668c7fa0e97a9cae873fbee7bfee1"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseMatricesCSR]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "38677ca58e80b5cad2382e5a1848f93b054ad28d"
uuid = "a0a7dd2c-ebf4-11e9-1f05-cf50bc540ca1"
version = "0.6.7"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "0adf069a2a490c47273727e029371b31d44b72b2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.5"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VTKBase]]
git-tree-sha1 = "c2d0db3ef09f1942d08ea455a9e252594be5f3b6"
uuid = "4004b06d-e244-455f-a6ce-a5f9919cc534"
version = "1.0.1"

[[deps.WriteVTK]]
deps = ["Base64", "CodecZlib", "FillArrays", "LightXML", "TranscodingStreams", "VTKBase"]
git-tree-sha1 = "41f0dc2a8f6fd860c266b91fd5cdf4fead65ae69"
uuid = "64499a7a-5c06-52f2-abe2-ccb03c286192"
version = "1.18.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "24b81b59bd35b3c42ab84fa589086e19be919916"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.11.5+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─32bfbf9c-fa5b-4f99-ad83-c612c9b0efe2
# ╟─81d7ad67-3362-4154-be4e-25f1a3faa35a
# ╟─cb8e8bc1-c2c5-4c62-8124-1d0fc520f181
# ╟─a177b3b1-e4a3-446b-94b5-aff5634e6571
# ╟─ec129a2a-0072-4a4a-b83c-21bda90176b3
# ╠═209f1308-7ee8-11ee-37a4-efac0b62e389
# ╟─1e6c7b24-524e-423c-9992-1d028abd0309
# ╠═c635e775-77d2-4f4a-8472-29aa7b4acb81
# ╠═034cdf88-ca69-4938-96b4-904487530cc2
# ╠═d8c681df-f4b8-482a-b1d5-49b0f9f02be8
# ╠═c56564a6-8d92-4b73-825f-a245d7950592
# ╠═f9c935e6-aae7-4eb3-90fd-2657b1a4f86e
# ╟─d9b279ea-1f36-45d4-9661-3cd87582426e
# ╠═64f868d6-e5c5-4b82-83a8-901c83ee1371
# ╠═34dbe9ce-a135-4880-8fdc-37e98c784f27
# ╠═ca063757-098a-4051-9f7b-75b4ca136fca
# ╟─a73ec73e-be63-450c-bb0f-ec43cca65621
# ╠═6dfd5b56-3162-4e9f-82fd-b53c88ec9a48
# ╠═2b4de08b-10a9-4d1b-8c92-b515a91f7702
# ╠═ded709b0-f743-4ba0-8eb7-9ecbbe28b2c6
# ╠═bafefb66-34a5-48e8-914c-9391953a310d
# ╠═b19e69cc-0fd2-4afb-8284-f9a18e2fc765
# ╟─fc85ccd9-0969-4688-80e1-3b9d072a318b
# ╟─6f78b279-879c-46ae-a7ca-ba06f3a847b0
# ╠═5029fc67-4660-4e38-b4ba-87cffe858e75
# ╠═ce06b4cd-28b8-4617-96b4-0d0eb0e4794b
# ╠═9a3844f0-f958-4732-93c6-be434c19149f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
