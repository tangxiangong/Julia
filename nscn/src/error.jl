using Gridap


function norm_error(u₁, u₂, dΩ, scarlar_or_vector="scalar")
    @assert scarlar_or_vector in ["scalar", "vector"]
    e = u₁ - u₂
    if scarlar_or_vector == "scalar"
        L²_error = ∫(e * e)dΩ |> sum |> sqrt
    else
        L²_error = ∫(e ⋅ e)dΩ |> sum |> sqrt
    end
    L²_error
end
