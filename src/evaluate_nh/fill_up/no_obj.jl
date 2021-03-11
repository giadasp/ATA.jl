#warmup feasibility
function fill_upᵥ(
    NH::Neighbourhood,
    v::Int64,
    iu::ItemUse,
    fs::FriendSets,
    constraints::Constraint,
    x_forced0ₜ::Vector{Bool},
    n_items::Int64,
    n_fs::Int64,
    ol_maxₜ::Vector{Float64},
)
    f₀, infeas₀ = Inf, Inf
    T = size(NH.x, 2)
    x₀ = copy(NH.x)
    xₜ = x₀[:, v]
    idxₜ₂ = findall(iszero.(xₜ)) #i₂ = 1, ..., I
    idxₜ₂ = Random.shuffle!(idxₜ₂) #removed
    i⁺ = copy(idxₜ₂[1])
    iu⁺ = 0
    iu₀ = sum(x₀, dims = 2)[:, 1] - iu.max
    ol₀ₜ = 0
    for i₂ in idxₜ₂
        if x_forced0ₜ[i₂]
            x₁, infeas₁, iu = copy(x₀), copy(infeas₀), copy(iu₀)
            x₁[i₂, v] = one(Float64)
            iu[i₂] = iu[i₂] + 1
            iu = sum_pos(iu)
            xₜ = copy(x₁[:, v])
            infeas₁, x_Iₜ = check_feas(fs, constraints, xₜ, n_fs, n_items)
            ol = eval_overlapᵥ(xₜ, x₁, fs.counts, ol_maxₜ, v)
            f₁ = infeas₁ + iu + ol
            if (f₁ < f₀)
                f₀, infeas₀ = copy(f₁), copy(infeas₁)
                iu⁺ = copy(iu)
                ol₀ₜ = copy(ol)
                i⁺ = copy(i₂)
            end
        end
    end
    NH.ol[v] = copy(ol₀ₜ)
    NH.infeas[v] = copy(infeas₀)
    NH.iu = copy(iu⁺)
    NH.x[i⁺, v] = one(Float64)
    return NH::Neighbourhood
end
