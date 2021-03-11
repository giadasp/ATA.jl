#warmup MINIMAX
function fill_upᵥ(
    NH::Neighbourhood,
    coreₜ::MinimaxObjectiveCore,
    opt_feas::Float64,
    v::Int64,
    iu::ItemUse,
    fs::FriendSets,
    constraints::Constraint,
    x_forced0ₜ::Vector{Bool},
    n_items::Int64,
    n_fs::Int64,
    ol_maxₜ::Vector{Float64},
)
    f₀, TIF₀, infeas₀ = Inf, copy(NH.obj), Inf
    T = size(NH.x, 2)
    x₀ = copy(NH.x)
    idxₜ₂ = findall(iszero.(x₀[:, v])) #i₂ = 1, ..., I
    idxₜ₂ = Random.shuffle!(idxₜ₂) # !removed
    i⁺ = 1
    iu₀ = sum(NH.x, dims = 2)[:, 1] - iu.max
    ol₀ₜ = 0
    iu⁺ = 0
    for i₂ in idxₜ₂
        if x_forced0ₜ[i₂]
            x₁, TIF₁, infeas₁, iu = copy(x₀), copy(TIF₀), copy(infeas₀), copy(iu₀)
            x₁[i₂, v] = one(Float64)
            xₜ = copy(x₁[:, v])
            iu[i₂] = iu[i₂] + 1
            iu = sum_pos(iu)
            infeas₁, x_Iₜ = check_feas(fs, constraints, xₜ, n_fs, n_items)
            TIF₁[v] = eval_TIFₜ(x_Iₜ, coreₜ)
            ol = eval_overlapᵥ(xₜ, x₁, fs.counts, ol_maxₜ, v)
            f₁ = (1 - opt_feas) * (infeas₁ + iu + ol) - opt_feas * minimum(TIF₁)
            if (f₁ < f₀)
                i⁺ = copy(i₂)
                f₀, TIF₀, infeas₀ = copy(f₁), copy(TIF₁), copy(infeas₁)
                iu⁺ = copy(iu)
                ol₀ₜ = copy(ol)
            end
        end
    end
    NH.ol[v] = copy(ol₀ₜ)#0
    NH.infeas[v] = copy(infeas₀)
    NH.iu = copy(iu⁺)
    NH.obj = copy(TIF₀)
    NH.x[i⁺, v] = one(Float64)
    return NH::Neighbourhood
end
