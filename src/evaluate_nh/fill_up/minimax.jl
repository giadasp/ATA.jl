#warmup MINIMAX
function find_best_itemᵥ(
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
    to_apply::Vector{Bool},
)
    f₀, TIF₀, infeas₀ = Inf, copy(NH.obj), Inf
    T = size(NH.x, 2)
    x₀ = copy(NH.x)
    idxₜ₂ = findall(iszero.(x₀[:, v])) #i₂ = 1, ..., I
    idxₜ₂ = Random.shuffle!(idxₜ₂) # !removed
    i⁺ = 1
    if to_apply[1]
        iu_max₀ = sum(NH.x, dims = 2)[:, 1] - iu.max
    else
        iu_max₀ = zeros(n_fs)
    end
    if to_apply[2] 
        iu_min₀ = - sum(NH.x, dims = 2)[:, 1] + iu.min
    else
        iu_min₀ = zeros(n_fs)
    end
    ol₀ₜ = 0
        iu_max⁺ = 0
    iu_min⁺ = 0
    for i₂ in idxₜ₂
        if x_forced0ₜ[i₂]
            x₁, TIF₁, infeas₁, iu_min, iu_max = copy(x₀), copy(TIF₀), copy(infeas₀), copy(iu_min₀), copy(iu_max₀)
            x₁[i₂, v] = one(Float64)
            xₜ = copy(x₁[:, v])
            if to_apply[1]
                iu_max[i₂] = iu_max[i₂] + 1    
                iu_max = sum_pos(iu_max)
            end        
            if to_apply[2]
                iu_min[i₂] = iu_min[i₂] - 1
                iu_min = sum_pos(iu_min)
            end 
            iu = sum_pos(iu)
            infeas₁, x_Iₜ = check_feas(fs, constraints, xₜ, n_fs, n_items)
            TIF₁[v] = eval_TIFₜ(x_Iₜ, coreₜ)
            if to_apply[3]
                ol = eval_overlapᵥ(xₜ, x₁, fs.counts, ol_maxₜ, v)
            else
                ol = 0
            end
            f₁ = (1 - opt_feas) * (infeas₁ + iu_min + iu_max + ol) - opt_feas * minimum(TIF₁)
            if (f₁ < f₀)
                i⁺ = copy(i₂)
                f₀, TIF₀, infeas₀ = copy(f₁), copy(TIF₁), copy(infeas₁)
                iu_max⁺ = copy(iu_max)
                iu_min⁺ = copy(iu_min)
                ol₀ₜ = copy(ol)
            end
        end
    end
    NH.ol[v] = copy(ol₀ₜ)#0
    NH.infeas[v] = copy(infeas₀)
    NH.iu = iu_min⁺ + iu_max⁺
    NH.obj = copy(TIF₀)
    NH.x[i⁺, v] = one(Float64)
    return NH::Neighbourhood
end
