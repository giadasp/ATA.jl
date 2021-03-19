#warmup custom
function find_best_itemᵥ(
    NH::Neighbourhood,
    obj::CustomObjective,
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
    f₀, obj_fun₀, infeas₀ = Inf, copy(NH.obj), Inf
    T = size(NH.x, 2)
    x₀ = copy(NH.x)
    idxₜ₂ = findall(iszero.(x₀[:, v])) #i₂ = 1, ..., I
    idxₜ₂ = Random.shuffle!(idxₜ₂) # !removed
    i⁺ = 1 
    if to_apply[1]
        iu_max₀ = sum(NH.x, dims = 2)[:, 1] - iu.max
    else
        iu_max₀ = 0
    end
    if to_apply[2]
        iu_min₀ = - sum(NH.x, dims = 2)[:, 1] + iu.min
    else
        iu_min₀ = 0
    end
    iu_max⁺ = 0
    iu_min⁺ = 0
    ol₀ₜ = 0
    for i₂ in idxₜ₂
        if x_forced0ₜ[i₂]
            x₁, obj_fun₁, infeas₁, iu_min, iu_max = copy(x₀), copy(obj_fun₀), copy(infeas₀), copy(iu_min₀), copy(iu_max₀)
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
            infeas₁, x_Iₜ = check_feas(fs, constraints, xₜ, n_fs, n_items)
            obj_fun₁ = obj.fun(x_Iₜ, obj.args)
            if to_apply[3]
                ol = eval_overlapᵥ(xₜ, x₁, fs.counts, ol_maxₜ, v)
            else
                ol = 0
            end
            f₁ = (1 - opt_feas) * (infeas₁ + iu_min + iu_max + ol) - opt_feas * minimum(obj_fun₁)
            if (f₁ < f₀)
                i⁺ = copy(i₂)
                f₀, obj_fun₀, infeas₀ = copy(f₁), copy(obj_fun₁), copy(infeas₁)
                iu_max⁺ = copy(iu_max)
                iu_min⁺ = copy(iu_min)
                ol₀ₜ = copy(ol)
            end
        end
    end
    NH.ol[v] = copy(ol₀ₜ)#0
    NH.infeas[v] = copy(infeas₀)
    NH.iu = iu_min⁺ + iu_max⁺
    NH.obj = copy(obj_fun₀)
    NH.x[i⁺, v] = one(Float64)
    return NH::Neighbourhood
end
