function eval_overlapᵥ(
    xₜ::Vector{Float64},
    x::Matrix{Float64},
    fs_counts::Vector{Int64},
    ol_maxₜ::Vector{Float64},
    v::Int64,
)
    xₜ = xₜ .* fs_counts
    cons = LinearAlgebra.BLAS.gemv('T', x, xₜ) - ol_maxₜ
    cons[v] = 0
    return sum_pos(cons)::Float64
end

function eval_overlap(
    x::Matrix{Float64},
    fs_counts::Matrix{Float64},
    ol_max::Matrix{Float64},
    T::Int64,
    ol2::Vector{Float64},
)
    fs_counts = fs_counts .* x
    ol = _gemmblas(x, fs_counts, ol_max)
    for v2 = 1:T
        olₜ = ol[:, v2]
        olₜ[v2] = zero(Float64)
        ol2[v2] = sum_pos(olₜ)
    end
    return ol2::Vector{Float64}
end

function eval_exp_scoreₜ(x_Iₜ::Vector{Float64}, expected_scoreₜ::ExpectedScore) #x = I
    es = (expected_scoreₜ.val * x_Iₜ) ./ sum(x_Iₜ)
    # _gemvblasT(
    #     expected_scoreₜ.val,
    #     x_Iₜ,
    #     size(expected_scoreₜ.val, 2)
    # ) ./ sum(x_Iₜ)
    es = vcat((es .- expected_scoreₜ.max)..., (-es .+ expected_scoreₜ.min)...)
    return es
end

function eval_cons(xₜ::Vector{Float64}, Aₜ::Matrix{Float64}, bₜ::Vector{Float64})
    cons = copy(bₜ)
    cons = _gemvblas(Aₜ, xₜ, cons, size(xₜ, 1))
    return cons::Vector{Float64}
end

function check_feas(
    fs::FriendSets,
    consₜ::Constraint,
    xₜ::Vector{Float64},
    n_fs::Int64,
    n_items::Int64,
)
    es = 0
    cons = 0
    if size(consₜ.constr_A, 1) > 0
        cons = eval_cons(xₜ, consₜ.constr_A, consₜ.constr_b)
    end
    if n_items > n_fs
        x_Iₜ = fs_to_items(xₜ, n_items, fs.items)
    else
        x_Iₜ = copy(xₜ)
    end
    if size(consₜ.expected_score.val, 1) > 0
        cons = vcat(cons, eval_exp_scoreₜ(x_Iₜ, consₜ.expected_score))
    end
    return sum_pos(cons)::Float64, x_Iₜ::Vector{Float64}
end
