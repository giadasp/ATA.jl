function eval_overlap_v(
    xₜ::Vector{Float64},
    x::Matrix{Float64},
    FS_counts::Vector{Int64},
    ol_maxₜ::Vector{Float64},
    v::Int64,
)
    xₜ = xₜ .* FS_counts
    cons = LinearAlgebra.BLAS.gemv('T', x, xₜ) - ol_maxₜ
    cons[v] = 0
    return sum(cons[cons.>0])
end

function eval_overlap(
    x::Matrix{Float64},
    FS_counts::Matrix{Float64},
    ol_max::Matrix{Float64},
    T::Int64,
    ol2::Vector{Float64},
)
    FS_counts = FS_counts .* x
    ol = _gemmblas(x, FS_counts, ol_max)
    for v2 = 1:T
        olₜ = view(ol, :, v2)
        olₜ[v2] = zero(Float64)
        olₜ = olₜ[olₜ.>0]
        if size(olₜ, 1) == 0
            ol2[v2] = 0
        else
            ol2[v2] = sum(olₜ)
        end
    end
    return ol2::Vector{Float64}
end

function eval_TIF_CCₜ(xₜ::Vector{Float64}, IIF::Array{Float64,3}; α = 0.1) # x = I
    K, I, R = size(IIF)
    alphaR = Int(ceil(R * (α)))
    αQle = Inf
    if K > 1
        for k = 1:K
            αQle = min(αQle, sort(LinearAlgebra.BLAS.gemv('T', IIF[k, :, :], xₜ))[alphaR])
        end
    else
        αQle = sort(_gemvblasT(IIF[1, :, :], xₜ, R))[alphaR]
    end
    return αQle
end

function eval_Exp_Scoreₜ(x_Iₜ::Vector{Float64}, expected_scoreₜ::ExpectedScore) #x = I
    # es =
    #     _gemvblas(
    #         expected_scoreₜ.val,
    #         x_Iₜ,
    #         zeros(size(expected_scoreₜ.val, 1)),
    #         size(x_Iₜ, 1),
    #     ) / sum(x_Iₜ)
    es =
        _gemvblasT(
            expected_scoreₜ.val,
            x_Iₜ,
            size(expected_scoreₜ.val, 2)
        ) / sum(x_Iₜ)
    es = vcat((es .- expected_scoreₜ.max)..., (-es .+ expected_scoreₜ.min)...)
    return es
end

function eval_TIF_MMₜ(xₜ::Vector{Float64}, IIF::Matrix{Float64}) # x = I
    return minimum(LinearAlgebra.BLAS.gemv('N', IIF, xₜ))::Float64
end

function eval_TIF_mmₜ(xₜ::Vector{Float64}, IIF::Matrix{Float64}, targetsₜ::Vector{Float64}) # x = I
    return -maximum(abs.(LinearAlgebra.BLAS.gemv('N', IIF, xₜ) - targetsₜ))::Float64
end
# ! not used anymore, item use is evaluated item by item
# function eval_IU(x::Matrix{Float64}, T::Int64) #x = I*T
# 	ItemUse = x * ones(T)
# 	return ItemUse
# end

function eval_cons(xₜ::Vector{Float64}, Aₜ::Matrix{Float64}, bₜ::Vector{Float64})
    cons = copy(bₜ)
    cons = _gemvblas(Aₜ, xₜ, cons, size(xₜ, 1))
    return cons::Vector{Float64}
end

function check_feas(
    FS::FS,
    Consₜ::Constraint,
    xₜ::Vector{Float64},
    n_FS::Int64,
    n_items::Int64,
)
    es = 0
    cons = 0
    if size(Consₜ.constr_A, 1) > 0
        cons = eval_cons(xₜ, Consₜ.constr_A, Consₜ.constr_b)
    end
    if n_items > n_FS
        x_Iₜ = FS_to_items(xₜ, n_items, FS.items)
    else
        x_Iₜ = copy(xₜ)
    end
    if size(Consₜ.expected_score.val, 1) > 0
        cons = vcat(cons, eval_Exp_Scoreₜ(x_Iₜ, Consₜ.expected_score))
    end
    cons = cons[cons.>0]
    if size(cons, 1) > 0
        cons = sum(cons)
    else
        cons = zero(Float64)
    end
    return cons::Float64, x_Iₜ::Vector{Float64}
end


#warmup feasibility
function fill_up_feas(
    NH::Neighbourhood,
    v::Int64,
    IU::IU,
    FS::FS,
    constraints::Constraint,
    x_forced0ₜ::Vector{Bool},
    n_items::Int64,
    n_FS::Int64,
    ol_maxₜ::Vector{Float64},
)
    f₀, infeas₀ = Inf, Inf
    T = size(NH.x, 2)
    x₀ = copy(NH.x)
    xₜ = x₀[:, v]
    idxₜ₂ = setdiff(collect(1:n_FS), findall(xₜ .== one(Float64))) #i₂ = 1, ..., I
    #idxₜ₂ = Random.shuffle!(idxₜ₂) #removed
    i⁺ = copy(idxₜ₂[1])
    iu⁺ = 0
    iu₀ = sum(x₀, dims = 2) - IU.max
    ol₀ₜ = 0
    for i₂ in idxₜ₂
        if x_forced0ₜ[i₂]
            x₁, infeas₁, iu = copy(x₀), copy(infeas₀), copy(iu₀)
            x₁[i₂, v] = one(Float64)
            iu[i₂] = iu[i₂] + 1
            iu = iu[iu.>0]
            if size(iu, 1) == 0
                iu = 0
            else
                iu = sum(iu)
            end
            xₜ = copy(x₁[:, v])
            infeas₁, x_Iₜ = check_feas(FS, constraints, xₜ, n_FS, n_items)
            ol = eval_overlap_v(xₜ, x₁, FS.counts, ol_maxₜ, v)
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

#warmup MAXIMIN
function fill_up_MAXIMIN(
    NH::Neighbourhood,
    IIFₜ::Matrix{Float64},
    opt_feas::Float64,
    v::Int64,
    IU::IU,
    FS::FS,
    constraints::Constraint,
    x_forced0ₜ::Vector{Bool},
    n_items::Int64,
    n_FS::Int64,
    ol_maxₜ::Vector{Float64},
)
    f₀, TIF₀, infeas₀ = Inf, copy(NH.obj), Inf
    T = size(NH.x, 2)
    x₀ = copy(NH.x)
    idxₜ₂ = setdiff(collect(1:n_FS), findall(x₀[:, v] .== one(Float64))) #i₂ = 1, ..., I
    #idxₜ₂ = Random.shuffle!(idxₜ₂) # !removed
    i⁺ = 1
    iu₀ = sum(NH.x, dims = 2) - IU.max
    ol₀ₜ = 0
    iu⁺ = 0
    for i₂ in idxₜ₂
        if x_forced0ₜ[i₂]
            x₁, TIF₁, infeas₁, iu = copy(x₀), copy(TIF₀), copy(infeas₀), copy(iu₀)
            x₁[i₂, v] = one(Float64)
            xₜ = copy(x₁[:, v])
            iu[i₂] = iu[i₂] + 1
            iu = iu[iu.>0]
            if size(iu, 1) == 0
                iu = 0
            else
                iu = sum(iu)
            end
            infeas₁, x_Iₜ = check_feas(FS, constraints, xₜ, n_FS, n_items)
            TIF₁[v] = eval_TIF_MMₜ(x_Iₜ, IIFₜ)
            ol = eval_overlap_v(xₜ, x₁, FS.counts, ol_maxₜ, v)
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

#warmup MAXIMIN
function fill_up_MINIMAX(
    NH::Neighbourhood,
    IIFₜ::Matrix{Float64},
    opt_feas::Float64,
    v::Int64,
    IU::IU,
    FS::FS,
    constraints::Constraint,
    x_forced0ₜ::Vector{Bool},
    n_items::Int64,
    n_FS::Int64,
    ol_maxₜ::Vector{Float64},
)
    f₀, TIF₀, infeas₀ = Inf, copy(NH.obj), Inf
    T = size(NH.x, 2)
    x₀ = copy(NH.x)
    idxₜ₂ = setdiff(collect(1:n_FS), findall(x₀[:, v] .== one(Float64))) #i₂ = 1, ..., I
    #idxₜ₂ = Random.shuffle!(idxₜ₂) # !removed
    i⁺ = 1
    iu₀ = sum(NH.x, dims = 2) - IU.max
    ol₀ₜ = 0
    iu⁺ = 0
    for i₂ in idxₜ₂
        if x_forced0ₜ[i₂]
            x₁, TIF₁, infeas₁, iu = copy(x₀), copy(TIF₀), copy(infeas₀), copy(iu₀)
            x₁[i₂, v] = one(Float64)
            xₜ = copy(x₁[:, v])
            iu[i₂] = iu[i₂] + 1
            iu = iu[iu.>0]
            if size(iu, 1) == 0
                iu = 0
            else
                iu = sum(iu)
            end
            infeas₁, x_Iₜ = check_feas(FS, constraints, xₜ, n_FS, n_items)
            TIF₁[v] = eval_TIF_mmₜ(x_Iₜ, IIFₜ)
            ol = eval_overlap_v(xₜ, x₁, FS.counts, ol_maxₜ, v)
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

#warmup MAXIMIN CC
function fill_up_CC(
    NH::Neighbourhood,
    α::Float64,
    IIFₜ::Array{Float64,3},
    opt_feas::Float64,
    v::Int64,
    IU::IU,
    FS::FS,
    constraints::Constraint,
    x_forced0ₜ::Vector{Bool},
    n_items::Int64,
    n_FS::Int64,
    ol_maxₜ::Vector{Float64},
)
    f₀, TIF₀, infeas₀ = Inf, copy(NH.obj), Inf
    T = size(NH.x, 2)
    x₀ = copy(NH.x)
    idxₜ₂ = setdiff(collect(1:n_FS), findall(x₀[:, v] .== one(Float64))) #i₂ = 1, ..., I
    #idxₜ₂ = Random.shuffle!(idxₜ₂) #! removed
    i⁺ = 1
    ol₀ₜ = 0
    iu₀ = sum(NH.x, dims = 2) - IU.max
    iu⁺ = 0
    for i₂ in idxₜ₂
        if x_forced0ₜ[i₂]
            x₁, TIF₁, infeas₁, iu = copy(x₀), copy(TIF₀), copy(infeas₀), copy(iu₀)
            x₁[i₂, v] = one(Float64)
            xₜ = copy(x₁[:, v])
            iu[i₂] = iu[i₂] + 1
            iu = iu[iu.>0]
            if size(iu, 1) == 0
                iu = 0
            else
                iu = sum(iu)
            end
            infeas₁, x_Iₜ = check_feas(FS, constraints, xₜ, n_FS, n_items)
            TIF₁[v] = eval_TIF_CCₜ(x_Iₜ, IIFₜ; α = α)
            ol = eval_overlap_v(xₜ, x₁, FS.counts, ol_maxₜ, v)
            f₁ = (1 - opt_feas) * (infeas₁ + iu + ol) - opt_feas * minimum(TIF₁)
            if (f₁ < f₀)
                i⁺ = copy(i₂)
                f₀, TIF₀, infeas₀ = copy(f₁), copy(TIF₁), copy(infeas₁)
                iu⁺ = copy(iu)
                ol₀ₜ = copy(ol)
            end
        end
    end
    NH.ol[v] = copy(ol₀ₜ)
    NH.infeas[v] = copy(infeas₀)
    NH.iu = copy(iu⁺)
    NH.obj = copy(TIF₀)
    NH.x[i⁺, v] = one(Float64)
    return NH::Neighbourhood
end

#warmup custom
function fill_up_custom(
    NH::Neighbourhood,
    obj::Obj,
    opt_feas::Float64,
    v::Int64,
    IU::IU,
    FS::FS,
    constraints::Constraint,
    x_forced0ₜ::Vector{Bool},
    n_items::Int64,
    n_FS::Int64,
    ol_maxₜ::Vector{Float64},
)
    f₀, obj_fun₀, infeas₀ = Inf, copy(NH.obj), Inf
    T = size(NH.x, 2)
    x₀ = copy(NH.x)
    idxₜ₂ = setdiff(collect(1:n_FS), findall(x₀[:, v] .== one(Float64))) #i₂ = 1, ..., I
    #idxₜ₂ = Random.shuffle!(idxₜ₂) # !removed
    i⁺ = 1
    iu₀ = sum(NH.x, dims = 2) - IU.max
    ol₀ₜ = 0
    iu⁺ = 0
    for i₂ in idxₜ₂
        if x_forced0ₜ[i₂]
            x₁, obj_fun₁, infeas₁, iu = copy(x₀), copy(obj_fun₀), copy(infeas₀), copy(iu₀)
            x₁[i₂, v] = one(Float64)
            xₜ = copy(x₁[:, v])
            iu[i₂] = iu[i₂] + 1
            iu = iu[iu.>0]
            if size(iu, 1) == 0
                iu = 0
            else
                iu = sum(iu)
            end
            infeas₁, x_Iₜ = check_feas(FS, constraints, xₜ, n_FS, n_items)
            obj_fun₁ = obj.fun(x_Iₜ, obj.args)
            ol = eval_overlap_v(xₜ, x₁, FS.counts, ol_maxₜ, v)
            f₁ = (1 - opt_feas) * (infeas₁ + iu + ol) - opt_feas * minimum(obj_fun₁)
            if (f₁ < f₀)
                i⁺ = copy(i₂)
                f₀, obj_fun₀, infeas₀ = copy(f₁), copy(obj_fun₁), copy(infeas₁)
                iu⁺ = copy(iu)
                ol₀ₜ = copy(ol)
            end
        end
    end
    NH.ol[v] = copy(ol₀ₜ)#0
    NH.infeas[v] = copy(infeas₀)
    NH.iu = copy(iu⁺)
    NH.obj = copy(obj_fun₀)
    NH.x[i⁺, v] = one(Float64)
    return NH::Neighbourhood
end
#MAXIMIN CC neighbourhood
function analyse_NH(
    NH_start::Neighbourhood,
    ATAmodel::Model,
    IIF::Vector{Array{Float64,3}};
    fF = true,
    n_fill = 1,
    opt_feas = 0.9,
    max_conv = 1,
    start_temp = 1000.0,
    geom_temp = 0.1,
    n_item_sample = 1,
    n_test_sample = 1,
    verbosity = 1,
    start_time = 0,
    max_time = 1000,
)
    if fF == true
        NH_start.obj = zeros(Float64, ATAmodel.settings.T)
    end
    NH₁ = Neighbourhood()
    NH₁ = _mycopy(NH_start, NH₁)
    NH₀ = Neighbourhood()
    NH⁺ = Neighbourhood()
    f_star = ones(2) .* Inf
    f_evals = 0
    t = copy(start_temp)
    T = ATAmodel.settings.T
    n_items = ATAmodel.settings.n_items
    n_FS = ATAmodel.settings.n_FS
    FS_counts = ATAmodel.settings.FS.counts * ones(Float64, T)'
    switches = 0
    removes = 0
    adds = 0
    #warm up
    println("Fill-up starting...")
    round = 1

    if n_fill > 0
        for round = 1:n_fill
            NH₁ = _mycopy(NH_start, NH₁)
            warmup = fill(true, T)
            while any(warmup) #filling forms
                v = findfirst(warmup .== true)
                constraints = ATAmodel.constraints[v]
                #println("ol = ", NH₁.ol)
                fₜ = (1 - opt_feas) * (NH₁.infeas + NH₁.ol) - (opt_feas * NH₁.obj)
                mm = fₜ[v]
                for v2 in findall(warmup .== true)
                    if fₜ[v2] > mm
                        mm = copy(fₜ[v2])
                        v = copy(v2)
                    end
                end
                #filling
                n_t = LinearAlgebra.dot(NH₁.x[:, v], ATAmodel.settings.FS.counts)
                #try to add other items, the first time it goes over n_max it stops
                if n_t < constraints.length_max
                    if opt_feas == 0 || fF == true
                        NH_add = fill_up_feas(
                            NH₁,
                            v,
                            ATAmodel.settings.IU,
                            ATAmodel.settings.FS,
                            constraints,
                            ATAmodel.settings.forced0[v],
                            n_items,
                            n_FS,
                            ATAmodel.settings.ol_max[:, v],
                        )
                        NH_add.f = opt_feas * NH_add.f
                    else
                        NH_add = fill_up_CC(
                            NH₁,
                            ATAmodel.obj.aux_float,
                            IIF[v],
                            opt_feas,
                            v,
                            ATAmodel.settings.IU,
                            ATAmodel.settings.FS,
                            constraints,
                            ATAmodel.settings.forced0[v],
                            n_items,
                            n_FS,
                            ATAmodel.settings.ol_max[:, v],
                        )
                    end
                    NH_add.ol = eval_overlap(
                        NH_add.x,
                        FS_counts,
                        ATAmodel.settings.ol_max,
                        T,
                        NH_add.ol,
                    )
                    Printf.@printf "."
                    nₜ = LinearAlgebra.dot(NH_add.x[:, v], ATAmodel.settings.FS.counts)
                    #println("length for test ", v, ": ", nₜ)
                    if n_t <= constraints.length_max
                        NH₁ = _mycopy(NH_add, NH₁)
                    else
                        warmup[v] = false
                        println("f₁:", NH₁.f)
                        println("-Test ", v, " filled up with ", n_t, " items, ")
                    end
                else
                    warmup[v] = false
                    println("-Test ", v, " filled up with ", n_t, " items, ")
                end
            end
        end #end of round, round = nRound
        NH₀ = _mycopy(NH₁, NH₀)
        NH₀.f = _comp_f(NH₀, opt_feas)
    end
    NH⁺ = _mycopy(NH₀, NH⁺)
    println("End of fill-up")
    if sum(NH₀.infeas + NH₀.ol) + NH₀.iu == 0
        println("Feasible solution found in fill-up")
        fF = false
    else
        println("Feasible solution not found in fill-up")
    end
    println("f = ", NH⁺.f)
    println("infeas = ", NH⁺.infeas)
    println("ol = ", NH⁺.ol)
    println("iu = ", NH⁺.iu)
    # statistics to repogeom_temp at each temp change, set back to zero
    coverage_ok = 0
    if n_test_sample > T
        n_test_sample = T
    end
    convergence = 0

    while coverage_ok == 0
        NH₁ = _mycopy(NH₀, NH₁)
        weights = (1 - opt_feas) .* (NH₀.infeas + NH₀.ol) - opt_feas .* (NH₀.obj)
        weights = (weights .- minimum(weights)) .+ 1
        weights = weights ./ sum(weights)
        weights = StatsBase.ProbabilityWeights(weights)
        #determine test order
        #testOrder = StatsBase.sample(collect(1:T), weights, n_test_sample, replace = false)
        #testOrder = Random.shuffle!(collect(1:ATAmodel.settings.T))
        testOrder = sortperm(weights, rev = true)
        #println("test order = ", Int.(testOrder))
        v₂ = 0
        exit = 0
        #iteratorTestItem = vec(collect(Iterators.product(collect(1:n_item_sample), testOrder[1:n_test_sample])))
        #println(iteratorTestItem[1:nItemoStatsBase.sample])
        #it = 0
        #xnew = copy(NH₀.x)
        while exit == 0 && v₂ < n_test_sample #it<size(iteratorTestItem, 1) #
            #it+= 1
            v₂ += 1
            exit = 0
            v = testOrder[v₂]
            x_forced0ₜ = ATAmodel.settings.forced0[v]
            Constraintsₜ = ATAmodel.constraints[v]
            IIFₜ = IIF[v]
            ol_maxₜ = ATAmodel.settings.ol_max[:, v]
            #it<size(iteratorTestItem, 1) #
            #v = iteratorTestItem[it][2]
            #NH₀.x = copy(xnew)
            taken_items = findall(NH₀.x[:, v] .== 1)
            if n_item_sample > size(taken_items, 1)
                nI = Int(size(taken_items, 1))
            else
                nI = n_item_sample
            end
            #taken_items = Random.shuffle!(taken_items) #reset #removed
            # exit2 = 0
            # if iteratorTestItem[it][1]>size(taken_items, 1)
            # 	exit2 = 1
            # end
            #if exit2 == 0
            # if v!= v2
            # 	v = copy(v2)
            #
            # end
            #h = taken_items[iteratorTestItem[it][1]]
            #println("test ", v₂, " of ", size(testOrder, 1))
            #fix test features
            h₂ = 0
            while exit == 0 && h₂ < nI
                add_remove = 0
                h₂ += 1
                #println("item ", h₂, " of ", size(taken_items, 1))
                h = taken_items[h₂]
                while exit == 0 && add_remove < 2 
                    add_remove += 1
                    NH₁ = _mycopy(NH₀, NH₁)
                    #try to remove h
                    #iu = copy(iu₀)
                    #if rand()>0.0
                    if add_remove == 1 #remove
                        NH₁.x[h, v] = zero(Float64)
                    end#else add (not remove)
                    #iu[h]-= 1
                    taken = 0
                    #else
                    #	taken = 1
                    #end
                    if (
                        add_remove == 1 &&
                        sum(NH₁.x[:, v] .* ATAmodel.settings.FS.counts) >=
                        Constraintsₜ.length_min
                    ) || (
                        add_remove == 2 &&
                        sum(NH₁.x[:, v] .* ATAmodel.settings.FS.counts) <
                        Constraintsₜ.length_max
                    )
                        NH₁.infeas[v], x_Iₜ = check_feas(
                            ATAmodel.settings.FS,
                            Constraintsₜ,
                            NH₁.x[:, v],
                            n_FS,
                            n_items,
                        )
                        iu = sum(NH₁.x, dims = 2) - ATAmodel.settings.IU.max
                        iu = iu[iu.>0]
                        if size(iu, 1) == 0
                            NH₁.iu = 0
                        else
                            NH₁.iu = sum(iu)
                        end
                        NH₁.ol = eval_overlap(
                            NH₁.x,
                            FS_counts,
                            ATAmodel.settings.ol_max,
                            T,
                            NH₁.ol,
                        )
                        #NH₁.ol[v] = eval_overlap_v(NH₁.x[:, v], NH₁.x, ATAmodel.settings.FS.counts, ol_maxₜ, v)
                        if fF == false
                            NH₁.obj[v] =
                                eval_TIF_CCₜ(x_Iₜ, IIFₜ; α = ATAmodel.obj.aux_float)
                        end
                        NH₁.f = _comp_f(NH₁, opt_feas)
                        f_evals += 1
                        if (NH₁.f <= NH₀.f)
                            #switch item
                            #NH₀.ol = eval_overlap(NH₀.x, FS_counts, ATAmodel.settings.ol_max, T, NH₀.ol)
                            #NH₀.f = _comp_f(NH₀, opt_feas)
                            NH₀ = _mycopy(NH₁, NH₀)
                            # println("better to remove, new f₀ = ", NH₀.f)
                            if (NH₁.f < NH⁺.f)
                                exit = 1
                                convergence = 0
                                NH⁺ = _mycopy(NH₁, NH⁺)
                                if verbosity == 2
                                    println("f⁺ = ", NH⁺.f, ", t = ", v)
                                else
                                    Printf.@printf "+"
                                end
                            end
                        else
                            p = _sig_c(-(NH₁.f - NH₀.f) / t) # -(NH₁.f - NH₀.f) / t always negative
                            if p < 1e-10
                                p = zeros(Float64)[1]
                            end

                            if (rand() < p)
                                #remove
                                #exit = 1
                                NH₀ = _mycopy(NH₁, NH₀)
                                #NH₀.ol = eval_overlap(NH₀.x, FS_counts, ATAmodel.settings.ol_max, T, NH₀.ol)
                                #NH₀.f = _comp_f(NH₀, opt_feas)
                                if verbosity == 2
                                    println("SA: f₀ = ", NH₀.f, ", t = ", v)
                                else
                                    Printf.@printf "_"
                                end
                            end
                        end
                    end
                    idxₜ₂ = findall(NH₁.x[:, v] .== zero(Float64))#Random.shuffle!(findall(NH₁.x[:, v] .== zero(Float64)))#i₂ = 1, ..., I
                    betterFound = 0
                    i₃ = 0
                    x₋₁ = copy(NH₁.x)
                    while betterFound == 0 && i₃ < size(idxₜ₂, 1)
                        i₃ += 1
                        i₂ = idxₜ₂[i₃]
                        if x_forced0ₜ[i₂] && i₂ != h
                            NH₁ = _mycopy(NH₀, NH₁)
                            NH₁.x = copy(x₋₁)
                            #NH₁.x[h, v] = zero(Float64)
                            NH₁.x[i₂, v] = one(Float64)
                            # iu = copy(iu₀)
                            # iu[i₂]+= 1
                            # iu[h]-= 1
                            iu = sum(NH₁.x, dims = 2) - ATAmodel.settings.IU.max
                            iu = iu[iu.>0]
                            if size(iu, 1) == 0
                                NH₁.iu = 0
                            else
                                NH₁.iu = sum(iu)
                            end
                            NH₁.infeas[v], x_Iₜ = check_feas(
                                ATAmodel.settings.FS,
                                Constraintsₜ,
                                NH₁.x[:, v],
                                n_FS,
                                n_items,
                            )
                            if fF == false
                                NH₁.obj[v] =
                                    eval_TIF_CCₜ(x_Iₜ, IIFₜ; α = ATAmodel.obj.aux_float)
                            end
                            NH₁.ol = eval_overlap(
                                NH₁.x,
                                FS_counts,
                                ATAmodel.settings.ol_max,
                                T,
                                NH₁.ol,
                            )
                            #NH₁.ol[v] = eval_overlap_v(NH₁.x[:, v], NH₁.x, ATAmodel.settings.FS.counts, ol_maxₜ, v)
                            NH₁.f = _comp_f(NH₁, opt_feas)
                            if (NH₁.f <= NH₀.f)
                                #switch item
                                #NH₀ = _mycopy(NH₁, NH₀)
                                # NH₀.ol = eval_overlap(NH₀.x, FS_counts, ATAmodel.settings.ol_max, T, NH₀.ol)
                                # NH₀.f = _comp_f(NH₀, opt_feas)
                                NH₀ = _mycopy(NH₁, NH₀)
                                # println("better to switch, new f₀ = ", NH₀.f)
                                if (NH₁.f < NH⁺.f)
                                    exit = 1
                                    betterFound = 1
                                    convergence = 0
                                    NH⁺ = _mycopy(NH₁, NH⁺)
                                    if verbosity == 2
                                        println("f⁺ = ", NH⁺.f, ", t = ", v)
                                    else
                                        Printf.@printf "+"
                                    end
                                end
                            else
                                p = _sig_c(-(NH₁.f - NH₀.f) / t)
                                if p < 1e-10
                                    p = zeros(Float64)[1]
                                end
                                if (rand() < p)
                                    #remove
                                    NH₀ = _mycopy(NH₁, NH₀)
                                    # NH₀.ol = eval_overlap(NH₀.x, FS_counts, ATAmodel.settings.ol_max, T, NH₀.ol)
                                    # NH₀.f = _comp_f(NH₀, opt_feas)
                                    if verbosity == 2
                                        println("SA: f₀ = ", NH₀.f, ", t = ", v)
                                    else
                                        Printf.@printf "_"
                                    end
                                end
                            end
                        end
                    end #end of betterFound (betterFound = 1)
                end #end of add_remove
            end #end of itemorder h₂ (exit = 1)

        end #end of testorder v₂ (exit = 1)
        if sum(NH₀.infeas + NH₀.ol) + NH₀.iu <= 0 && fF == true
            fF = false
            for v = 1:T
                x_Iₜ = FS_to_items(
                    NH₀.x[:, v],
                    ATAmodel.settings.n_items,
                    ATAmodel.settings.FS.items,
                )
                NH₀.obj[v] = eval_TIF_CCₜ(x_Iₜ, IIF[v]; α = ATAmodel.obj.aux_float)
            end
            NH₀.f = _comp_f(NH₀, opt_feas)
            NH⁺ = _mycopy(NH₀, NH⁺)
        end
        f_star[1] = copy(NH₀.f)
        #if exit == 0
        if f_star[2] == f_star[1]
            convergence += 1
            Printf.@printf("	%16.10f", convergence)
        end
        #println("convergence is ", convergence)
        # ? how many are equal f₀ in the last iterations?
        if (f_star[2] == f_star[1] && convergence == max_conv) ||
           time() - start_time >= max_time
            coverage_ok = 1
            #t += 1
        else
            t *= geom_temp
            coverage_ok = 0
            if f_star[1] <= f_star[2]
                f_star[2] = copy(f_star[1])
            end
        end
        if verbosity == 2
            println("\n")
            Printf.@printf("\n  f⁺:	%16.10f", NH⁺.f)
            Printf.@printf("\n  f₀:	%16.10f", NH₀.f)
            Printf.@printf("\n  Local convergence:	%16.10f", convergence)
            println("\n")
        end
    end #end of NH coverage (coverage_ok = 1)
    return NH⁺::Neighbourhood
end

#MAXIMIN or MINIMAX neighbourhood
function analyse_NH(
    NH_start::Neighbourhood,
    ATAmodel::Model,
    IIF::Vector{Array{Float64,2}};
    fF = true,
    n_fill = 1,
    opt_feas = 0.9,
    max_conv = 1,
    start_temp = 1000.0,
    geom_temp = 0.1,
    n_item_sample = 1,
    n_test_sample = 1,
    verbosity = 1,
    start_time = 0,
    max_time = 1000,
)
    if fF == true
        NH_start.obj = zeros(Float64, ATAmodel.settings.T)
    end
    NH₁ = Neighbourhood()
    NH₁ = _mycopy(NH_start, NH₁)
    NH₀ = Neighbourhood()
    NH⁺ = Neighbourhood()
    f_star = ones(2) .* Inf
    f_evals = 0
    t = copy(start_temp)
    T = ATAmodel.settings.T
    n_items = ATAmodel.settings.n_items
    n_FS = ATAmodel.settings.n_FS
    FS_counts = ATAmodel.settings.FS.counts * ones(Float64, T)'
    #Fill up
    println("Fill-up starting...")
    round = 1

    if n_fill > 0
        for round = 1:n_fill
            NH₁ = _mycopy(NH_start, NH₁)
            warmup = fill(true, T)
            while any(warmup) #filling forms
                v = findfirst(warmup .== true)
                constraints = ATAmodel.constraints[v]
                #println("ol = ", NH₁.ol)
                fₜ = (1 - opt_feas) * (NH₁.infeas + NH₁.ol) - (opt_feas * NH₁.obj)
                mm = fₜ[v]
                for v2 in findall(warmup .== true)
                    if fₜ[v2] > mm
                        mm = copy(fₜ[v2])
                        v = copy(v2)
                    end
                end
                #filling
                n_t = LinearAlgebra.dot(NH₁.x[:, v], ATAmodel.settings.FS.counts)
                #println("test ", v, " chosen, length was: ", n_t)
                #try to add other items, the first time it goes over n_max it stops
                if n_t < constraints.length_max
                    if opt_feas == 0 || fF == true
                        NH_add = fill_up_feas(
                            NH₁,
                            v,
                            ATAmodel.settings.IU,
                            ATAmodel.settings.FS,
                            constraints,
                            ATAmodel.settings.forced0[v],
                            n_items,
                            n_FS,
                            ATAmodel.settings.ol_max[:, v],
                        )
                        NH_add.f = opt_feas * NH_add.f
                    elseif ATAmodel.obj.type == "MAXIMIN"
                        NH_add = fill_up_MAXIMIN(
                            NH₁,
                            IIF[v],
                            opt_feas,
                            v,
                            ATAmodel.settings.IU,
                            ATAmodel.settings.FS,
                            constraints,
                            ATAmodel.settings.forced0[v],
                            n_items,
                            n_FS,
                            ATAmodel.settings.ol_max[:, v],
                        )
                    elseif ATAmodel.obj.type == "MINIMAX"
                        NH_add = fill_up_MINIMAX(
                            NH₁,
                            IIF[v],
                            opt_feas,
                            v,
                            ATAmodel.settings.IU,
                            ATAmodel.settings.FS,
                            constraints,
                            ATAmodel.settings.forced0[v],
                            n_items,
                            n_FS,
                            ATAmodel.settings.ol_max[:, v],
                        )
                    end
                    NH_add.ol = eval_overlap(
                        NH_add.x,
                        FS_counts,
                        ATAmodel.settings.ol_max,
                        T,
                        NH_add.ol,
                    )
                    Printf.@printf "."
                    nₜ = LinearAlgebra.dot(NH_add.x[:, v], ATAmodel.settings.FS.counts)
                    #println("length for test ", v, ": ", nₜ)
                    if n_t <= constraints.length_max
                        NH₁ = _mycopy(NH_add, NH₁)
                    else
                        warmup[v] = false
                        #NH₁ = _mycopy(NH_add, NH₁)
                        println("f₁:", NH₁.f)
                        println("-Test ", v, " filled up with ", n_t, " items, ")
                    end
                else
                    warmup[v] = false
                    println("-Test ", v, " filled up with ", n_t, " items, ")
                end
            end
        end #end of round, round = nRound
        NH₀ = _mycopy(NH₁, NH₀)
        NH₀.f = _comp_f(NH₀, opt_feas)
    end
    NH⁺ = _mycopy(NH₀, NH⁺)
    println("End of fill-up")
    if sum(NH₀.infeas + NH₀.ol) + NH₀.iu == 0
        println("Feasible solution found in fill-up")
        fF = false
    else
        println("Feasible solution not found in fill-up")
    end
    println("f = ", NH⁺.f)
    println("infeas = ", NH⁺.infeas)
    println("ol = ", NH⁺.ol)
    println("iu = ", NH⁺.iu)
    coverage_ok = 0
    if n_test_sample > T
        n_test_sample = T
    end
    convergence = 0

    while coverage_ok == 0
        NH₁ = _mycopy(NH₀, NH₁)
        weights = (1 - opt_feas) .* (NH₀.infeas + NH₀.ol) - opt_feas .* (NH₀.obj)
        weights = (weights .- minimum(weights)) .+ 1
        weights = weights ./ sum(weights)
        weights = StatsBase.ProbabilityWeights(weights)
        #determine test order
        #testOrder = StatsBase.sample(collect(1:T), weights, n_test_sample, replace = false)
        #testOrder = Random.shuffle!(collect(1:ATAmodel.settings.T))
        testOrder = sortperm(weights, rev = true)
        #println("test order = ", Int.(testOrder))
        v₂ = 0
        exit = 0
        #iteratorTestItem = vec(collect(Iterators.product(collect(1:n_item_sample), testOrder[1:n_test_sample])))
        #println(iteratorTestItem[1:nItemoStatsBase.sample])
        #it = 0
        #xnew = copy(NH₀.x)
        while exit == 0 && v₂ < n_test_sample #it<size(iteratorTestItem, 1) #
            #it+= 1
            v₂ += 1
            exit = 0
            v = testOrder[v₂]
            x_forced0ₜ = ATAmodel.settings.forced0[v]
            Constraintsₜ = ATAmodel.constraints[v]
            IIFₜ = IIF[v]
            ol_maxₜ = ATAmodel.settings.ol_max[:, v]
            #it<size(iteratorTestItem, 1) #
            #v = iteratorTestItem[it][2]
            #NH₀.x = copy(xnew)
            taken_items = findall(NH₀.x[:, v] .== 1)
            if n_item_sample > size(taken_items, 1)
                nI = Int(size(taken_items, 1))
            else
                nI = n_item_sample
            end
            # #Random.shuffle!(taken_items) #reset
            # exit2 = 0
            # if iteratorTestItem[it][1]>size(taken_items, 1)
            # 	exit2 = 1
            # end
            #if exit2 == 0
            # if v!= v2
            # 	v = copy(v2)
            #
            # end
            #h = taken_items[iteratorTestItem[it][1]]
            #println("test ", v₂, " of ", size(testOrder, 1))
            #fix test features
            h₂ = 0
            while exit == 0 && h₂ < nI
                NH₁ = _mycopy(NH₀, NH₁)
                h₂ += 1
                #println("item ", h₂, " of ", size(taken_items, 1))
                #try to remove h
                h = taken_items[h₂]
                #iu = copy(iu₀)
                #if rand()>0.0
                NH₁.x[h, v] = zero(Float64)
                #iu[h]-= 1
                taken = 0
                #else
                #	taken = 1
                #end
                if sum(NH₁.x[:, v] .* ATAmodel.settings.FS.counts) >=
                   Constraintsₜ.length_min
                    NH₁.infeas[v], x_Iₜ = check_feas(
                        ATAmodel.settings.FS,
                        Constraintsₜ,
                        NH₁.x[:, v],
                        n_FS,
                        n_items,
                    )
                    iu = sum(NH₁.x, dims = 2) - ATAmodel.settings.IU.max
                    iu = iu[iu.>0]
                    if size(iu, 1) == 0
                        NH₁.iu = 0
                    else
                        NH₁.iu = sum(iu)
                    end
                    NH₁.ol =
                        eval_overlap(NH₁.x, FS_counts, ATAmodel.settings.ol_max, T, NH₁.ol)
                    #NH₁.ol[v] = eval_overlap_v(NH₁.x[:, v], NH₁.x, ATAmodel.settings.FS.counts, ol_maxₜ, v)
                    if fF == false
                        NH₁.obj[v] = eval_TIF_MMₜ(x_Iₜ, IIFₜ)
                    end
                    NH₁.f = _comp_f(NH₁, opt_feas)
                    f_evals += 1
                    if (NH₁.f <= NH₀.f)
                        #switch item
                        #NH₀.ol = eval_overlap(NH₀.x, FS_counts, ATAmodel.settings.ol_max, T, NH₀.ol)
                        #NH₀.f = _comp_f(NH₀, opt_feas)
                        NH₀ = _mycopy(NH₁, NH₀)
                        # println("better to remove, new f₀ = ", NH₀.f)
                        if (NH₁.f < NH⁺.f)
                            exit = 1
                            convergence = 0
                            NH⁺ = _mycopy(NH₁, NH⁺)
                            if verbosity == 2
                                println("f⁺ = ", NH⁺.f, ", t = ", v)
                            else
                                Printf.@printf "+"
                            end
                        end
                    else
                        p = _sig_c(-(NH₁.f - NH₀.f) / t) # -(NH₁.f - NH₀.f) / t always negative
                        if p < 1e-10
                            p = zeros(Float64)[1]
                        end
                        if (rand() < p)
                            #remove
                            #exit = 1
                            NH₀ = _mycopy(NH₁, NH₀)
                            #NH₀.ol = eval_overlap(NH₀.x, FS_counts, ATAmodel.settings.ol_max, T, NH₀.ol)
                            #NH₀.f = _comp_f(NH₀, opt_feas)
                            if verbosity == 2
                                println("SA: f₀ = ", NH₀.f, ", t = ", v)
                            else
                                Printf.@printf "_"
                            end
                        end
                    end
                end
                #try to switch
                idxₜ₂ = findall(NH₁.x[:, v] .== zero(Float64)) #Random.shuffle!(findall(NH₁.x[:, v] .== zero(Float64)))#i₂ = 1, ..., I
                betterFound = 0
                i₃ = 0
                x₋₁ = copy(NH₁.x)
                while betterFound == 0 && i₃ < size(idxₜ₂, 1)
                    i₃ += 1
                    i₂ = idxₜ₂[i₃]
                    if x_forced0ₜ[i₂] && i₂ != h
                        NH₁ = _mycopy(NH₀, NH₁)
                        NH₁.x = copy(x₋₁)
                        #NH₁.x[h, v] = zero(Float64)
                        NH₁.x[i₂, v] = one(Float64)
                        # iu = copy(iu₀)
                        # iu[i₂]+= 1
                        # iu[h]-= 1
                        iu = sum(NH₁.x, dims = 2) - ATAmodel.settings.IU.max
                        iu = iu[iu.>0]
                        if size(iu, 1) == 0
                            NH₁.iu = 0
                        else
                            NH₁.iu = sum(iu)
                        end
                        NH₁.infeas[v], x_Iₜ = check_feas(
                            ATAmodel.settings.FS,
                            Constraintsₜ,
                            NH₁.x[:, v],
                            n_FS,
                            n_items,
                        )
                        if fF == false
                            NH₁.obj[v] = eval_TIF_MMₜ(x_Iₜ, IIFₜ)
                        end
                        NH₁.ol = eval_overlap(
                            NH₁.x,
                            FS_counts,
                            ATAmodel.settings.ol_max,
                            T,
                            NH₁.ol,
                        )
                        #NH₁.ol[v] = eval_overlap_v(NH₁.x[:, v], NH₁.x, ATAmodel.settings.FS.counts, ol_maxₜ, v)
                        NH₁.f = _comp_f(NH₁, opt_feas)
                        if (NH₁.f <= NH₀.f)
                            #switch item
                            #NH₀ = _mycopy(NH₁, NH₀)
                            # NH₀.ol = eval_overlap(NH₀.x, FS_counts, ATAmodel.settings.ol_max, T, NH₀.ol)
                            # NH₀.f = _comp_f(NH₀, opt_feas)
                            NH₀ = _mycopy(NH₁, NH₀)
                            # println("better to switch, new f₀ = ", NH₀.f)
                            if (NH₁.f < NH⁺.f)
                                exit = 1
                                betterFound = 1
                                convergence = 0
                                NH⁺ = _mycopy(NH₁, NH⁺)
                                if verbosity == 2
                                    println("f⁺ = ", NH⁺.f, ", t = ", v)
                                else
                                    Printf.@printf "+"
                                end
                            end
                        else
                            p = _sig_c(-(NH₁.f - NH₀.f) / t)
                            if p < 1e-10
                                p = zeros(Float64)[1]
                            end
                            if (rand() < p)
                                #remove
                                NH₀ = _mycopy(NH₁, NH₀)
                                # NH₀.ol = eval_overlap(NH₀.x, FS_counts, ATAmodel.settings.ol_max, T, NH₀.ol)
                                # NH₀.f = _comp_f(NH₀, opt_feas)
                                if verbosity == 2
                                    println("SA: f₀ = ", NH₀.f, ", t = ", v)
                                else
                                    Printf.@printf "_"
                                end
                            end
                        end
                    end
                end #end of betterFound (betterFound = 1)
            end #end of itemorder h₂ (exit = 1)

        end #end of testorder v₂ (exit = 1)
        if sum(NH₀.infeas + NH₀.ol) + NH₀.iu <= 0 && fF == true
            fF = false
            for v = 1:T
                x_Iₜ = FS_to_items(
                    NH₀.x[:, v],
                    ATAmodel.settings.n_items,
                    ATAmodel.settings.FS.items,
                )
                NH₀.obj[v] = eval_TIF_MMₜ(x_Iₜ, IIF[v])
            end
            NH₀.f = _comp_f(NH₀, opt_feas)
            NH⁺ = _mycopy(NH₀, NH⁺)
        end
        f_star[1] = copy(NH₀.f)
        #if exit == 0
        if f_star[2] == f_star[1]
            convergence += 1
            println(convergence)
        end
        #println("convergence is ", convergence)
        #how many equal f₀ in the last iterations?
        if (f_star[2] == f_star[1] && convergence == max_conv) ||
           time() - start_time >= max_time
            coverage_ok = 1
            #t += 1
        else
            t *= geom_temp
            coverage_ok = 0
            if f_star[1] <= f_star[2]
                f_star[2] = copy(f_star[1])
            end
        end
        if verbosity == 2
            println("\n")
            Printf.@printf("\n  f⁺:	%16.10f", NH⁺.f)
            Printf.@printf("\n  f₀:	%16.10f", NH₀.f)
            Printf.@printf("\n  Convergence:	%16.10f", convergence)
            println("\n")
        end
    end #end of NH coverage (coverage_ok = 1)
    return NH⁺::Neighbourhood
end

#custom neighborhood
function analyse_NH(
    NH_start::Neighbourhood,
    ATAmodel::Model,
    IIF::Vector{Float64};
    fF = true,
    n_fill = 1,
    opt_feas = 0.9,
    max_conv = 1,
    start_temp = 1000.0,
    geom_temp = 0.1,
    n_item_sample = 1,
    n_test_sample = 1,
    verbosity = 1,
    start_time = 0,
    max_time = 1000,
)
    if fF == true
        NH_start.obj = zeros(Float64, ATAmodel.settings.T)
    end
    NH₁ = Neighbourhood()
    NH₁ = _mycopy(NH_start, NH₁)
    NH₀ = Neighbourhood()
    NH⁺ = Neighbourhood()
    f_star = ones(2) .* Inf
    f_evals = 0
    t = copy(start_temp)
    T = ATAmodel.settings.T
    n_items = ATAmodel.settings.n_items
    n_FS = ATAmodel.settings.n_FS
    FS_counts = ATAmodel.settings.FS.counts * ones(Float64, T)'
    #Fill up
    println("Fill-up starting...")
    round = 1

    if n_fill > 0
        for round = 1:n_fill
            NH₁ = _mycopy(NH_start, NH₁)
            warmup = fill(true, T)
            while any(warmup) #filling forms
                v = findfirst(warmup .== true)
                constraints = ATAmodel.constraints[v]
                #println("ol = ", NH₁.ol)
                fₜ = (1 - opt_feas) * (NH₁.infeas + NH₁.ol) - (opt_feas * NH₁.obj)
                mm = fₜ[v]
                for v2 in findall(warmup .== true)
                    if fₜ[v2] > mm
                        mm = copy(fₜ[v2])
                        v = copy(v2)
                    end
                end
                #filling
                n_t = LinearAlgebra.dot(NH₁.x[:, v], ATAmodel.settings.FS.counts)
                #println("test ", v, " chosen, length was: ", n_t)
                #try to add other items, the first time it goes over n_max it stops
                if n_t < constraints.length_max
                    if opt_feas == 0 || fF == true
                        NH_add = fill_up_feas(
                            NH₁,
                            v,
                            ATAmodel.settings.IU,
                            ATAmodel.settings.FS,
                            constraints,
                            ATAmodel.settings.forced0[v],
                            n_items,
                            n_FS,
                            ATAmodel.settings.ol_max[:, v],
                        )
                        NH_add.f = opt_feas * NH_add.f
                    else
                        NH_add = fill_up_custom(
                            NH₁,
                            ATAmodel.obj,
                            opt_feas,
                            v,
                            ATAmodel.settings.IU,
                            ATAmodel.settings.FS,
                            constraints,
                            ATAmodel.settings.forced0[v],
                            n_items,
                            n_FS,
                            ATAmodel.settings.ol_max[:, v],
                        )
                    end
                    NH_add.ol = eval_overlap(
                        NH_add.x,
                        FS_counts,
                        ATAmodel.settings.ol_max,
                        T,
                        NH_add.ol,
                    )
                    Printf.@printf "."
                    nₜ = LinearAlgebra.dot(NH_add.x[:, v], ATAmodel.settings.FS.counts)
                    #println("length for test ", v, ": ", nₜ)
                    if n_t <= constraints.length_max
                        NH₁ = _mycopy(NH_add, NH₁)
                    else
                        warmup[v] = false
                        #NH₁ = _mycopy(NH_add, NH₁)
                        println("f₁:", NH₁.f)
                        println("-Test ", v, " filled up with ", n_t, " items, ")
                    end
                else
                    warmup[v] = false
                    println("-Test ", v, " filled up with ", n_t, " items, ")
                end
            end
        end #end of round, round = nRound
        NH₀ = _mycopy(NH₁, NH₀)
        NH₀.f = _comp_f(NH₀, opt_feas)
    end
    NH⁺ = _mycopy(NH₀, NH⁺)
    println("End of fill-up")
    if sum(NH₀.infeas + NH₀.ol) + NH₀.iu == 0
        println("Feasible solution found in fill-up")
        fF = false
    else
        println("Feasible solution not found in fill-up")
    end
    println("f = ", NH⁺.f)
    println("infeas = ", NH⁺.infeas)
    println("ol = ", NH⁺.ol)
    println("iu = ", NH⁺.iu)
    coverage_ok = 0
    if n_test_sample > T
        n_test_sample = T
    end
    convergence = 0

    while coverage_ok == 0
        NH₁ = _mycopy(NH₀, NH₁)
        weights = (1 - opt_feas) .* (NH₀.infeas + NH₀.ol) - opt_feas .* (NH₀.obj)
        weights = (weights .- minimum(weights)) .+ 1
        weights = weights ./ sum(weights)
        weights = StatsBase.ProbabilityWeights(weights)
        #determine test order
        #testOrder = StatsBase.sample(collect(1:T), weights, n_test_sample, replace = false)
        #testOrder = Random.shuffle!(collect(1:ATAmodel.settings.T))
        testOrder = sortperm(weights, rev = true)
        #println("test order = ", Int.(testOrder))
        v₂ = 0
        exit = 0
        #iteratorTestItem = vec(collect(Iterators.product(collect(1:n_item_sample), testOrder[1:n_test_sample])))
        #println(iteratorTestItem[1:nItemoStatsBase.sample])
        #it = 0
        #xnew = copy(NH₀.x)
        while exit == 0 && v₂ < n_test_sample #it<size(iteratorTestItem, 1) #
            #it+= 1
            v₂ += 1
            exit = 0
            v = testOrder[v₂]
            x_forced0ₜ = ATAmodel.settings.forced0[v]
            Constraintsₜ = ATAmodel.constraints[v]
            IIFₜ = IIF[v]
            ol_maxₜ = ATAmodel.settings.ol_max[:, v]
            #it<size(iteratorTestItem, 1) #
            #v = iteratorTestItem[it][2]
            #NH₀.x = copy(xnew)
            taken_items = findall(NH₀.x[:, v] .== 1)
            if n_item_sample > size(taken_items, 1)
                nI = Int(size(taken_items, 1))
            else
                nI = n_item_sample
            end
            # Random.shuffle!(taken_items) #reset
            # exit2 = 0
            # if iteratorTestItem[it][1]>size(taken_items, 1)
            # 	exit2 = 1
            # end
            #if exit2 == 0
            # if v!= v2
            # 	v = copy(v2)
            #
            # end
            #h = taken_items[iteratorTestItem[it][1]]
            #println("test ", v₂, " of ", size(testOrder, 1))
            #fix test features
            h₂ = 0
            while exit == 0 && h₂ < nI
                NH₁ = _mycopy(NH₀, NH₁)
                h₂ += 1
                #println("item ", h₂, " of ", size(taken_items, 1))
                #try to remove h
                h = taken_items[h₂]
                #iu = copy(iu₀)
                #if rand()>0.0
                NH₁.x[h, v] = zero(Float64)
                #iu[h]-= 1
                taken = 0
                #else
                #	taken = 1
                #end
                if sum(NH₁.x[:, v] .* ATAmodel.settings.FS.counts) >=
                   Constraintsₜ.length_min
                    NH₁.infeas[v], x_Iₜ = check_feas(
                        ATAmodel.settings.FS,
                        Constraintsₜ,
                        NH₁.x[:, v],
                        n_FS,
                        n_items,
                    )
                    iu = sum(NH₁.x, dims = 2) - ATAmodel.settings.IU.max
                    iu = iu[iu.>0]
                    if size(iu, 1) == 0
                        NH₁.iu = 0
                    else
                        NH₁.iu = sum(iu)
                    end
                    NH₁.ol =
                        eval_overlap(NH₁.x, FS_counts, ATAmodel.settings.ol_max, T, NH₁.ol)
                    #NH₁.ol[v] = eval_overlap_v(NH₁.x[:, v], NH₁.x, ATAmodel.settings.FS.counts, ol_maxₜ, v)
                    if fF == false
                        NH₁.obj[v] = eval_TIF_MMₜ(x_Iₜ, IIFₜ)
                    end
                    NH₁.f = _comp_f(NH₁, opt_feas)
                    f_evals += 1
                    if (NH₁.f <= NH₀.f)
                        #switch item
                        #NH₀.ol = eval_overlap(NH₀.x, FS_counts, ATAmodel.settings.ol_max, T, NH₀.ol)
                        #NH₀.f = _comp_f(NH₀, opt_feas)
                        NH₀ = _mycopy(NH₁, NH₀)
                        # println("better to remove, new f₀ = ", NH₀.f)
                        if (NH₁.f < NH⁺.f)
                            exit = 1
                            convergence = 0
                            NH⁺ = _mycopy(NH₁, NH⁺)
                            if verbosity == 2
                                println("f⁺ = ", NH⁺.f, ", t = ", v)
                            else
                                Printf.@printf "+"
                            end
                        end
                    else
                        p = _sig_c(-(NH₁.f - NH₀.f) / t) # -(NH₁.f - NH₀.f) / t always negative
                        if p < 1e-10
                            p = zeros(Float64)[1]
                        end
                        if (rand() < p)
                            #remove
                            #exit = 1
                            NH₀ = _mycopy(NH₁, NH₀)
                            #NH₀.ol = eval_overlap(NH₀.x, FS_counts, ATAmodel.settings.ol_max, T, NH₀.ol)
                            #NH₀.f = _comp_f(NH₀, opt_feas)
                            if verbosity == 2
                                println("SA: f₀ = ", NH₀.f, ", t = ", v)
                            else
                                Printf.@printf "_"
                            end
                        end
                    end
                end
                #try to switch
                idxₜ₂ = findall(NH₁.x[:, v] .== zero(Float64)) #Random.shuffle!(findall(NH₁.x[:, v] .== zero(Float64)))#i₂ = 1, ..., I
                betterFound = 0
                i₃ = 0
                x₋₁ = copy(NH₁.x)
                while betterFound == 0 && i₃ < size(idxₜ₂, 1)
                    i₃ += 1
                    i₂ = idxₜ₂[i₃]
                    if x_forced0ₜ[i₂] && i₂ != h
                        NH₁ = _mycopy(NH₀, NH₁)
                        NH₁.x = copy(x₋₁)
                        #NH₁.x[h, v] = zero(Float64)
                        NH₁.x[i₂, v] = one(Float64)
                        # iu = copy(iu₀)
                        # iu[i₂]+= 1
                        # iu[h]-= 1
                        iu = sum(NH₁.x, dims = 2) - ATAmodel.settings.IU.max
                        iu = iu[iu.>0]
                        if size(iu, 1) == 0
                            NH₁.iu = 0
                        else
                            NH₁.iu = sum(iu)
                        end
                        NH₁.infeas[v], x_Iₜ = check_feas(
                            ATAmodel.settings.FS,
                            Constraintsₜ,
                            NH₁.x[:, v],
                            n_FS,
                            n_items,
                        )
                        if fF == false
                            NH₁.obj[v] = eval_TIF_MMₜ(x_Iₜ, IIFₜ)
                        end
                        NH₁.ol = eval_overlap(
                            NH₁.x,
                            FS_counts,
                            ATAmodel.settings.ol_max,
                            T,
                            NH₁.ol,
                        )
                        #NH₁.ol[v] = eval_overlap_v(NH₁.x[:, v], NH₁.x, ATAmodel.settings.FS.counts, ol_maxₜ, v)
                        NH₁.f = _comp_f(NH₁, opt_feas)
                        if (NH₁.f <= NH₀.f)
                            #switch item
                            #NH₀ = _mycopy(NH₁, NH₀)
                            # NH₀.ol = eval_overlap(NH₀.x, FS_counts, ATAmodel.settings.ol_max, T, NH₀.ol)
                            # NH₀.f = _comp_f(NH₀, opt_feas)
                            NH₀ = _mycopy(NH₁, NH₀)
                            # println("better to switch, new f₀ = ", NH₀.f)
                            if (NH₁.f < NH⁺.f)
                                exit = 1
                                betterFound = 1
                                convergence = 0
                                NH⁺ = _mycopy(NH₁, NH⁺)
                                if verbosity == 2
                                    println("f⁺ = ", NH⁺.f, ", t = ", v)
                                else
                                    Printf.@printf "+"
                                end
                            end
                        else
                            p = _sig_c(-(NH₁.f - NH₀.f) / t)
                            if p < 1e-10
                                p = zeros(Float64)[1]
                            end
                            if (rand() < p)
                                #remove
                                NH₀ = _mycopy(NH₁, NH₀)
                                # NH₀.ol = eval_overlap(NH₀.x, FS_counts, ATAmodel.settings.ol_max, T, NH₀.ol)
                                # NH₀.f = _comp_f(NH₀, opt_feas)
                                if verbosity == 2
                                    println("SA: f₀ = ", NH₀.f, ", t = ", v)
                                else
                                    Printf.@printf "_"
                                end
                            end
                        end
                    end
                end #end of betterFound (betterFound = 1)
            end #end of itemorder h₂ (exit = 1)

        end #end of testorder v₂ (exit = 1)
        if sum(NH₀.infeas + NH₀.ol) + NH₀.iu <= 0 && fF == true
            fF = false
            for v = 1:T
                x_Iₜ = FS_to_items(
                    NH₀.x[:, v],
                    ATAmodel.settings.n_items,
                    ATAmodel.settings.FS.items,
                )
                NH₀.obj[v] = eval_TIF_MMₜ(x_Iₜ, IIF[v])
            end
            NH₀.f = _comp_f(NH₀, opt_feas)
            NH⁺ = _mycopy(NH₀, NH⁺)
        end
        f_star[1] = copy(NH₀.f)
        #if exit == 0
        if f_star[2] == f_star[1]
            convergence += 1
            println(convergence)
        end
        #println("convergence is ", convergence)
        #how many equal f₀ in the last iterations?
        if (f_star[2] == f_star[1] && convergence == max_conv) ||
           time() - start_time >= max_time
            coverage_ok = 1
            #t += 1
        else
            t *= geom_temp
            coverage_ok = 0
            if f_star[1] <= f_star[2]
                f_star[2] = copy(f_star[1])
            end
        end
        if verbosity == 2
            println("\n")
            Printf.@printf("\n  f⁺:	%16.10f", NH⁺.f)
            Printf.@printf("\n  f₀:	%16.10f", NH₀.f)
            Printf.@printf("\n  Convergence:	%16.10f", convergence)
            println("\n")
        end
    end #end of NH coverage (coverage_ok = 1)
    return NH⁺::Neighbourhood
end
