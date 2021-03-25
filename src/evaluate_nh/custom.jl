
#custom neighborhood
function analyse_NH(
    NH_start::Neighbourhood,
    ata_model::CustomModel;
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
        NH_start.obj = zeros(Float64, ata_model.settings.T)
    end
    NH₁ = Neighbourhood()
    NH₀ = Neighbourhood()
    NH⁺ = Neighbourhood()
    NH₁ = _mycopy(NH_start, NH₁)
    NH₀ = _mycopy(NH_start, NH₀)
    NH⁺ = _mycopy(NH_start, NH⁺)

    f_star = ones(2) .* Inf
    f_evals = 0
    t = copy(start_temp)
    T = ata_model.settings.T
    n_items = ata_model.settings.n_items
    n_fs = ata_model.settings.n_fs
    fs_counts = ata_model.settings.fs.counts * ones(Float64, T)'

    if n_fill > 0
        #Fill up
        println("Fill-up starting...")
        round = 1
        for round = 1:n_fill
            NH₁ = _mycopy(NH_start, NH₁)
            warmup = fill(true, T)
            while any(warmup) #filling forms
                v = findfirst(warmup .== true)
                constraints = ata_model.constraints[v]
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
                n_t = LinearAlgebra.dot(NH₁.x[:, v], ata_model.settings.fs.counts)
                #println("test ", v, " chosen, length was: ", n_t)
                #try to add other items, the first time it goes over n_max it stops
                if n_t < constraints.length_max
                    if opt_feas == 0 || fF == true
                        NH_add = find_best_itemᵥ(
                            NH₁,
                            v,
                            ata_model.settings.iu,
                            ata_model.settings.fs,
                            constraints,
                            ata_model.settings.forced0[v],
                            n_items,
                            n_fs,
                            ata_model.settings.ol_max[:, v],
                            ata_model.settings.to_apply,
                        )
                        NH_add.f = opt_feas * NH_add.f
                    else
                        NH_add = find_best_itemᵥ(
                            NH₁,
                            ata_model.obj,
                            opt_feas,
                            v,
                            ata_model.settings.iu,
                            ata_model.settings.fs,
                            constraints,
                            ata_model.settings.forced0[v],
                            n_items,
                            n_fs,
                            ata_model.settings.ol_max[:, v],
                            ata_model.settings.to_apply,
                        )
                    end
                    if ata_model.settings.to_apply[3]
                        NH_add.ol = eval_overlap(
                            NH_add.x,
                            fs_counts,
                            ata_model.settings.ol_max,
                            T,
                            NH_add.ol,
                        )
                    end
                    Printf.@printf "."
                    n_t = LinearAlgebra.dot(NH_add.x[:, v], ata_model.settings.fs.counts)
                    #println("length for test ", v, ": ", n_t)
                    if n_t <= constraints.length_max
                        NH₁ = _mycopy(NH_add, NH₁)
                    else
                        warmup[v] = false
                        #NH₁ = _mycopy(NH_add, NH₁)
                        Printf.@printf("\n f₁:  ", NH₁.f)
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
        NH⁺ = _mycopy(NH₀, NH⁺)
        println("End of fill-up")
        if sum(NH₀.infeas + NH₀.ol) + NH₀.iu == 0
            println("Feasible solution found in fill-up")
            fF = false
        else
            println("Feasible solution not found in fill-up")
        end
    end

    if verbosity > 1
        print_neighbourhood(NH⁺)
    end
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
        #Random.shuffle!(testOrder)
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
            x_forced0ₜ = ata_model.settings.forced0[v]
            Constraintsₜ = ata_model.constraints[v]
            ol_maxₜ = ata_model.settings.ol_max[:, v]
            #it<size(iteratorTestItem, 1) #
            #v = iteratorTestItem[it][2]
            #NH₀.x = copy(xnew)
            taken_items = findall(isone.(NH₀.x[:, v]))
            if n_item_sample > size(taken_items, 1)
                nI = Int(size(taken_items, 1))
            else
                nI = n_item_sample
            end
            taken_items = Random.shuffle!(taken_items) # !removed
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
                if sum(NH₁.x[:, v] .* ata_model.settings.fs.counts) >=
                   Constraintsₜ.length_min
                    NH₁.infeas[v], x_Iₜ = check_feas(
                        ata_model.settings.fs,
                        Constraintsₜ,
                        NH₁.x[:, v],
                        n_fs,
                        n_items,
                    )
                    if ata_model.settings.to_apply[1]
                        iu_max = sum(NH₁.x, dims = 2)[:, 1] - ata_model.settings.iu.max
                        iu_max = sum_pos(iu_max)
                    else
                        iu_max = 0
                    end
                    if ata_model.settings.to_apply[2]
                        iu_min = -sum(NH₁.x, dims = 2)[:, 1] + ata_model.settings.iu.min
                        iu_min = sum_pos(iu_min)
                    else
                        iu_min = 0
                    end
                    NH₁.iu = iu_max + iu_min
                    if ata_model.settings.to_apply[3]
                        NH₁.ol = eval_overlap(
                            NH₁.x,
                            fs_counts,
                            ata_model.settings.ol_max,
                            T,
                            NH₁.ol,
                        )
                    end
                    #NH₁.ol[v] = eval_overlapᵥ(NH₁.x[:, v], NH₁.x, ata_model.settings.fs.counts, ol_maxₜ, v)
                    if fF == false
                        NH₁.obj = ata_model.obj.fun(NH₁.x, ata_model.obj.args)
                    end
                    NH₁.f = _comp_f(NH₁, opt_feas)
                    f_evals += 1
                    if (NH₁.f <= NH₀.f)
                        #switch item
                        #NH₀.ol = eval_overlap(NH₀.x, fs_counts, ata_model.settings.ol_max, T, NH₀.ol)
                        #NH₀.f = _comp_f(NH₀, opt_feas)
                        NH₀ = _mycopy(NH₁, NH₀)
                        # println("better to remove, new f₀ = ", NH₀.f)
                        if (NH₁.f < NH⁺.f)
                            exit = 1
                            convergence = 0
                            NH⁺ = _mycopy(NH₁, NH⁺)
                            if verbosity == 2
                                Printf.@printf("\n f⁺ : %5.3f", NH⁺.f)
                                Printf.@printf(" in test : %5d", v)
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
                            #NH₀.ol = eval_overlap(NH₀.x, fs_counts, ata_model.settings.ol_max, T, NH₀.ol)
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
                idxₜ₂ = findall(iszero.(NH₁.x[:, v]))
                idxₜ₂ = Random.shuffle!(idxₜ₂)#i₂ = 1, ..., I !removed
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
                        if ata_model.settings.to_apply[1]
                            iu_max = sum(NH₁.x, dims = 2)[:, 1] - ata_model.settings.iu.max
                            iu_max = sum_pos(iu_max)
                        else
                            iu_max = 0
                        end
                        if ata_model.settings.to_apply[2]
                            iu_min = -sum(NH₁.x, dims = 2)[:, 1] + ata_model.settings.iu.min
                            iu_min = sum_pos(iu_min)
                        else
                            iu_min = 0
                        end
                        NH₁.iu = iu_max + iu_min
                        NH₁.infeas[v], x_Iₜ = check_feas(
                            ata_model.settings.fs,
                            Constraintsₜ,
                            NH₁.x[:, v],
                            n_fs,
                            n_items,
                        )
                        if fF == false
                            NH₁.obj = ata_model.obj.fun(NH₁.x, ata_model.obj.args)
                        end
                        if ata_model.settings.to_apply[3]
                            NH₁.ol = eval_overlap(
                                NH₁.x,
                                fs_counts,
                                ata_model.settings.ol_max,
                                T,
                                NH₁.ol,
                            )
                        end
                        #NH₁.ol[v] = eval_overlapᵥ(NH₁.x[:, v], NH₁.x, ata_model.settings.fs.counts, ol_maxₜ, v)
                        NH₁.f = _comp_f(NH₁, opt_feas)
                        if (NH₁.f <= NH₀.f)
                            #switch item
                            #NH₀ = _mycopy(NH₁, NH₀)
                            # NH₀.ol = eval_overlap(NH₀.x, fs_counts, ata_model.settings.ol_max, T, NH₀.ol)
                            # NH₀.f = _comp_f(NH₀, opt_feas)
                            NH₀ = _mycopy(NH₁, NH₀)
                            # println("better to switch, new f₀ = ", NH₀.f)
                            if (NH₁.f < NH⁺.f)
                                exit = 1
                                betterFound = 1
                                convergence = 0
                                NH⁺ = _mycopy(NH₁, NH⁺)
                                if verbosity == 2
                                    Printf.@printf("\n f⁺ : %5.3f", NH⁺.f)
                                    Printf.@printf(" in test : %5d", v)
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
                                # NH₀.ol = eval_overlap(NH₀.x, fs_counts, ata_model.settings.ol_max, T, NH₀.ol)
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
                x_Iₜ = fs_to_items(
                    NH₀.x[:, v],
                    ata_model.settings.n_items,
                    ata_model.settings.fs.items,
                )
            end
            NH₀.obj = ata_model.obj.fun(NH₀.x, ata_model.obj.args)
            NH₀.f = _comp_f(NH₀, opt_feas)
            NH⁺ = _mycopy(NH₀, NH⁺)
        end
        f_star[1] = copy(NH₀.f)
        #if exit == 0
        if f_star[2] == f_star[1]
            convergence += 1
            Printf.@printf(" %5d", convergence)
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
            Printf.@printf("\n  f⁺:	%5.3f", NH⁺.f)
            Printf.@printf("\n  f₀:	%5.3f", NH₀.f)
            Printf.@printf("\n  Convergence:	%5d", convergence)
            println("\n")
        end
    end #end of NH coverage (coverage_ok = 1)
    print_neighbourhood(NH⁺)
    return NH⁺::Neighbourhood
end
