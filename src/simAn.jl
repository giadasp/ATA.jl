include("evaluate_nh/eval.jl")

# optimize with simulated annealing
function siman!(
    ata_model::AbstractModel;
    starting_design = Matrix{Float64}(undef, 0, 0),
    results_folder = "RESULTS",
    start_temp = 0.1,
    geom_temp = 0.1,
    max_time = 1000.00,
    max_conv = 2,
    feas_nh = 0,
    opt_nh = 5,
    n_item_sample = 1000,
    n_test_sample = 2,
    opt_feas = 0.0,
    n_fill = 1,
    verbosity = 2,
    kwargs...,
)
    message = ""
    if !(results_folder in readdir())
        mkdir(results_folder)
    else
        message *= string(
            "There is already a folder with this name, files in ",
            results_folder,
            " will be overwritten.\n",
        )
    end

    ##
    nNH = feas_nh + opt_nh
    opt_nh = feas_nh + 1
    if feas_nh == 0
        fF = false
    else
        fF = true
    end
    ata_model.settings.n_fs = size(ata_model.settings.fs.items, 1)
    T = ata_model.settings.T
    n_fs = ata_model.settings.n_fs
    if n_fs == 0
        n_fs = ata_model.settings.n_items
    end
    n_items = ata_model.settings.n_items
    fs_counts = ata_model.settings.fs.counts * ones(Float64, T)'
    iu⁺ = 0
    start_time = copy(time())
    nacc = 0 # total accepted trials
    t = copy(start_temp) # temperature - will initially rise or fall to cover parameter space. Then it will fall
    finish = 0 # convergence indicator 0 (failure), 1 (normal success), or 2 (convergence but near bounds)
    f_evals = 0
    hline = " = "^80
    f_star = typemax(Float64) * ones(2)
    NH = 1

    # starting design check
    if size(starting_design, 1) > 0
        if (
            size(starting_design, 1) != n_fs ||
            size(starting_design, 2) != ata_model.settings.T
        )
            push!(
                ata_model.output.infos,
                ["danger", "- Starting design must be of size: (n_items x T).\n"],
            )
            return nothing
        end
        if (any(starting_design != 0 && starting_design != 1))
            push!(
                ata_model.output.infos,
                ["danger", "- Starting design must contain only 1 or 0.\n"],
            )
            return nothing
        end
        x₀ = Float64.(starting_design)
    else
        x₀ = zeros(Float64, n_fs, ata_model.settings.T)
    end

    NH⁺ = Neighbourhood(x₀, Inf, zeros(T), 1e6 * ones(T), 1e6 * ones(T), zero(Float64))
    bestNH = 1
    if size(x₀, 1) == 0
        NH⁺.x = zeros(Float64, n_fs, T)
    end
    # compute f
    iu = sum(NH⁺.x; dims = 2)[:, 1] - ata_model.settings.iu.max
    NH⁺.iu = sum_pos(iu)
    t = copy(start_temp)
    for v2 = 1:T
        NH⁺.infeas[v2], x_Iₜ = check_feas(
            ata_model.settings.fs,
            ata_model.constraints[v2],
            NH⁺.x[:, v2],
            n_fs,
            n_items,
        )
        NH⁺.ol = eval_overlap(NH⁺.x, fs_counts, ata_model.settings.ol_max, T, NH⁺.ol)
        if fF == false
            if ata_model.obj.name in ["MAXIMIN", "MINIMAX", "CCMAXIMIN"]
                NH⁺.obj[v2] = eval_TIFₜ(x_Iₜ, ata_model.obj.cores[v2])
            elseif ata_model.obj.name == "custom"
                NH⁺.obj = ata_model.obj.fun(x_Iₜ, ata_model.obj.args)
            end
        end
    end
    NH⁺.f = _comp_f(NH⁺, opt_feas)
    if sum(NH⁺.x) <= 1.0
        f⁺ = Inf
    else
        f⁺ = copy(NH⁺.f)
    end
    println("Starting solution: ")
    print_neighbourhood(NH⁺)
    if (Distributed.nprocs() > 1)
        workers = collect(1:(Distributed.nprocs()-1))
    else
        workers = [one(Int64)]
    end
    NHs = copy(workers)
    ata_model.output.neighbourhoods = [Neighbourhood() for n = 1:NHs[end]]
    for nh = 1:NHs[end]
        ata_model.output.neighbourhoods[nh] =
            _mycopy(NH⁺, ata_model.output.neighbourhoods[nh])
    end
    round = 1
    while finish == 0
        #for each proc in nprocs analyse nh
        Distributed.@sync Distributed.@distributed for p in workers
            nh_tot = (NHs[end] * (round - 1)) + NHs[p]
            # println("f was ", ata_model.output.neighbourhoods[NHs[p]].f)
            NH_proc = Neighbourhood()
            NH_proc = _mycopy(ata_model.output.neighbourhoods[NHs[p]], NH_proc)
            if round > 1
                #fill-up phase just in first round (first n_procs neighbourhoods)
                n_fill = 0
            end
            NH_proc = analyse_NH(
                NH_proc,
                ata_model;
                fF = fF,
                n_fill = n_fill,
                opt_feas = opt_feas,
                max_conv = max_conv,
                start_temp = t,
                geom_temp = geom_temp,
                n_item_sample = n_item_sample,
                n_test_sample = n_test_sample,
                verbosity = verbosity,
                start_time = start_time,
                max_time = max_time,
            )
            open(string(results_folder, "/neigh_", nh_tot, "_x.csv"), "w") do io
                return DelimitedFiles.writedlm(io, NH_proc.x)
            end
            open(string(results_folder, "/neigh_", nh_tot, "_infeas.csv"), "w") do io
                return DelimitedFiles.writedlm(io, NH_proc.infeas)
            end
            open(string(results_folder, "/neigh_", nh_tot, "_obj.csv"), "w") do io
                return DelimitedFiles.writedlm(io, NH_proc.obj)
            end
            open(string(results_folder, "/neigh_", nh_tot, "_ol.csv"), "w") do io
                return DelimitedFiles.writedlm(io, NH_proc.ol)
            end
            open(string(results_folder, "/neigh_", nh_tot, "_iu.csv"), "w") do io
                return DelimitedFiles.writedlm(io, NH_proc.iu)
            end
            open(string(results_folder, "/neigh_", nh_tot, "_f.csv"), "w") do io
                return DelimitedFiles.writedlm(io, NH_proc.f)
            end
            println("\nNeighbourhood ", nh_tot, " fully explored, increase temperature")
        end
        #save last nprocs nhs in text files
        for p in workers
            nh_tot = (NHs[end] * (round - 1)) + NHs[p]
            ata_model.output.neighbourhoods[NHs[p]].x = DelimitedFiles.readdlm(
                string(results_folder, "/neigh_", nh_tot, "_x.csv"),
                '\t',
            )
            ata_model.output.neighbourhoods[NHs[p]].obj = DelimitedFiles.readdlm(
                string(results_folder, "/neigh_", nh_tot, "_obj.csv"),
                '\t',
            )[
                :,
                1,
            ]
            ata_model.output.neighbourhoods[NHs[p]].infeas = DelimitedFiles.readdlm(
                string(results_folder, "/neigh_", nh_tot, "_infeas.csv"),
                '\t',
            )[
                :,
                1,
            ]
            ata_model.output.neighbourhoods[NHs[p]].ol = DelimitedFiles.readdlm(
                string(results_folder, "/neigh_", nh_tot, "_ol.csv"),
                '\t',
            )[
                :,
                1,
            ]
            ata_model.output.neighbourhoods[NHs[p]].iu = DelimitedFiles.readdlm(
                string(results_folder, "/neigh_", nh_tot, "_iu.csv"),
                '\t',
            )[
                1,
                1,
            ]
            ata_model.output.neighbourhoods[NHs[p]].f = DelimitedFiles.readdlm(
                string(results_folder, "/neigh_", nh_tot, "_f.csv"),
                '\t',
            )[
                1,
                1,
            ]
        end
        # NH+= 1
        # find best nh
        f_loc, nh_loc = findmin([ata_model.output.neighbourhoods[n].f for n = 1:NHs[end]])
        # copy first best result as benchmark (avoid to keep empty design as best solution in large scale models)
        if f_loc <= f⁺
            bestNH = (NHs[end] * (round - 1)) + nh_loc
            f⁺ = copy(f_loc)
            NH⁺ = _mycopy(ata_model.output.neighbourhoods[nh_loc], NH⁺)
        end
        # NHs = NHs.+Distributed.nprocs()
        if NHs[end] * round + 1 > feas_nh
            fF = false
        end
        if NHs[end] * round + 1 > nNH
            finish = 1
        else
            # perturbate the solution and go to warmup
        end
        round += 1
        # f_evals+= f_evals_NH
        # if f_evals >= max_evals
        # 	finish = 2
        # end
        if verbosity >= 1
            println(hline)
            println("Results")
            if (finish == 1)
                println(" == > Maximum number of neighbourhoods explored <== ")
                # elseif (finish == 2) #time max reached
                # 	println(" == > Maximum number of evaluations reached <== ")
            elseif (finish == 3) # evals max reached
                println(" == > Maximum time reached <== ")
            end
            Printf.@printf("\n obj. value:	%5.1f", NH⁺.f)
            Printf.@printf("\n Optimality:	%5.1f", minimum(NH⁺.obj))
            Printf.@printf("\n Infeasibility:	%5.1f", sum(NH⁺.infeas + NH⁺.ol) + NH⁺.iu)
            Printf.@printf("\n Elapsed Time:	%5.1f", time() - start_time)
            Printf.@printf("\n Best Neighbourhood:	%5d", bestNH)
            println("\n")
            println(hline)
        end
        if time() - start_time >= max_time
            ata_model.output.elapsed_time = time() - start_time
            finish = 3
        end
        if finish > 0
            # f⁺, NH = findmin([ata_model.output.neighbourhoods[n].f for n = 1:nNH])
            # x⁺ = ata_model.output.neighbourhoods[NH].x
            # infeas⁺ = ata_model.output.neighbourhoods[NH].infeas+ata_model.output.neighbourhoods[NH].ol
            # iu⁺ = ata_model.output.neighbourhoods[NH].iu
            # TIF⁺ = ata_model.output.neighbourhoods[NH].obj
            ata_model.output.design = NH⁺.x
            ata_model.output.feas = NH⁺.infeas + NH⁺.ol
            ata_model.output.f = NH⁺.f
            ata_model.output.TIF = NH⁺.obj
        end
        # j = 1
        # for v = 1:ata_model.settings.T
        # 	for i = 1:ata_model.settings.n_fs
        # 		if NH⁺.x[i, v] == one(Float64)
        # 			if j == 1
        # 				NH⁺.x[i, v] = zero(Float64)
        # 				j+= 1
        # 			else
        # 				j - = 1
        # 			end
        # 		end
        # 	end
        # end
    end # end of finish
    JLD2.@save string(results_folder, "/ata_model.jld2") ata_model
    DelimitedFiles.writedlm(
        string(results_folder, "/design.csv"),
        reshape(ata_model.output.design, n_fs, T),
    )
    open(string(results_folder, "/ResultsATA.txt"), "w") do io
        write(io, "tests")
        write(io, "\r\n")
        DelimitedFiles.writedlm(io, collect(1:T)', ", ")
        write(io, "\r\n")
        write(io, "f⁺")
        write(io, "\r\n")
        DelimitedFiles.writedlm(io, ata_model.output.f)
        write(io, "\r\n")
        write(io, "infeasibility")
        write(io, "\r\n")
        DelimitedFiles.writedlm(io, ata_model.output.feas', ", ")
        write(io, "\r\n")
        write(io, "Item use")
        write(io, "\r\n")
        write(io, string(NH⁺.iu))
        write(io, "\r\n")
        write(io, "TIF")
        write(io, "\r\n")
        DelimitedFiles.writedlm(io, ata_model.output.TIF', ", ")
        write(io, "\r\n")
        write(io, "elapsed Time for optimization")
        write(io, "\r\n")
        write(io, string(ata_model.output.elapsed_time))
        return write(io, "\r\n")
    end
    return nothing
end
