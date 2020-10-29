"""
	print_results(ATAmodel; group_by_fs = false, results_folder = "RESULTS")

# Description

Print the features of the assembled tests.

# Arguments

- **`ATAmodel::Model`** : Required. The model built with `ATA` fuctions, `ATAmodel.design` matrix must be `IxT` or `nFSxT` if the items are grouped by friend sets. 
- **`group_by_fs`** : Optional. Default: `false`. Set to `true` if items have been grouped by friend sets by [`group_by_friends!`](#ATA.group_by_friends!-Tuple{ATA.Model}).
- **`results_folder`** : Optional. Default: "RESULTS". The folder in which the output is stored.
"""
function print_results(ATAmodel; group_by_fs = false, results_folder = "RESULTS")
    if !(results_folder in readdir())
        mkdir(results_folder)
    else
        println(string(
            "You have already a folder with this name, files in ",
            results_folder,
            " will be overwritten.\n",
        ))
    end
    if size(ATAmodel.output.design, 1) > 0
        T = ATAmodel.settings.T
        n_items = ATAmodel.settings.n_items
        length_max = [ATAmodel.constraints[t].length_max for t = 1:T]
        categories = ATAmodel.output.categories
        n_FS = size(ATAmodel.settings.FS.items, 1)
        if n_FS < ATAmodel.settings.n_items
            n = n_FS * T
        else
            n = n_items * T
        end

        if ATAmodel.obj.type == "CC"
            JLD2.@load "OPT/IIF_CC.jld2" IIF
            JLD2.@load "OPT/ICF_CC.jld2" ICF
            IIF_CC = copy(IIF)
            ICF_CC = copy(ICF)
        end
        if ATAmodel.obj.type == "MAXIMIN" || ATAmodel.obj.type == "CC"
            JLD2.@load "OPT/IIF.jld2" IIF
        end
	if isfile("OPT/ICF.jld2")
	    JLD2.@load "OPT/ICF.jld2" ICF
	end
        if group_by_fs == true
            design = reshape(ATAmodel.output.design, n_FS, T)
            new_design = zeros(Float64, n_items, T)
            for t = 1:T
                new_design[
                    vcat(ATAmodel.settings.FS.items[findall(
                        design[:, t] .== one(Float64),
                    )]...),
                    t,
                ] .= one(Float64)
            end
            design = new_design
            DelimitedFiles.writedlm(string(results_folder, "/designItemLevel.csv"), design)
        else
            design = reshape(ATAmodel.output.design, n_items, T)
        end
        #overlap
        olMatrixOUT = zeros(Int64, T, T)
        for t = 1:T
            for v = 1:T
                olMatrixOUT[t, v] = design[:, t]' * design[:, v]
            end
        end
        #DelimitedFiles.writedlm(string(results_folder,"/olMatrixOUT.csv"),olMatrixOUT)
        #TIF e ICF
        if ATAmodel.obj.type == "MAXIMIN" || ATAmodel.obj.type == "CC"
            if isfile("simPool.csv")
                simPool = CSV.read("simPool.csv")
            else
                simPool = Float64[]
            end
        end
        #save values
        if ATAmodel.obj.type == "CC"
            min = Float64[]
            First = Float64[]
            Median = Float64[]
            Third = Float64[]
            max = Float64[]
            Alpha = Float64[]
            Estimated = Float64[]
            if isfile("simPool.csv")
                True = Float64[]
            end
            Distributed.@sync for t = 1:T
                K = size(ATAmodel.obj.points[t], 1)
                for k = 1:K
                    IIF_t = (IIF_CC[t][k, :, :]' * design[:, t])
                    min = vcat(min, minimum(IIF_t))
                    First = vcat(First, StatsBase.quantile(IIF_t, [0.25]))
                    Median = vcat(Median, StatsBase.quantile(IIF_t, [0.5]))
                    Third = vcat(Third, StatsBase.quantile(IIF_t, [0.75]))
                    max = vcat(max, maximum(IIF_t))
                    Alpha = vcat(Alpha, StatsBase.quantile(IIF_t, [ATAmodel.obj.aux_float]))
                    Estimated = vcat(Estimated, IIF[t][k, :]' * design[:, t])
                    if isfile("simPool.csv")
                        True = vcat(
                            True,
                            item_info(
                                simPool,
                                ATAmodel.obj.points[t][k],
                                model = ATAmodel.settings.IRT.model,
                            )' * design[:, t],
                        )
                    end
                end
            end
            TIFatTheta_k =
                DataFrames.DataFrame(hcat(min, First, Median, Third, max, Alpha, Estimated))

            if isfile("simPool.csv")
                TIFatTheta_k[:True] = True
                DataFrames.DataFrames.names!(
                    TIFatTheta_k,
                    Symbol.([
                        "min",
                        "0.25-Qle",
                        "Median",
                        "0.75-Qle",
                        "max",
                        string(ATAmodel.obj.aux_float, "-Qle"),
                        "Estimated",
                        "True",
                    ]),
                )
            else
                DataFrames.DataFrames.names!(
                    TIFatTheta_k,
                    Symbol.([
                        "min",
                        "0.25-Qle",
                        "Median",
                        "0.75-Qle",
                        "max",
                        string(ATAmodel.obj.aux_float, "-Qle"),
                        "Estimated",
                    ]),
                )
            end
            CSV.write(string(results_folder, "/TIFatTheta_k.csv"), TIFatTheta_k)

	elseif ATAmodel.obj.type == "MAXIMIN"
            Estimated = Float64[]
            if isfile("simPool.csv")
                True = Float64[]
            end
            for t = 1:T
                K = size(ATAmodel.obj.points[t], 1)
                for k = 1:K
                    Estimated = vcat(Estimated, IIF[t][k, :]' * design[:, t])
                    if isfile("simPool.csv")
                        True = vcat(
                            True,
                            item_info(
                                simPool,
                                ATAmodel.obj.points[t][k],
                                model = ATAmodel.settings.IRT.model,
                            )' * design[:, t],
                        )
                    end
                end
            end
            TIFatTheta_k = DataFrames.DataFrame(Estimated = Estimated)
            if isfile("simPool.csv")
                TIFatTheta_k[:True] = True
                DataFrames.DataFrames.names!(TIFatTheta_k, Symbol.(["Estimated", "True"]))
            end
            CSV.write(string(results_folder, "/TIFatTheta_k.csv"), TIFatTheta_k)

        end

        #expected score
        if isfile("OPT/ICF.jld2")
		ESprintIRT = Vector{Vector{Float64}}(undef, T)
		for t = 1:T
		    ESprintIRT[t] = zeros(size(ICF[t], 1))
		    for k = 1:size(ICF[t], 1)
			ESprintIRT[t][k] = ((ICF[t][k, :]'*design[:, t]))[1] / sum(design[:, t])
		    end
		end
        end
									
        #number of items
        n = sum(design, dims = 1)
        design_items = zeros(Int64, Int(maximum(n)), T)
        #design_items=Array{String,2}(nothing,n_max,T)
        for t = 1:T
            design_items[collect(1:Int(n[t])), t] .= sort!(findall(design[:, t] .== 1))#sort!(CODEVar[findall(design[:,t].==1)])
        end
        open(string(results_folder, "/Results.txt"), "w") do io
            write(io, "Length")
            write(io, "\r\n")
            DelimitedFiles.writedlm(io, [1:T, Int.(n)'], "\t")
            write(io, "\r\n")
            write(io, "design")
            write(io, "\r\n")
            DelimitedFiles.writedlm(io, design_items, "\t")
            write(io, "\r\n")
            write(io, "expected score")
            write(io, "\r\n")
	    if isfile("OPT/ICF.jld2")
		    for t = 1:T
			DelimitedFiles.writedlm(
			    io,
			    vcat(t, round.(ESprintIRT[t], digits = 3))',
			    "\t",
			)
			write(io, "\r\n")
		    end
            end
            if size(categories, 1) > 0
                write(io, "Categorical variables")
                write(io, "\r\n")
                for cats in categories
                    write(io, cats)
                    write(io, "\r\n")
                    vals1 =
                        setdiff(unique(ATAmodel.settings.bank[!, Symbol(cats)]), [missing])
                    vals = copy(vals1)
                    for t = 1:T
                        valscat = Int64[]
                        for val in vals1
                            valscat_t =
                                ATAmodel.settings.bank[!, Symbol(cats)][design[
                                    :,
                                    t,
                                ].==1] .== val
                            valscat_t = valscat_t[.!ismissing.(valscat_t)]
                            if size(valscat_t, 1) > 0
                                valscat = vcat(valscat, sum(valscat_t))
                            else
                                valscat = vcat(valscat, 0)
                            end
                        end
                        vals = hcat(vals, valscat)
                    end
                    DelimitedFiles.writedlm(io, vals, "\t")
                    #DelimitedFiles.writedlm(io,hcat(vals,[sum(skipmissing(ATAmodel.settings.bank[Symbol(cats)][design[:,t].==1].== i))   for i in vals,t=1:T]),", ")
                    write(io, "\r\n")
                end
                write(io, "Sum variables")
                write(io, "\r\n")
                for var in unique(vcat([ATAmodel.constraints[t].sum_vars for t = 1:T]...))
                    write(io, var)
                    write(io, "\r\n")
                    DelimitedFiles.writedlm(
                        io,
                        round.(
                            [
                                sum(skipmissing(
                                    ATAmodel.settings.bank[var] .* design[:, t],
                                )) for t = 1:T
                            ],
                            digits = 3,
                        )',
                        "\t",
                    )
                    write(io, "\r\n")
                end
                write(io, "\r\n")
                write(io, "Overlap matrix")
                write(io, "\r\n")
                olMatrixOUT =
                    hcat(vcat("", string.(collect(1:T))), vcat(collect(1:T)', olMatrixOUT))
                DelimitedFiles.writedlm(io, olMatrixOUT, "\t")
                write(io, "\r\n")
            end
            # if size(summarize,1)>0
            # write(io,"quantitative variables")
            # write(io,"\r\n")
            #   for summ in 1:size(summarize[:,1],1)
            # vals=[summarize[summ][2](skipmissing(ATAmodel.settings.bank[Symbol(summarize[summ][1])][design[:,t].==1]))   for t=1:T]
            # write(io,string(summarize[summ][1]))
            # write(io,",")
            # write(io,string(summarize[summ][2]))
            # write(io,",")
            # DelimitedFiles.writedlm(io,vals',",")
            # write(io,"\n")
            # end
            # end
            write(io, "Item use")
            write(io, "\r\n")
            iu = [sum(design[i, :]) for i = 1:ATAmodel.settings.n_items]
            ids = sortperm(iu; rev = true)
            sort!(iu; rev = true)
            for i = 1:ATAmodel.settings.n_items
                write(io, string(ids[i]))
                write(io, "\t")
                write(io, string(iu[i]))
                write(io, "\r\n")
            end
        end
    else
        println("No solution found")
    end
end
