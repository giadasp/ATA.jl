"""
    print_results(
        ata_model::CcMaximinModel;
        group_by_fs = false,
        results_folder = "results",
        sim_pool::DataFrames.DataFrame = DataFrame(),
    )

# Description

Print the features of the assembled tests.

# Arguments

- **`ata_model::CcMaximinModel`** : Required. The model built with `ATA` fuctions, `ata_model.design` matrix must be `IxT` or `nfsxT` if the items are grouped by friend sets. 
- **`group_by_fs`** : Optional. Default: `false`. Set to `true` if items have been grouped by friend sets by [`group_by_friends!`](#ATA.group_by_friends!-Tuple{ATA.Model}).
- **`results_folder`** : Optional. Default: "results". The folder in which the output is stored.
- **`sim_pool::DataFrames.DataFrame`** : Optional. Default: DataFrame(). The pool with true item paramaters. For simulation studies.

"""
function print_results(
    ata_model::CcMaximinModel;
    group_by_fs = false,
    results_folder = "results",
    sim_pool::DataFrames.DataFrame = DataFrame(),
)
    if !isdir(results_folder)
        mkdir(results_folder)
    else
        println(
            string(
                "There is already a folder with this name, files in ",
                results_folder,
                " will be overwritten.\n",
            ),
        )
    end
    if size(ata_model.output.design, 1) > 0
        T = ata_model.settings.T
        n_items = ata_model.settings.n_items
        length_max = [ata_model.constraints[t].length_max for t = 1:T]
        categories = ata_model.output.categories
        n_fs = size(ata_model.settings.fs.items, 1)
        if n_fs < ata_model.settings.n_items
            n = n_fs * T
        else
            n = n_items * T
        end
        if group_by_fs == true
            design = reshape(ata_model.output.design, n_fs, T)
            new_design = zeros(Float64, n_items, T)
            for t = 1:T
                new_design[
                    vcat(
                        ata_model.settings.fs.items[findall(
                            design[:, t] .== one(Float64),
                        )]...,
                    ),
                    t,
                ] .= one(Float64)
            end
            design = new_design
            DelimitedFiles.writedlm(string(results_folder, "/designItemLevel.csv"), design)
        else
            design = reshape(ata_model.output.design, n_items, T)
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

        #save values
        zero_qle = Float64[]
        first_qle = Float64[]
        second_qle = Float64[]
        third_qle = Float64[]
        fourth_qle = Float64[]
        alpha = Float64[]
        estimated = Float64[]
        if size(sim_pool, 1) > 0
            True = Float64[]
        end
        for t = 1:T
            K = size(ata_model.obj.cores[t].points, 1)
            IIF_t = ata_model.obj.cores[t].IIF
            for k = 1:K
                IIF_t_k = (IIF_t[k, :, :]' * design[:, t])
                zero_qle = vcat(zero_qle, minimum(IIF_t_k))
                first_qle = vcat(first_qle, StatsBase.quantile(IIF_t_k, [0.25]))
                second_qle = vcat(second_qle, StatsBase.quantile(IIF_t_k, [0.5]))
                third_qle = vcat(third_qle, StatsBase.quantile(IIF_t_k, [0.75]))
                fourth_qle = vcat(fourth_qle, maximum(IIF_t_k))
                alpha = vcat(alpha, StatsBase.quantile(IIF_t_k, ata_model.obj.cores[t].alpha))
                estimated = vcat(
                    estimated,
                    item_info(
                        ata_model.settings.irt.parameters,
                        ata_model.obj.cores[t].points[k],
                        model = ata_model.settings.irt.model,
                        parametrization = ata_model.settings.irt.parametrization,
                        D = ata_model.settings.irt.D,
                    )' * design[:, t],
                )
                if size(sim_pool, 1) > 0
                    True = vcat(
                        True,
                        item_info(
                            sim_pool,
                            ata_model.obj.cores[t].points[k],
                            model = ata_model.settings.irt.model,
                            parametrization = ata_model.settings.irt.parametrization,
                            D = ata_model.settings.irt.D,
                        )' * design[:, t],
                    )
                end
            end
        end
        tif_at_theta_pts =
            DataFrames.DataFrame(hcat(zero_qle, first_qle, second_qle, third_qle, fourth_qle, alpha, estimated))

        if size(sim_pool, 1) > 0
            tif_at_theta_pts[:true] = True
            DataFrames.DataFrames.names!(
                tif_at_theta_pts,
                Symbol.([
                    "min",
                    "0.25-qle",
                    "median",
                    "0.75-qle",
                    "max",
                    string(ata_model.obj.cores[1].alpha, "-qle"),
                    "estimated",
                    "true",
                ]),
            )
        else
            DataFrames.DataFrames.names!(
                tif_at_theta_pts,
                Symbol.([
                    "min",
                    "0.25-qle",
                    "median",
                    "0.75-qle",
                    "max",
                    string(ata_model.obj.cores[1].alpha, "-qle"),
                    "estimated",
                ]),
            )
        end
        CSV.write(string(results_folder, "/tif_at_theta_pts.csv"), tif_at_theta_pts)

        ICF = map(c -> c.expected_score.val, ata_model.constraints)
        for t = 1:T
            if size(ICF, 1) == 0
                ICF[t] =
                    item_char(
                        ata_model.settings.irt.parameters,
                        ata_model.obj.cores[t].points[k],
                        model = ata_model.settings.irt.model,
                        parametrization = ata_model.settings.irt.parametrization,
                        D = ata_model.settings.irt.D,
                    )' * design[:, t]
            end
        end

        #expected score
        es_print_irt = Vector{Vector{Float64}}(undef, T)
        for t = 1:T
            es_print_irt[t] = zeros(size(ICF[t], 1))
            for k = 1:size(ICF[t], 1)
                es_print_irt[t][k] = ((ICF[t][k, :]'*design[:, t]))[1] / sum(design[:, t])
            end
        end

        #number of items
        n = sum(design, dims = 1)
        design_items = zeros(Int64, Int(maximum(n)), T)
        #design_items=Array{String,2}(nothing,n_max,T)
        for t = 1:T
            design_items[collect(1:Int(n[t])), t] .= sort!(findall(design[:, t] .== 1))#sort!(CODEVar[findall(design[:,t].==1)])
        end
        open(string(results_folder, "/results.txt"), "w") do io
            write(io, "Length")
            write(io, "\r\n")
            DelimitedFiles.writedlm(io, [1:T, Int.(n)'], "\t")
            write(io, "\r\n")
            write(io, "design")
            write(io, "\r\n")
            DelimitedFiles.writedlm(io, design_items, "\t")
            write(io, "\r\n")
            write(io, "Expected score")
            write(io, "\r\n")
            for t = 1:T
                DelimitedFiles.writedlm(
                    io,
                    vcat(t, round.(es_print_irt[t], digits = 3))',
                    "\t",
                )
                write(io, "\r\n")
            end
            if size(categories, 1) > 0
                write(io, "Categorical variables")
                write(io, "\r\n")
                for cats in categories
                    write(io, cats)
                    write(io, "\r\n")
                    vals1 =
                        setdiff(unique(ata_model.settings.bank[!, Symbol(cats)]), [missing])
                    vals = copy(vals1)
                    for t = 1:T
                        valscat = Int64[]
                        for val in vals1
                            valscat_t =
                                ata_model.settings.bank[!, Symbol(cats)][design[
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
                    #DelimitedFiles.writedlm(io,hcat(vals,[sum(skipmissing(ata_model.settings.bank[Symbol(cats)][design[:,t].==1].== i))   for i in vals,t=1:T]),", ")
                    write(io, "\r\n")
                end
                write(io, "Sum variables")
                write(io, "\r\n")
                for var in unique(vcat([ata_model.constraints[t].sum_vars for t = 1:T]...))
                    write(io, var)
                    write(io, "\r\n")
                    DelimitedFiles.writedlm(
                        io,
                        round.(
                            [
                                sum(
                                    skipmissing(
                                        ata_model.settings.bank[!, var] .* design[:, t],
                                    ),
                                ) for t = 1:T
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
            write(io, "Item use")
            write(io, "\r\n")
            iu = [sum(design[i, :]) for i = 1:ata_model.settings.n_items]
            ids = sortperm(iu; rev = true)
            sort!(iu; rev = true)
            for i = 1:ata_model.settings.n_items
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

