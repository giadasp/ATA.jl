function soyster_load_parameters_chain!(
    ata_model::SoysterMaximinModel;
    items_file = "items.jld2",
    items::Vector{Psychometrics.Item} = Psychometrics.Item[],
    kwargs...,
)


    try
        T = ata_model.settings.T
        n_items = ata_model.settings.n_items
        irt_parameters = ata_model.settings.irt.parameters
        irt_model = ata_model.settings.irt.model
        irt_D = ata_model.settings.irt.D
        irt_parametrization = ata_model.settings.irt.parametrization
        R = ata_model.obj.R
        n_items = ata_model.settings.n_items
        if !isfile(items_file) || (size(items, 1) > 0)
            error!(
                ata_model,
                string(
                    "",
                    items_file,
                    " does not exist.\nProvide a `Psychometrics.Item` vector through the `items` argument, or the name of a file which can be loaded by `FileIO` that contains that vector through the `items_file` argument.",
                ),
            )
        elseif isfile(items_file)
            items = FileIO.load(items_file)
            items = items[collect(keys(items))[1]]
        end
        IIF = Vector{Array{Float64,2}}(undef, T)
        # ICF = Vector{Array{Float64,2}}(undef, T)
        K = zeros(Int, T)
        for t = 1:T
            K[t] = size(ata_model.obj.cores[t].points, 1)
            IIF[t] = zeros(K[t], n_items)
            # ICF[t] = zeros(K[t], n_items)
        end
        for i = 1:n_items
            item = items[i]
            if size(item.parameters.chain, 1) != R
                error!(ata_model, string("Item 1 does not have R=", R, " chains."))
            end
            if irt_model == "1PL"
                if !(item.parameters isa Psychometrics.Parameters1PL)
                    error!(ata_model, string("Item 1 is not of type ", irt_model, "."))
                end
                df = DataFrames.DataFrame(b = [item.parameters.chain[r][1]])
            elseif irt_model == "2PL"
                if !(item.parameters isa Psychometrics.Parameters2PL)
                    error!(ata_model, string("Item 1 is not of type ", irt_model, "."))
                end
                df = DataFrames.DataFrame(
                    a = [item.parameters.chain[r][1] for r = 1:R],
                    b = [item.parameters.chain[r][2] for r = 1:R],
                )
            elseif irt_model == "3PL"
                if !(item.parameters isa Psychometrics.Parameters3PL)
                    error!(ata_model, string("Item 1 is not of type ", irt_model, "."))
                end
                df = DataFrames.DataFrame(
                    a = [item.parameters.chain[r][1] for r = 1:R],
                    b = [item.parameters.chain[r][2] for r = 1:R],
                    c = [item.parameters.chain[r][3] for r = 1:R],
                )
            end
            for t = 1:T
                IIF_i = zeros(R, K[t])
                IIF_i_t_min = Inf .* ones(K[t])
                #check if chain has length R
                for k = 1:K[t]
                    IIF_i[:, k] = item_info(
                        df,
                        ata_model.obj.cores[t].points[k];
                        model = irt_model,
                        parametrization = irt_parametrization,
                        D = irt_D,
                    )
                    # ICF_i[:, k] = item_char(
                    #     df,
                    #     ata_model.obj.cores[t].points[k];
                    #     model = irt_model,
                    #     parametrization = irt_parametrization,
                    #     D = irt_D,
                    # )[1][
                    #     :,
                    #     :,
                    #     1,
                    # ]
                    IIF[t][k, i] = max(
                        0.0,
                        StatsBase.mean(IIF_i[:, k]) - 3 * StatsBase.std(IIF_i[:, k]),
                    )
                    # ICF[t][k, i] = StatsBase.mean(IIF_i[:, k]) - StatsBase.std(IIF_i[:, k])
                end
            end
        end
        for t = 1:T
            ata_model.obj.cores[t].IIF = IIF[t]
        end
    catch e
        error!(ata_model, string(sprint(showerror, e)))

    end
    return nothing
end
