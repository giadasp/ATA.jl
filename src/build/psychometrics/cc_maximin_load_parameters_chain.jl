function cc_maximin_load_parameters_chain!(
    ata_model::CCMaximinModel;
    items_file = "items.jld2",
    items::Vector{Psychometrics.Item} = Psychometrics.Item[],
    kwargs...,
)
    message = ["", ""]

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
            error(
                string(
                    "- ",
                    items_file,
                    " does not exist.\nProvide a `Psychometrics.Item` vector through the `items` argument, or the name of a file which can be loaded by `FileIO` that contains that vector through the `items_file` argument.\n",
                ),
            )
        elseif isfile(items_file)
            items = FileIO.load(items_file)
            items = items[collect(keys(items))[1]]
        end
        IIF = Vector{Array{Float64,3}}(undef, T)
        ICF = Vector{Array{Float64,3}}(undef, T)
        K = zeros(Int, T)
        for t = 1:T
            K[t] = size(ata_model.obj.cores[t].points, 1)
            IIF[t] = zeros(K[t], n_items, R)
            ICF[t] = zeros(K[t], n_items, R)
            #check if chain has length R
            if size(items[1].parameters.chain, 1) != R
                error("Item 1 does not have R=", R, " chains.")
            end
            for r = 1:R
                if irt_model == "1PL"
                    if !(items[1].parameters isa Psychometrics.Parameters1PL)
                        error("Item 1 is not of type ", irt_model, ".")
                    end
                    df = DataFrames.DataFrame(b = map(i -> i.parameters.chain[r], items))
                elseif irt_model == "2PL"
                    if !(items[1].parameters isa Psychometrics.Parameters2PL)
                        error("Item 1 is not of type ", irt_model, ".")
                    end
                    df = DataFrames.DataFrame(
                        a = map(i -> i.parameters.chain[r][1], items),
                        b = map(i -> i.parameters.chain[r][2], items),
                    )
                elseif irt_model == "3PL"
                    if !(items[1].parameters isa Psychometrics.Parameters3PL)
                        error("Item 1 is not of type ", irt_model, ".")
                    end
                    df = DataFrames.DataFrame(
                        a = map(i -> i.parameters.chain[r][1], items),
                        b = map(i -> i.parameters.chain[r][2], items),
                        c = map(i -> i.parameters.chain[r][3], items),
                    )
                end
                for k = 1:K[t]
                    IIF[t][k, :, r] = item_info(
                        df,
                        ata_model.obj.cores[t].points[k];
                        model = irt_model,
                        parametrization = irt_parametrization,
                        D = irt_D,
                    ) # K[t] x I x R
                    ICF[t][k, :, r] = item_char(
                        df,
                        ata_model.obj.cores[t].points[k];
                        model = irt_model,
                        parametrization = irt_parametrization,
                        D = irt_D,
                    )[1][
                        :,
                        :,
                        1,
                    ] # K[t] x I x R
                end
            end
            ata_model.obj.cores[t].IIF = IIF[t]
        end
        JLD2.@save "opt/IIF_CC.jld2" IIF
        JLD2.@save "opt/ICF_CC.jld2" ICF
    catch e
        message[1] = "danger"
        message[2] = message[2] * string("- ", sprint(showerror, e), "\n")
        push!(ata_model.output.infos, message)
    end
    return nothing
end
