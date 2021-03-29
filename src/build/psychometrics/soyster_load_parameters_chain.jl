
function soyster_load_parameters_chain!(
    ata_model::SoysterMaximinModel;
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
        IIF = Vector{Array{Float64,2}}(undef, T)
        ICF = Vector{Array{Float64,2}}(undef, T)
        K = zeros(Int, T)
        for t = 1:T
            K[t] = size(ata_model.obj.cores[t].points, 1)
            IIF[t] = zeros(K[t], n_items)
            ICF[t] = zeros(K[t], n_items)
            IIF_t_min = Inf .* ones(K[t], n_items)
            ICF_t_min = Inf .* ones(K[t], n_items)
            #check if chain has length R
            if size(items[1].parameters.chain, 1) != R
                error("Item 1 does not have R=", R, " chains.")
            end
            for i = 1:n_items
                IIF_i = zeros(R, K[t])
                ICF_i = zeros(R, K[t])
                for r = 1:R
                    if irt_model == "1PL"
                        if !(items[1].parameters isa Psychometrics.Parameters1PL)
                            error("Item 1 is not of type ", irt_model, ".")
                        end
                        df = DataFrames.DataFrame(b = [items[i].parameters.chain[r][1]])
                    elseif irt_model == "2PL"
                        if !(items[1].parameters isa Psychometrics.Parameters2PL)
                            error("Item 1 is not of type ", irt_model, ".")
                        end
                        df = DataFrames.DataFrame(
                            a = [items[i].parameters.chain[r][1]],
                            b = [items[i].parameters.chain[r][2]],
                        )
                    elseif irt_model == "3PL"
                        if !(items[1].parameters isa Psychometrics.Parameters3PL)
                            error("Item 1 is not of type ", irt_model, ".")
                        end
                        df = DataFrames.DataFrame(
                            a = [items[i].parameters.chain[r][1]],
                            b = [items[i].parameters.chain[r][2]],
                            c = [items[i].parameters.chain[r][3]],
                        )
                    end
                    for k = 1:K[t]
                        IIF_i[r, k] = item_info(
                            df,
                            ata_model.obj.cores[t].points[k];
                            model = irt_model,
                            parametrization = irt_parametrization,
                            D = irt_D,
                        )[1]
                        ICF_i[r, k] = item_char(
                            df,
                            ata_model.obj.cores[t].points[k];
                            model = irt_model,
                            parametrization = irt_parametrization,
                            D = irt_D,
                        )[1][
                            :,
                            :,
                            1,
                        ][1]
                    end
                end
                for k = 1:K[t]
                    IIF[t][k, i] = minimum(IIF_i[:, k])
                    ICF[t][k, i] = minimum(ICF_i[:, k])
                end
                ata_model.obj.cores[t].IIF = IIF[t]
            end
        end
    catch e
        message[1] = "danger"
        message[2] = message[2] * string("- ", sprint(showerror, e), "\n")
        push!(ata_model.output.infos, message)
    end
    return nothing
end
