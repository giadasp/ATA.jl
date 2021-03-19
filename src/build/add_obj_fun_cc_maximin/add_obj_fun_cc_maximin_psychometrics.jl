@require Psychometrics = "ce8202d9-98c3-4990-890a-8616ce2c06f9" include("add_obj_fun_cc_maximin_psychometrics.jl")

function load_item_parameters_chain(
        ata_model::AbstractModel;
        items_file = "items.jld2",
        kwargs...
    )
    IRT_parameters = ata_model.settings.IRT.parameters
    IRT_model = ata_model.settings.IRT.model
    IRT_D = ata_model.settings.IRT.D
    IRT_parametrization = ata_model.settings.IRT.parametrization
    R = ata_model.obj.cores[1].R
    T = ata_model.settings.T
    n_items = ata_model.settings.n_items
    if !isfile(items_file)
        error(string("- ", items_file, " does not exist.\n"))
    end
    JLD2.@load items_file items
    K = zeros(Int, T)
    for t = 1:T
        K[t] = size(ata_model.obj.cores[t].points, 1)
        IIF[t] = zeros(K[t], n_items, R)
        ICF[t] = zeros(K[t], n_items, R)
        #check if chain has length R
        if size(items[1].parameters.chain,1) != R
            error("Item 1 does not have R=", R," chains.")
        end
        for r in 1:R
            if IRT_model == "1PL"
                if !(items[1].parameters isa Psychometrics.Parameters1PL)
                    error("Item 1 is not of type ", IRT_model,".")
                end
                df = DataFrames.DataFrame(
                    b = map(i -> i.parameters.chain[r], items)
                    ) 
            elseif IRT_model == "2PL"
                if !(items[1].parameters isa Psychometrics.Parameters2PL)
                    error("Item 1 is not of type ", IRT_model,".")
                end
                df = DataFrames.DataFrame(
                    a = map(i -> i.parameters.chain[r][1], items),
                    b = map(i -> i.parameters.chain[r][2], items)
                    ) 
            elseif IRT_model == "3PL"
                if !(items[1].parameters isa Psychometrics.Parameters3PL)
                    error("Item 1 is not of type ", IRT_model,".")
                end
                df = DataFrames.DataFrame(
                    a = map(i -> i.parameters.chain[r][1], items),
                    b = map(i -> i.parameters.chain[r][2], items),
                    c = map(i -> i.parameters.chain[r][3], items)
                    )
            end
            for k = 1:K[t]
                IIF[t][k, :, r] = item_info(
                    df,
                    ata_model.obj.cores[t].points[k];
                    model = IRT_model,
                    parametrization = IRT_parametrization,
                    D = IRT_D,
                ) # K[t] x I x R
                ICF[t][k, :, r] = item_char(
                    df,
                    ata_model.obj.cores[t].points[k];
                    model = IRT_model,
                    parametrization = IRT_parametrization,
                    D = IRT_D,
                )[1] # K[t] x I x R
            end
            ata_model.obj.cores[t].IIF = IIF[t]
        end
    JLD2.@save "OPT/IIF_CC.jld2" IIF
    JLD2.@save "OPT/ICF_CC.jld2" ICF
    return nothing
end