"""
    add_obj_fun!(
        ata_model::CCMaximinModel;
        psychometrics = false,
        items_file = "",
        items::Vector{Psychometrics.Item} = Psychometrics.Item[],
        kwargs...
    )

# Description

It adds the objective function as specified in the `settings_file`. It requires the [`start_ata`](#ATA.start_ata) build step.  

Computes the IIFs at predefined ability points using the sampled item parameter estimates.
For each triplet of item, ability point, and item parameter sample \$(i, θ_k, r)\$, an IIF\$(θ_k)_{ir}\$ is computed.
   
The vector of items can be passed either through the argument `items` (`Vector{Psychometrics.Item}`) or through the argument `items_file` (string file name).
In both cases, the items in the vector must have the same order of the items in the item bank and they must contain the R sampled item parameters in the field `parameters.chain`.  

# Arguments

- **`ata_model::CCMaximinModel`** : Required. 
"""
function add_obj_fun!(
    ata_model::CCMaximinModel;
    psychometrics = false,
    items_file = "",
    items::Vector{Psychometrics.Item} = Psychometrics.Item[],
    kwargs...,
)
    message = ["", ""]
    try
        if !psychometrics
            R = ata_model.obj.R
            K = zeros(Int, T)
            IIF = Vector{Array{Float64,3}}(undef, T)
            ICF = Vector{Array{Float64,3}}(undef, T)
            if !isfile("BSPar.jld2")
                push!(
                    ata_model.output.infos,
                    [
                        "danger",
                        "- ccmaximin objective requires a jld2 file named \"BSPar.jld2\" with sampled values for the item parameters.\nAlternatively, a vector of `Psychometrics.Item` with samples in chain can be provided.",
                    ],
                )
                return nothing
            end
            JLD2.@load "BSPar.jld2" BSPar
            BSa = Matrix(BSPar[2]) |> x -> x[:, 2:end]
            BSb = Matrix(BSPar[1]) |> x -> x[:, 2:end]
            for t = 1:T
                K[t] = size(ata_model.obj.cores[t].points, 1)
                IIF[t] = zeros(K[t], n_items, R)
                ICF[t] = zeros(K[t], n_items, R)
                for r = 1:R
                    if irt_model == "1PL"
                        df = DataFrames.DataFrame(b = BSb[:, r]) #nqp values in interval\r\n",
                    elseif irt_model == "2PL"
                        df = DataFrames.DataFrame(a = BSa[:, r], b = BSb[:, r]) #nqp values in interval\r\n",
                    elseif irt_model == "3PL"
                        df = DataFrames.DataFrame(
                            a = BSa[:, r],
                            b = BSb[:, r],
                            c = BSc[:, r],
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
                    ata_model.obj.cores[t].IIF = IIF[t]
                end
            end
            JLD2.@save "opt/IIF_CC.jld2" IIF
            JLD2.@save "opt/ICF_CC.jld2" ICF
        else
            cc_maximin_load_parameters_chain!(
                ata_model;
                items_file = items_file,
                items = items,
                kwargs...,
            )
        end
        message = [
            "success",
            "- IIFs for all item parameters samples computed.\n",
        ]
        open("opt/Settings.jl", "a") do f
            write(f, "K = $K\n\n")
        end

        push!(ata_model.output.infos, message)
    catch e
        message[1] = "danger"
        message[2] = message[2] * string("- ", sprint(showerror, e), "\n")
        push!(ata_model.output.infos, message)
    end
    return nothing
end
