"""
    add_obj_fun!(
        ata_model::CCMaximinModel;
        psychometrics = false,
        items_file = "",
        items::Vector{Psychometrics.Item} = Psychometrics.Item[],
        kwargs...
    )

# Description

Add the objective function as specified in the `settings_file`. It requires the [`start_ata`](#ATA.start_ata) build step.  
Computes the IIFs at predefined ability points using `R` sampled item parameters.

# Arguments

- **`ata_model::CCMaximinModel`** : Required. The model built with `start_ata()` and with settings loaded by [`start_ata`](#ATA.start_ata) function.

"""
function add_obj_fun!(
    ata_model::CCMaximinModel;
    psychometrics = false,
    items_file = "",
    items::Vector{Psychometrics.Item} = Psychometrics.Item[],
    kwargs...
    )
    message = ["", ""]
    try
        T = ata_model.settings.T
        n_items = ata_model.settings.n_items
        irt_parameters = ata_model.settings.irt.parameters
        irt_model = ata_model.settings.irt.model
        irt_D = ata_model.settings.irt.D
        irt_parametrization = ata_model.settings.irt.parametrization
        IIF = Vector{Array{Float64,2}}(undef, T)
        ICF = Vector{Array{Float64,2}}(undef, T)
        K = zeros(Int, T)
        for t = 1:T
            K[t] = size(ata_model.obj.cores[t].points, 1)
            IIF[t] = zeros(K[t], n_items)
            ICF[t] = zeros(K[t], n_items)
            for k = 1:K[t]
                IIF[t][k, :] = item_info(
                    irt_parameters,
                    ata_model.obj.cores[t].points[k],
                    model = irt_model,
                    parametrization = irt_parametrization,
                    D = irt_D,
                )# K[t] x I
                ICF[t][k, :] = item_char(
                    irt_parameters,
                    ata_model.obj.cores[t].points[k],
                    model = irt_model,
                    parametrization = irt_parametrization,
                    D = irt_D,
                )[1][
                    :,
                    :,
                    1,
                ] # K[t] x I
            end
        end
        JLD2.@save "OPT/IIF.jld2" IIF
        JLD2.@save "OPT/ICF.jld2" ICF
        R = ata_model.obj.cores[1].R
        K = zeros(Int, T)
        IIF = Vector{Array{Float64,3}}(undef, T)
        ICF = Vector{Array{Float64,3}}(undef, T)
        if !psychometrics
            if !isfile("BSPar.jld2")
                push!(
                    ata_model.output.infos,
                    [
                        "danger",
                        "- CCMAXIMIN objective requires a jld2 file \"BSPar.jld2\" with sampled values for the item parameters.",
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
                        )[1][:,:,1] # K[t] x I x R
                    end
                    ata_model.obj.cores[t].IIF = IIF[t]
                end
            end
            JLD2.@save "OPT/IIF_CC.jld2" IIF
            JLD2.@save "OPT/ICF_CC.jld2" ICF
        else
            load_item_parameters_chain!(
                ata_model;
                items_file = items_file,
                items = items,
                kwargs...
            )
        end
        message = [
            "success",
            "- Objective function loaded.\n- IIFs and ICFs computed.\n- IIFs and ICFs for all item parameters samples computed.\n",
        ]
        open("OPT/Settings.jl", "a") do f
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
