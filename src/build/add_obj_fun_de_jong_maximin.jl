"""
    add_obj_fun!(
        ata_model::DeJongMaximinModel;
        psychometrics = false,
        items::Vector{Psychometrics.Item} = Psychometrics.Item[],
        items_file = "",
        kwargs...
    )

# Description

Add the objective function as specified in the `settings_file`. It requires the [`start_ata`](#ATA.start_ata) build step.  
Computes the IIFs at predefined ability points using `R` sampled item parameters, and set the mean minus the standard deviation
across all the samples as the IIFs.

# Arguments

- **`ata_model::DeJongMaximinModel`** : Required. The model built with `start_ata()` and with settings loaded by [`start_ata`](#ATA.start_ata) function.

"""
function add_obj_fun!(
    ata_model::DeJongMaximinModel;
    psychometrics = false,
    items::Vector{Psychometrics.Item} = Psychometrics.Item[],
    items_file = "",
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
        R = ata_model.obj.R
        K = zeros(Int, T)
        IIF = Vector{Array{Float64,2}}(undef, T)
        ICF = Vector{Array{Float64,2}}(undef, T)
        if psychometrics
            de_jong_load_parameters_chain!(
                ata_model;
                items_file = items_file,
                items = items,
                kwargs...,
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
