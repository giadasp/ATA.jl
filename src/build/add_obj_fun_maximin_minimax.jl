
"""
add_obj_fun!(ata_model::Union{MaximinModel, MinimaxModel})

# Description

Add the objective function as specified in the `settings_file`. It requires the [`start_ata`](#ATA.start_ata) build step.  
Computes the IIFs at predefined ability points.

# Arguments

- **`ata_model::Union{MaximinModel, MinimaxModel}`** : Required. The model built with `start_ata()` and with settings loaded by [`start_ata`](#ATA.start_ata) function.

"""
function add_obj_fun!(ata_model::Union{MaximinModel,MinimaxModel})
    message = ["",""]
    try
        T = ata_model.settings.T
        n_items = ata_model.settings.n_items
        IRT_parameters = ata_model.settings.IRT.parameters
        IRT_model = ata_model.settings.IRT.model
        IRT_D = ata_model.settings.IRT.D
        IRT_parametrization = ata_model.settings.IRT.parametrization
        IIF = Vector{Array{Float64,2}}(undef, T)
        ICF = Vector{Array{Float64,2}}(undef, T)
        K = zeros(Int, T)
        for t = 1:T
            K[t] = size(ata_model.obj.cores[t].points, 1)
            IIF[t] = zeros(K[t], n_items)
            ICF[t] = zeros(K[t], n_items)
            for k = 1:K[t]
                IIF[t][k, :] = item_info(
                    IRT_parameters,
                    ata_model.obj.cores[t].points[k],
                    model = IRT_model,
                    parametrization = IRT_parametrization,
                    D = IRT_D,
                )# K[t] x I
                ICF[t][k, :] = item_char(
                    IRT_parameters,
                    ata_model.obj.cores[t].points[k],
                    model = IRT_model,
                    parametrization = IRT_parametrization,
                    D = IRT_D,
                )[1][
                    :,
                    :,
                    1,
                ] # K[t] x I
            end
            ata_model.obj.cores[t].IIF = IIF[t]
        end
        message = ["success", "- Objective function loaded.\n- IIFs and ICFs computed.\n"]
        JLD2.@save "OPT/IIF.jld2" IIF
        if !isfile("OPT/ICF.jld2")
            JLD2.@save "OPT/ICF.jld2" ICF
        end
        push!(ata_model.output.infos, message)
    catch e
        message[1] = "danger"
        message[2] = message[2] * string("- ",sprint(showerror, e),"\n")
        push!(ata_model.output.infos, message)
    end
    return nothing
end
