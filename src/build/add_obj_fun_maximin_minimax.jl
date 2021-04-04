
"""
    add_obj_fun!(ata_model::Union{MaximinModel,MinimaxModel}; kwargs...)

# Description

It adds the objective function as specified in the `settings_file`. It requires the [`start_ata`](#ATA.start_ata) build step.  

Computes the IIFs at predefined ability points using the item parameter estimates in the item bank.
For each pair of item, and ability point \$(i, θ_k)\$, an IIF\$(θ_k)_{i}\$ is computed.
           
# Arguments

- **`ata_model::Union{MaximinModel, MinimaxModel}`** : Required. 
"""
function add_obj_fun!(ata_model::Union{MaximinModel,MinimaxModel}; kwargs...)
    message = ["", ""]
    try
        compute_estimated_iif!(ata_model)
        message = ["success", "- Estimated IIFs computed.\n"]
        push!(ata_model.output.infos, message)
    catch e
        message[1] = "danger"
        message[2] = message[2] * string("- ", sprint(showerror, e), "\n")
        push!(ata_model.output.infos, message)
    end
    return nothing
end
