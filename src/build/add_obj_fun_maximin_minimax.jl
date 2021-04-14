
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

    try
        compute_estimated_iif!(ata_model)
        success!(ata_model, "Estimated IIFs computed.")

    catch e
        error!(ata_model, string(sprint(showerror, e)))

    end
    return nothing
end
