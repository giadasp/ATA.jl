"""
    add_obj_fun!(
        ata_model::DeJongMaximinModel;
        psychometrics = false,
        items::Vector{Psychometrics.Item} = Psychometrics.Item[],
        items_file = "",
        kwargs...
    )

# Description

It adds the objective function as specified in the `settings_file`. It requires the [`start_ata`](#ATA.start_ata) build step.  

Computes the IIFs at predefined ability points using the sampled item parameter estimates.
For each pair of item and ability point \$(i, θ_k)\$, the IIF\$(θ_k)_i\$ is computed as the mean minus the standard deviation
of the IIF\$(θ_k)_{ir}\$ computed on the R sampled item parameter estimates, where \$r = 1, \\ldots, R\$.

The vector of items can be passed either through the argument `items` (`Vector{Psychometrics.Item}`) or through the argument `items_file` (string file name).
In both cases, the items in the vector must have the same order of the items in the item bank and they must contain the R sampled item parameters in the field `parameters.chain`.
    
- **`ata_model::DeJongMaximinModel`** : Required.
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
        if psychometrics
            de_jong_load_parameters_chain!(
                ata_model;
                items_file = items_file,
                items = items,
                kwargs...,
            )
        end
        message = ["success", "- IIFs for all item parameters samples computed.\n"]
        # open("opt/Settings.jl", "a") do f
        #     write(f, "K = $K\n\n")
        # end

        push!(ata_model.output.infos, message)
    catch e
        message[1] = "danger"
        message[2] = message[2] * string("- ", sprint(showerror, e), "\n")
        push!(ata_model.output.infos, message)
    end
    return nothing
end
