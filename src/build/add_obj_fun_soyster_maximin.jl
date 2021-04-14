"""
    add_obj_fun!(
        ata_model::SoysterMaximinModel;
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
    ata_model::SoysterMaximinModel;
    psychometrics = false,
    items::Vector{Psychometrics.Item} = Psychometrics.Item[],
    items_file = "",
    kwargs...,
)

    try
        if psychometrics
            soyster_load_parameters_chain!(
                ata_model;
                items_file = items_file,
                items = items,
                kwargs...,
            )
            success!(ata_model, "IIFs for all item parameters samples computed.")

        end
    catch e
        error!(ata_model, string(sprint(showerror, e)))

    end
    return nothing
end
