"""
    add_obj_fun!(
        ata_model::RobustMaximinModel;
        psychometrics = false,
        items::Vector{Psychometrics.Item} = Psychometrics.Item[],
        items_file = "",
        kwargs...
    )

# Description

It adds the objective function as specified in the `settings_file`. It requires the [`start_ata`](#ATA.start_ata) build step.  
Computes the IIFs at predefined ability points using the estimated item parameters in the item bank.
It computes also the standard deviations of the IIFs given R samples of the item parameters provided by a vector of `Psychometrics.Item`.
The vector of items can be passed either through the argument `items` (`Vector{Psychometrics.Item}`) or through the argument `items_file` (string file name).
In both cases, the items in the vector must have the same order of the items in the item bank and they must contain the R sampled item parameters in the field `parameters.chain`.

- **`ata_model::RobustMaximinModel`** : Required.
"""
function add_obj_fun!(
    ata_model::RobustMaximinModel;
    psychometrics = false,
    items::Vector{Psychometrics.Item} = Psychometrics.Item[],
    items_file = "",
    kwargs...,
)

    try
        if psychometrics
            compute_estimated_iif!(ata_model)
            robust_load_parameters_chain!(
                ata_model;
                items_file = items_file,
                items = items,
                kwargs...,
            )
        end
        success!(ata_model, "Objective function loaded.")
        success!(ata_model, "IIFs computed.")
        success!(ata_model, "Standard deviations computed.")

        # open("opt/Settings.jl", "a") do f
        #     write(f, "K = $K\n\n")
        # end


    catch e
        error!(ata_model, string(sprint(showerror, e)))

    end
    return nothing
end
