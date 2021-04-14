"""
    load_design!(design::Matrix{Any}, ata_model::AbstractModel)

# Description

Load a 0-1 `IxT` (or `n_fsxT` if the items are grouped by friend sets) design matrix into the ATA model.
Useful for loading a custom starting design before the `assemble!` step or to print/plot the features of the tests produced by a custom design before running `print_results`

# Arguments

- **`ata_model::AbstractModel`** : Required. The model built with `start_ata()`, settings loaded by [`start_ata`](#ATA.start_ata).
"""
function load_design!(design::Matrix{Any}, ata_model::AbstractModel)

    try
        ata_model.output.design = Float64.(design)
        success!(ata_model, string("Starting design loaded."))
    catch e
        error!(ata_model, string(sprint(showerror, e)))

    end
    return nothing
end
