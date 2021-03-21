"""
load_design!(design::Matrix{Any}, ata_model::AbstractModel)

# Description

Load a 0-1 `IxT` (or `n_fsxT` if the items are grouped by friend sets) design matrix into the ATA model.
Useful for loading a custom starting design before the `assemble!` step or to print/plot the features of the tests produced by a custom design before running `print_results`

# Arguments

- **`ata_model::AbstractModel`** : Required. The model built with `start_ATA()`, settings loaded by [`start_ATA`](#ATA.start_ATA).
"""
function load_design!(design::Matrix{Any}, ata_model::AbstractModel)
    ata_model.output.design = Float64.(design)
    push!(ata_model.output.infos, ["success", string("- Starting design loaded.\n")])
    return nothing
end
