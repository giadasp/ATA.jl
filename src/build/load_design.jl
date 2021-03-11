"""
load_design!(design::Matrix{Any}, ATAmodel::AbstractModel)

# Description

Load a 0-1 `IxT` (or `n_FSxT` if the items are grouped by friend sets) design matrix into the ATA model.
Useful for loading a custom starting design before the `assemble!` step or to print/plot the features of the tests produced by a custom design before running `print_results`

# Arguments

- **`ATAmodel::AbstractModel`** : Required. The model built with `start_ATA()`, settings loaded by [`start_ATA`](#ATA.start_ATA).
"""
function load_design!(design::Matrix{Any}, ATAmodel::AbstractModel)
    ATAmodel.output.design = Float64.(design)
    return nothing
end
