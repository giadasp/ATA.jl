

"""
add_overlap!(ATAmodel::AbstractModel; overlap_file = "OverlapMatrix.csv", overlap_delim=";")

# Description

Add maximum overlap constraints to the `ATA.AbstractModel` as specified by the `overlap_file` which contains a `TxT` matrix defining the maximum number of shared items between each pair of tests.
It requires the [`start_ATA`](#ATA.start_ATA) step.

# Arguments

- **`ATAmodel::AbstractModel`** : Required. The model built with `start_ATA()` and with settings loaded by [`start_ATA`](#ATA.start_ATA) function.
- **`overlap_file`** : Optional. Default: "OverlapMatrix.csv". The path of the file containing the maximum overlap between each pair of tests in the form of a matrix with custom-separated values.
- **`overlap_delim`** : Optional. Default: ";". The custom-separator for the `overlap_file`.

"""
function add_overlap!(
    ATAmodel::AbstractModel;
    overlap_file = "OverlapMatrix.csv",
    overlap_delim = ";",
)
    T = ATAmodel.settings.T
    n_items = ATAmodel.settings.n_items
    if !isfile(overlap_file)
        push!(ATAmodel.output.infos, ["danger", string(overlap_file, " not found.\n")])
        return nothing
    else
        opMatrix = Matrix{Int64}(
            CSV.read(
                overlap_file,
                DataFrames.DataFrame,
                delim = overlap_delim,
                header = false,
            ),
        )
        if size(opMatrix, 1) > 0
            #ol_max = Vector{Vector{Int64}}(undef, T)
            # for t = 1:T
            # 	ol_max[t] = opMatrix[t, setdiff(collect(1:T), t)]
            # end
            opMatrix = opMatrix[1:T, 1:T]
            DelimitedFiles.writedlm("OPT/overlap.txt", opMatrix)
            ATAmodel.settings.ol_max = opMatrix
        else
            ATAmodel.settings.ol_max =
                ones(ATAmodel.settings.n_items, ATAmodel.settings.n_items) .*
                ATAmodel.settings.n_items
            push!(
                ATAmodel.output.infos,
                [
                    "success",
                    string(
                        "- No lines in ",
                        overlap_file,
                        ", Maximum overlap set at ",
                        ATAmodel.settings.n_items,
                        ".\n",
                    ),
                ],
            )
        end
        JLD2.@save "OPT/ATAmodel.jld2" ATAmodel
        push!(ATAmodel.output.infos, ["success", "- Maximum overlap constrained.\n"])
    end
    return nothing
end
