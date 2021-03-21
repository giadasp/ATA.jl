

"""
add_overlap!(ata_model::AbstractModel; overlap_file = "overlap_matrix.csv", overlap_delim=";")

# Description

Add maximum overlap constraints to the `ATA.AbstractModel` as specified by the `overlap_file` which contains a `TxT` matrix defining the maximum number of shared items between each pair of tests.
It requires the [`start_ATA`](#ATA.start_ATA) step.

# Arguments

- **`ata_model::AbstractModel`** : Required. The model built with `start_ATA()` and with settings loaded by [`start_ATA`](#ATA.start_ATA) function.
- **`overlap_file`** : Optional. Default: "overlap_matrix.csv". The path of the file containing the maximum overlap between each pair of tests in the form of a matrix with custom-separated values.
- **`overlap_delim`** : Optional. Default: ";". The custom-separator for the `overlap_file`.

"""
function add_overlap!(
    ata_model::AbstractModel;
    overlap_file = "overlap_matrix.csv",
    overlap_delim = ";",
)
    T = ata_model.settings.T
    n_items = ata_model.settings.n_items
    if !isfile(overlap_file)
        push!(ata_model.output.infos, ["danger", string(overlap_file, " not found.\n")])
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
            ata_model.settings.ol_max = opMatrix
            ata_model.settings.to_apply[3] = true
        else
            ata_model.settings.ol_max =
                ones(ata_model.settings.n_items, ata_model.settings.n_items) .*
                ata_model.settings.n_items
            ata_model.settings.to_apply[3] = false
            push!(
                ata_model.output.infos,
                [
                    "success",
                    string(
                        "- No lines in ",
                        overlap_file,
                        ", overlap constraints are not applied.\n",
                    ),
                ],
            )
        end
        JLD2.@save "OPT/ata_model.jld2" ata_model
        push!(ata_model.output.infos, ["success", "- Maximum overlap constrained.\n"])
    end
    return nothing
end
