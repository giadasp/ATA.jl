

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
    overlap::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),
    overlap_file = "overlap_matrix.csv",
    overlap_delim = ";",
)
    message = ["", ""]
    #load overlap matrix
    if size(overlap, 1) > 0
        overlap_matrix = overlap
        message[2] = message[2] * "- Overlap matrix loaded.\n"
    elseif isfile(overlap_file)
        try
            overlap_matrix = CSV.read(
                overlap_file,
                DataFrames.DataFrame,
                delim = overlap_delim,
                header = false,
            )
        catch e
            message[1] = "danger"
            message[2] = message[2] * "- Error in reading the overlap file.\n"
            push!(ata_model.output.infos, message)
            return nothing
        end
        message[2] = message[2] * "- Overlap file loaded.\n"
    else
        push!(
            ata_model.output.infos,
            [
                "danger",
                "Overlap file with name ",
                overlap_file,
                " does not exist. Provide a valid overlap matrix or a name of an existing file.",
            ],
        )
        return nothing
    end
    try
        T = ata_model.settings.T
        n_items = ata_model.settings.n_items

        if size(overlap_matrix, 1) > 0
            #ol_max = Vector{Vector{Int64}}(undef, T)
            # for t = 1:T
            # 	ol_max[t] = overlap_matrix[t, setdiff(collect(1:T), t)]
            # end
            overlap_matrix = Matrix{Float64}(overlap_matrix)[1:T, 1:T]
            DelimitedFiles.writedlm("OPT/overlap.txt", overlap_matrix)
            ata_model.settings.ol_max = overlap_matrix
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
    catch e
        message[1] = "danger"
        message[2] = message[2] * string("- ",sprint(showerror, e),"\n")
        push!(ata_model.output.infos, message)
    end
    return nothing
end
