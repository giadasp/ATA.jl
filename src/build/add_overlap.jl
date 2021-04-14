

"""
    add_overlap!(ata_model::AbstractModel; overlap_file = "overlap_matrix.csv", overlap_delim=";")

# Description

Add maximum overlap constraints to the `ATA.AbstractModel` as specified by the `overlap_file` which contains a `TxT` matrix defining the maximum number of shared items between each pair of tests.
Alternatively, the overlap matrix can be passed using the argument `overlap`.
It requires the [`start_ata`](#ATA.start_ata) step.

# Arguments

- **`ata_model::AbstractModel`** : Required. The model built with `start_ata()` and with settings loaded by [`start_ata`](#ATA.start_ata) function.
- **`overlap::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)`**: Optional. 
- **`overlap_file`** : Optional. Default: "overlap_matrix.csv". The path of the file containing the maximum overlap between each pair of tests in the form of a matrix with custom-separated values.
- **`overlap_delim`** : Optional. Default: ";". The custom-separator for the `overlap_file`.

"""
function add_overlap!(
    ata_model::AbstractModel;
    overlap::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),
    overlap_file = "overlap_matrix.csv",
    overlap_delim = ";",
)

    #load overlap matrix
    if size(overlap, 1) > 0
        overlap_matrix = overlap
        success!(ata_model, "Overlap matrix loaded.")
    elseif isfile(overlap_file)
        try
            overlap_matrix = CSV.read(
                overlap_file,
                DataFrames.DataFrame,
                delim = overlap_delim,
                header = false,
            )
        catch e
            error!(ata_model, "Error in reading the overlap file.")
            return nothing
        end
        success!(ata_model, "Overlap file loaded.")
    else
        error!(
            ata_model,
            string(
                "Overlap file with name ",
                overlap_file,
                " does not exist. Provide a valid overlap matrix or a name of an existing file.",
            ),
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
            DelimitedFiles.writedlm("opt/overlap.txt", overlap_matrix)
            ata_model.settings.ol_max = overlap_matrix
            ata_model.settings.to_apply[3] = true
        else
            ata_model.settings.ol_max =
                ones(ata_model.settings.n_items, ata_model.settings.n_items) .*
                ata_model.settings.n_items
            ata_model.settings.to_apply[3] = false
            success!(
                ata_model,
                string(
                    "No lines in ",
                    overlap_file,
                    ", overlap constraints are not applied.",
                ),
            )
        end
        JLD2.@save "opt/ata_model.jld2" ata_model
        success!(ata_model, "Maximum overlap constrained.")
    catch e
        error!(ata_model, string(sprint(showerror, e)))

    end
    return nothing
end
