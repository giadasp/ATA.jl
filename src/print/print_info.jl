"""
    print_infos(ata_model)

# Description

Print the results of each build step of the ATA model.

# Arguments

- **`ata_model::AbstractModel`** : Required. Model processed with build functions.
"""
function print_infos(ata_model::AbstractModel)
    for m in ata_model.output.infos
        if m[1] == "danger"
            printstyled(m[2]; color = :red)
        else
            printstyled(m[2]; color = :green)
        end
    end
end

"""
    print_last_info(ata_model)

# Description

Print info of the last build step.

# Arguments

- **`ata_model::AbstractModel`** : Required. Model processed with build functions.
"""
function print_last_info(ata_model::AbstractModel)
    m = ata_model.output.infos[end]
    if m[1] == "danger"
        printstyled(m[2]; color = :red)
    else
        printstyled(m[2]; color = :green)
    end
end
