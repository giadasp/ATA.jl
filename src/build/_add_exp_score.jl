"""
    _add_exp_score!(ata_model::AbstractModel)

# Description

Add the expected score constraints as specified in the `settings_file`.
If the names of the expected score columns of the item bank are not provided in the `settings_file` they are computed as the item characteristic functions following the IRT model always specified in the `settings_file`.
It requires the [`start_ata`](#ATA.start_ata) step.

# Arguments

- **`ata_model::AbstractModel`** : Required. The model built with `start_ata()` and with settings loaded by [`start_ata`](#ATA.start_ata) function.

"""
function _add_exp_score!(ata_model::AbstractModel)
    message = ["", ""]
    T = ata_model.settings.T
    n_items = ata_model.settings.n_items
    n_groups = ata_model.settings.n_groups
    expected_score_var =
        [ata_model.constraints[t].expected_score.var for t = 1:ata_model.settings.T]
    ICF = Vector{Matrix{Float64}}(undef, T)
    for t = 1:T
        if expected_score_var[t] == Symbol("")
            ICF[t] = item_char(
                ata_model.settings.irt.parameters,
                ata_model.constraints[t].expected_score.pts,
                model = ata_model.settings.irt.model,
                parametrization = ata_model.settings.irt.parametrization,
                D = ata_model.settings.irt.D,
            )[1][
                :,
                :,
                1,
            ]# K[t] x I
        else
            ICF[t] = zeros(Float64, 1, n_items)
            ICF[t][1, :] = ata_model.settings.bank[!, expected_score_var[t]]
        end
    end
    for t = 1:T
        ata_model.constraints[t].expected_score.val = ICF[t]
    end
    # for t = 1:ata_model.settings.T
    # 	DelimitedFiles.writedlm("opt/A_$t.csv", ata_model.constraints[t].constr_A)
    # 	DelimitedFiles.writedlm("opt/b_$t.csv", ata_model.constraints[t].constr_b)
    # end
    JLD2.@save "opt/ICF.jld2" ICF
    push!(ata_model.output.infos, ["success", "- Expected Score constrained.\n"])
    return nothing
end
