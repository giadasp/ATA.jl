
"""
add_exp_score!(ATAmodel::AbstractModel)

# Description

Add the expected score constraints as specified in the `settings_file`.
If the names of the expected score columns of the item bank are not provided in the `settings_file` they are computed as the item characteristic functions following the IRT model always specified in the `settings_file`.
It requires the [`start_ATA`](#ATA.start_ATA) step.

# Arguments

- **`ATAmodel::AbstractModel`** : Required. The model built with `start_ATA()` and with settings loaded by [`start_ATA`](#ATA.start_ATA) function.

"""
function add_exp_score!(ATAmodel::AbstractModel)
    T = ATAmodel.settings.T
    n_items = ATAmodel.settings.n_items
    n_groups = ATAmodel.settings.n_groups
    expected_score_var =
        [ATAmodel.constraints[t].expected_score.var for t = 1:ATAmodel.settings.T]
    ICF = Vector{Matrix{Float64}}(undef, T)
    for t = 1:T
        if expected_score_var[t] == Symbol("")
            ICF[t] = item_char(
                ATAmodel.settings.IRT.parameters,
                ATAmodel.constraints[t].expected_score.pts,
                model = ATAmodel.settings.IRT.model,
                parametrization = ATAmodel.settings.IRT.parametrization,
                D = ATAmodel.settings.IRT.D,
            )[1][
                :,
                :,
                1,
            ]# K[t] x I
        else
            ICF[t] = zeros(Float64, 1, n_items)
            ICF[t][1, :] = ATAmodel.settings.bank[!, expected_score_var[t]]
        end
    end
    for t = 1:T
        ATAmodel.constraints[t].expected_score.val = ICF[t]
    end
    # for t = 1:ATAmodel.settings.T
    # 	DelimitedFiles.writedlm("OPT/A_$t.csv", ATAmodel.constraints[t].constr_A)
    # 	DelimitedFiles.writedlm("OPT/b_$t.csv", ATAmodel.constraints[t].constr_b)
    # end
    JLD2.@save "OPT/ICF.jld2" ICF
    push!(ATAmodel.output.infos, ["success", "- Expected Score constrained.\n"])
    return nothing
end
