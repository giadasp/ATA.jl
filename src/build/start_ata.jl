"""
    start_ATA(;
        settings_file = "SettingsATA.jl",
        bank_file = "Bank.csv",
        bank_delim = ";",
    )

# Description

Initialize an empty ATA model and load the test assembly settings contained in the `settings_file`.
It requires the [`start_ATA`](#ATA.start_ATA) step.

# Arguments

- **`ATAmodel::AbstractModel`** : Required. The model built with `start_ATA()` function.
- **`settings_file`** : Optional. Default: "SettingsATA.jl". The path of the file containing the ATA settings in the form of an `InputSettings` struct.
- **`bank_file`** : Optional. Default: "Bank.csv". The path of the file containing the item pool/bank in the form of custom-separated values.
- **`bank_delim`** : Optional. Default: ";". The custom-separator for the bank_file.

# Output

- An ATA model.

"""
function start_ATA(;
    settings_file = "SettingsATA.jl",
    bank_file = "Bank.csv",
    bank_delim = ";",
)
    message = ["", ""]

    if !isfile(settings_file)
        error("Settings file not valid")
    else
        include(settings_file)

        if Inputs.obj_type != ""
            if Inputs.obj_type == "MAXIMIN"
                ATAmodel = MaximinModel()
            elseif Inputs.obj_type == "MINIMAX"
                ATAmodel = MinimaxModel()
            elseif Inputs.obj_type == "CCMAXIMIN"
                ATAmodel = CcMaximinModel()
            elseif Inputs.obj_type == "custom"
                ATAmodel = CustomModel()
            else
                error(
                    "Only \"MAXIMIN\", \"MINIMAX\", \"CCMAXIMIN\", \"custom\" and \"\" objective types are supported.",
                )
            end
            message[2] =
                message[2] * string("- ", ATAmodel.obj.name, " optimization type loaded.\n")
        else
            ATAmodel = NoObjModel()
        end
        ATAmodel.settings.n_items = Inputs.n_items
        ATAmodel.settings.T = Int(sum(Inputs.T))
        ATAmodel.settings.Tg = Inputs.T
        #initialize constraints
        ATAmodel.constraints = [Constraint() for t = 1:ATAmodel.settings.T]
        ATAmodel.settings.n_groups = size(Inputs.groups, 1)
        x_forced0 = Vector{Vector{Bool}}(undef, ATAmodel.settings.T)
        for t = 1:ATAmodel.settings.T
            x_forced0[t] = fill(true, ATAmodel.settings.n_items)
            ATAmodel.constraints[t].constr_A = zeros(Float64, 0, ATAmodel.settings.n_items)
        end
        #load bank
        if isfile(bank_file)
            ATAmodel.settings.bank =
                CSV.read(bank_file, DataFrames.DataFrame, delim = bank_delim)
            message[2] = message[2] * "- Item bank file read.\n"
        else
            push!(ATAmodel.output.infos, ["danger", "Not a valid file for bank."])
            return nothing
        end
        if !("OPT" in readdir())
            mkdir("OPT")
        end
        open("OPT/Settings.jl", "w") do f
            write(f, "#Settings \n\n")

            ATAmodel.settings.IRT.model = Inputs.IRT_model
            ATAmodel.settings.IRT.parameters = DataFrames.DataFrame(
                ATAmodel.settings.bank[!, Symbol.(Inputs.IRT_parameters)],
            )
            ATAmodel.settings.IRT.parametrization = Inputs.IRT_parametrization
            ATAmodel.settings.IRT.D = Inputs.IRT_D

            if ATAmodel.settings.IRT.model == "1PL"
                DataFrames.DataFrames.rename!(ATAmodel.settings.IRT.parameters, [:b])#nqp values in interval\r\n",
            elseif ATAmodel.settings.IRT.model == "2PL"
                DataFrames.DataFrames.rename!(ATAmodel.settings.IRT.parameters, [:a, :b]) #nqp values in interval\r\n",
            elseif ATAmodel.settings.IRT.model == "3PL"
                DataFrames.DataFrames.rename!(
                    ATAmodel.settings.IRT.parameters,
                    [:a, :b, :c],
                ) #nqp values in interval\r\n",
            else
                push!(
                    ATAmodel.output.infos,
                    ["danger", "Only 1PL, 2PL and 3PL IRT models are allowed."],
                )
                return nothing
            end
            CSV.write("OPT/IRT_parameters.csv", ATAmodel.settings.IRT.parameters)
            message[2] = message[2] * "- IRT item parameters loaded.\n"
            if Inputs.enemy_sets_var != String[]
                ATAmodel.settings.es.var = Symbol.(Inputs.enemy_sets_var)
                val = Symbol.(Inputs.enemy_sets_var)
                write(f, "enemy_sets_var = $val\n\n")
                message[2] = message[2] * "- Variable for Enemy sets loaded.\n"
            end
            if Inputs.friend_sets_var != String[]
                ATAmodel.settings.fs.var = Symbol.(Inputs.friend_sets_var)
                val = Symbol.(Inputs.friend_sets_var)
                write(f, "ATAmodel.settings.fs.var = $val\n\n")
                message[2] = message[2] * "- Variable for Friend sets loaded.\n"
            end

            if Inputs.length_min != Int64[]
                lengthmin = zeros(Int64, ATAmodel.settings.T)
                lengthweight = ones(Int64, ATAmodel.settings.T)
                t1 = 1
                for g = 1:ATAmodel.settings.n_groups
                    for t = 1:ATAmodel.settings.Tg[g]
                        lengthmin[t1] = Int(Inputs.length_min[g])
                        lengthweight[t1] = Inputs.length_weight[g]
                        t1 += 1
                    end
                end
                for t = 1:ATAmodel.settings.T
                    ATAmodel.constraints[t].constr_A = vcat(
                        ATAmodel.constraints[t].constr_A,
                        (-lengthweight[t]) .* ones(Float64, ATAmodel.settings.n_items)',
                    )
                    ATAmodel.constraints[t].constr_b = vcat(
                        ATAmodel.constraints[t].constr_b,
                        -lengthmin[t] * lengthweight[t],
                    )
                    ATAmodel.constraints[t].length_min = lengthmin[t]
                end
                write(f, "length_min = $lengthmin\n")
                message[2] = message[2] * "- Minimum length of tests constrained.\n"
            end

            if Inputs.length_max != Int64[]
                lengthmax = zeros(Int64, ATAmodel.settings.T)
                lengthweight = ones(Int64, ATAmodel.settings.T)
                t1 = 1
                for g = 1:ATAmodel.settings.n_groups
                    for t = 1:ATAmodel.settings.Tg[g]
                        lengthmax[t1] = Int(Inputs.length_max[g])
                        lengthweight[t1] = Inputs.length_weight[g]
                        t1 += 1
                    end
                end
                for t = 1:ATAmodel.settings.T
                    ATAmodel.constraints[t].constr_A = vcat(
                        ATAmodel.constraints[t].constr_A,
                        (lengthweight[t]) .* ones(ATAmodel.settings.n_items)',
                    )
                    ATAmodel.constraints[t].constr_b = vcat(
                        ATAmodel.constraints[t].constr_b,
                        lengthmax[t] * lengthweight[t],
                    )
                    ATAmodel.constraints[t].length_max = lengthmax[t]
                end
                write(f, "length_max = $lengthmax\n")
                message[2] = message[2] * "- Maximum length of tests constrained.\n"
            end

            if Inputs.expected_score_var != String[]
                t1 = 1
                for g = 1:ATAmodel.settings.n_groups
                    for t = 1:ATAmodel.settings.Tg[g]
                        ATAmodel.constraints[t1].expected_score.var =
                            Symbol.(Inputs.expected_score_var[g])
                        t1 += 1
                    end
                end
            end
            t1 = 1
            if Inputs.expected_score_pts != Float64[]
                for g = 1:ATAmodel.settings.n_groups
                    for t = 1:ATAmodel.settings.Tg[g]
                        ATAmodel.constraints[t1].expected_score.pts =
                            Inputs.expected_score_pts[g]
                        t1 += 1
                    end
                end
            end
            message[2] =
                message[2] * "- Expected score variable and expected_score_pts loaded.\n"

            if size(Inputs.expected_score_min, 1) > 0
                t1 = 1
                for g = 1:ATAmodel.settings.n_groups
                    for t = 1:ATAmodel.settings.Tg[g]
                        ATAmodel.constraints[t1].expected_score.min =
                            Inputs.expected_score_min[g]
                        t1 += 1
                    end
                end
                message[2] = message[2] * "- Minimum expected score constrained.\n"
            end
            if size(Inputs.expected_score_max, 1) > 0
                t1 = 1
                for g = 1:ATAmodel.settings.n_groups
                    for t = 1:ATAmodel.settings.Tg[g]
                        ATAmodel.constraints[t1].expected_score.max =
                            Inputs.expected_score_max[g]
                        t1 += 1
                    end
                end
                message[2] = message[2] * "- Maximum expected score constrained.\n"
            end
            if Inputs.sum_vars != Vector{Vector{String}}(undef, 0)
                t1 = 1
                for g = 1:ATAmodel.settings.n_groups
                    if !ismissing(sum_vars_min[g])
                        for v = 1:length(Inputs.sum_vars[g])
                            var = ATAmodel.settings.bank[Symbol(Inputs.sum_vars[g][v])]
                            var[ismissing.(var)] .= zero(Float64)
                            #min
                            if !ismissing(sum_vars_min[g][v])
                                for t = 1:ATAmodel.settings.Tg[g]
                                    if g > 1
                                        t1 = (g - 1) * Tg[g-1] + t
                                    else
                                        t1 = copy(t)
                                    end
                                    ATAmodel.constraints[t1].constr_A =
                                        vcat(ATAmodel.constraints[t1].constr_A, .-var')
                                    ATAmodel.constraints[t1].constr_b = vcat(
                                        ATAmodel.constraints[t1].constr_b,
                                        -sum_vars_min[g],
                                    )
                                end
                            end
                            #min
                            if !ismissing(sum_vars_max[g][v])
                                for t = 1:ATAmodel.settings.Tg[g]
                                    if g > 1
                                        t1 = (g - 1) * Tg[g-1] + t
                                    else
                                        t1 = copy(t)
                                    end
                                    ATAmodel.constraints[t1].constr_A =
                                        vcat(ATAmodel.constraints[t1].constr_A, var')
                                    ATAmodel.constraints[t1].constr_b = vcat(
                                        ATAmodel.constraints[t1].constr_b,
                                        sum_vars_max[g],
                                    )
                                end
                            end
                        end
                    end
                end
                message[2] =
                    message[2] *
                    string("- Sum of variables", Inputs.sum_vars, " constrained.\n")
            end
            if size(Inputs.item_use_min, 1) > 0
                ATAmodel.settings.iu.min = Inputs.item_use_min
                message[2] = message[2] * "- Minimum item use constrained.\n"
            end
            if size(Inputs.item_use_max, 1) > 0
                ATAmodel.settings.iu.max = Inputs.item_use_max
                for v = 1:ATAmodel.settings.T
                    x_forced0[v][findall(ATAmodel.settings.iu.max .< 1)] .= false
                end
                message[2] = message[2] * "- Maximum item use constrained.\n"
            end
            if Inputs.obj_type == "MAXIMIN"
                ATAmodel.obj.cores =
                    [MaximinObjectiveCore() for t = 1:sum(ATAmodel.settings.Tg)]
                if size(Inputs.obj_points, 1) > 0
                    t1 = 1
                    for g = 1:ATAmodel.settings.n_groups
                        for t = 1:ATAmodel.settings.Tg[g]
                            ATAmodel.obj.cores[t1].points = Inputs.obj_points[g]
                            t1 += 1
                        end
                    end
                    message[2] = message[2] * "- Optimization points loaded.\n"
                else
                    push!(
                        ATAmodel.output.infos,
                        [
                            "danger",
                            "error: MAXIMIN, CCMAXIMIN and MINIMAX objective types require optimization points. Use obj_points field in input settings.",
                        ],
                    )
                    return nothing
                end
            elseif Inputs.obj_type == "CCMAXIMIN"
                ATAmodel.obj.cores =
                    [CcMaximinObjectiveCore() for t = 1:sum(ATAmodel.settings.Tg)]
                if size(Inputs.obj_points, 1) > 0
                    t1 = 1
                    for g = 1:ATAmodel.settings.n_groups
                        for t = 1:ATAmodel.settings.Tg[g]
                            ATAmodel.obj.cores[t1].points = Inputs.obj_points[g]
                            ATAmodel.obj.cores[t1].alpha = Inputs.obj_aux_float
                            t1 += 1
                        end
                    end
                    message[2] = message[2] * "- Optimization points loaded.\n"
                else
                    push!(
                        ATAmodel.output.infos,
                        [
                            "danger",
                            "error: MAXIMIN, CCMAXIMIN and MINIMAX objective types require optimization points. Use obj_points field in input settings.",
                        ],
                    )
                    return nothing
                end
            elseif Inputs.obj_type == "MINIMAX"
                ATAmodel.obj.cores =
                    [MinimaxObjectiveCore() for t = 1:sum(ATAmodel.settings.Tg)]
                if size(Inputs.obj_points, 1) > 0
                    t1 = 1
                    for g = 1:ATAmodel.settings.n_groups
                        for t = 1:ATAmodel.settings.Tg[g]
                            ATAmodel.obj.cores[t1].points = Inputs.obj_points[g]
                            t1 += 1
                        end
                    end
                    message[2] = message[2] * "- Optimization points loaded.\n"
                else
                    push!(
                        ATAmodel.output.infos,
                        [
                            "danger",
                            "error: MAXIMIN, CCMAXIMIN and MINIMAX objective types require optimization points. Use obj_points field in input settings.",
                        ],
                    )
                    return nothing
                end
                if (size(Inputs.obj_targets, 1) > 0)
                    t1 = 1
                    for g = 1:ATAmodel.settings.n_groups
                        if size(Inputs.obj_targets[g], 1) != size(Inputs.obj_points[g], 1)
                            push!(
                                ATAmodel.output.infos,
                                [
                                    "danger",
                                    string(
                                        "error: for group ",
                                        g,
                                        " size of obj_targets (",
                                        size(Inputs.obj_targets[g], 1),
                                        ") is not the same as size of obj_points (",
                                        size(Inputs.obj_points[g], 1),
                                        "). /n",
                                    ),
                                ],
                            )
                            return nothing
                        end
                        for t = 1:ATAmodel.settings.Tg[g]
                            ATAmodel.obj.cores[t1].targets = Inputs.obj_targets[g]
                            t1 += 1
                        end
                    end
                    message[2] = message[2] * "- Targets loaded.\n"
                else
                    push!(
                        ATAmodel.output.infos,
                        [
                            "danger",
                            "error: MINIMAX objective type requires targets. Use obj_targets field in input settings.",
                        ],
                    )
                    return nothing
                end
            end
            #fictiuos friendSets
            ATAmodel.settings.fs.counts = ones(ATAmodel.settings.n_items)
            ATAmodel.settings.fs.sets = string.(collect(1:ATAmodel.settings.n_items))
            ATAmodel.settings.fs.items = [[i] for i = 1:ATAmodel.settings.n_items]
            ATAmodel.settings.forced0 = x_forced0
            if Inputs.categories != String[]
                val = Symbol.(Inputs.categories)
                ATAmodel.output.categories = copy(val)
                write(f, "categories = $val\n\n")
                message[2] = message[2] * "- Categories for output loaded.\n"
            end
        end
        message[2] =
            message[2] * string(
                "ASSEMBLE ",
                ATAmodel.settings.T,
                " FORMS DIVIDED INTO ",
                ATAmodel.settings.n_groups,
                " GROUPS.\n",
            )
        #update model
        JLD2.@save "OPT/ATAmodel.jld2" ATAmodel
    end
    message[1] = "success"
    push!(ATAmodel.output.infos, message)
    return ATAmodel
end
