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

- **`ata_model::AbstractModel`** : Required. The model built with `start_ATA()` function.
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
                ata_model = MaximinModel()
            elseif Inputs.obj_type == "MINIMAX"
                ata_model = MinimaxModel()
            elseif Inputs.obj_type == "CCMAXIMIN"
                ata_model = CcMaximinModel()
            elseif Inputs.obj_type == "custom"
                ata_model = CustomModel()
            else
                error(
                    "Only \"MAXIMIN\", \"MINIMAX\", \"CCMAXIMIN\", \"custom\" and \"\" objective types are supported.",
                )
            end
            message[2] =
                message[2] * string("- ", ata_model.obj.name, " optimization type loaded.\n")
        else
            ata_model = NoObjModel()
        end
        ata_model.settings.n_items = Inputs.n_items
        ata_model.settings.T = Int(sum(Inputs.T))
        ata_model.settings.Tg = Inputs.T
        #initialize constraints
        ata_model.constraints = [Constraint() for t = 1:ata_model.settings.T]
        ata_model.settings.n_groups = size(Inputs.groups, 1)
        x_forced0 = Vector{Vector{Bool}}(undef, ata_model.settings.T)
        for t = 1:ata_model.settings.T
            x_forced0[t] = fill(true, ata_model.settings.n_items)
            ata_model.constraints[t].constr_A = zeros(Float64, 0, ata_model.settings.n_items)
        end
        #load bank
        if isfile(bank_file)
            ata_model.settings.bank =
                CSV.read(bank_file, DataFrames.DataFrame, delim = bank_delim)
            message[2] = message[2] * "- Item bank file read.\n"
        else
            push!(ata_model.output.infos, ["danger", "Not a valid file for bank."])
            return nothing
        end
        if !("OPT" in readdir())
            mkdir("OPT")
        end
        open("OPT/Settings.jl", "w") do f
            write(f, "#Settings \n\n")

            ata_model.settings.IRT.model = Inputs.IRT_model
            ata_model.settings.IRT.parameters = DataFrames.DataFrame(
                ata_model.settings.bank[!, Symbol.(Inputs.IRT_parameters)],
            )
            ata_model.settings.IRT.parametrization = Inputs.IRT_parametrization
            ata_model.settings.IRT.D = Inputs.IRT_D

            if ata_model.settings.IRT.model == "1PL"
                DataFrames.DataFrames.rename!(ata_model.settings.IRT.parameters, [:b])#nqp values in interval\r\n",
            elseif ata_model.settings.IRT.model == "2PL"
                DataFrames.DataFrames.rename!(ata_model.settings.IRT.parameters, [:a, :b]) #nqp values in interval\r\n",
            elseif ata_model.settings.IRT.model == "3PL"
                DataFrames.DataFrames.rename!(
                    ata_model.settings.IRT.parameters,
                    [:a, :b, :c],
                ) #nqp values in interval\r\n",
            else
                push!(
                    ata_model.output.infos,
                    ["danger", "Only 1PL, 2PL and 3PL IRT models are allowed."],
                )
                return nothing
            end
            CSV.write("OPT/IRT_parameters.csv", ata_model.settings.IRT.parameters)
            message[2] = message[2] * "- IRT item parameters loaded.\n"
            if Inputs.enemy_sets_var != String[]
                ata_model.settings.es.var = Symbol.(Inputs.enemy_sets_var)
                val = Symbol.(Inputs.enemy_sets_var)
                write(f, "enemy_sets_var = $val\n\n")
                message[2] = message[2] * "- Variable for Enemy sets loaded.\n"
            end
            if Inputs.friend_sets_var != String[]
                ata_model.settings.fs.var = Symbol.(Inputs.friend_sets_var)
                val = Symbol.(Inputs.friend_sets_var)
                write(f, "ata_model.settings.fs.var = $val\n\n")
                message[2] = message[2] * "- Variable for Friend sets loaded.\n"
            end

            if Inputs.length_min != Int64[]
                lengthmin = zeros(Int64, ata_model.settings.T)
                lengthweight = ones(Int64, ata_model.settings.T)
                t1 = 1
                for g = 1:ata_model.settings.n_groups
                    for t = 1:ata_model.settings.Tg[g]
                        lengthmin[t1] = Int(Inputs.length_min[g])
                        lengthweight[t1] = Inputs.length_weight[g]
                        t1 += 1
                    end
                end
                for t = 1:ata_model.settings.T
                    ata_model.constraints[t].constr_A = vcat(
                        ata_model.constraints[t].constr_A,
                        (-lengthweight[t]) .* ones(Float64, ata_model.settings.n_items)',
                    )
                    ata_model.constraints[t].constr_b = vcat(
                        ata_model.constraints[t].constr_b,
                        -lengthmin[t] * lengthweight[t],
                    )
                    ata_model.constraints[t].length_min = lengthmin[t]
                end
                write(f, "length_min = $lengthmin\n")
                message[2] = message[2] * "- Minimum length of tests constrained.\n"
            end

            if Inputs.length_max != Int64[]
                lengthmax = zeros(Int64, ata_model.settings.T)
                lengthweight = ones(Int64, ata_model.settings.T)
                t1 = 1
                for g = 1:ata_model.settings.n_groups
                    for t = 1:ata_model.settings.Tg[g]
                        lengthmax[t1] = Int(Inputs.length_max[g])
                        lengthweight[t1] = Inputs.length_weight[g]
                        t1 += 1
                    end
                end
                for t = 1:ata_model.settings.T
                    ata_model.constraints[t].constr_A = vcat(
                        ata_model.constraints[t].constr_A,
                        (lengthweight[t]) .* ones(ata_model.settings.n_items)',
                    )
                    ata_model.constraints[t].constr_b = vcat(
                        ata_model.constraints[t].constr_b,
                        lengthmax[t] * lengthweight[t],
                    )
                    ata_model.constraints[t].length_max = lengthmax[t]
                end
                write(f, "length_max = $lengthmax\n")
                message[2] = message[2] * "- Maximum length of tests constrained.\n"
            end

            if Inputs.expected_score_var != String[]
                t1 = 1
                for g = 1:ata_model.settings.n_groups
                    for t = 1:ata_model.settings.Tg[g]
                        ata_model.constraints[t1].expected_score.var =
                            Symbol.(Inputs.expected_score_var[g])
                        t1 += 1
                    end
                end
            end
            t1 = 1
            if Inputs.expected_score_pts != Float64[]
                for g = 1:ata_model.settings.n_groups
                    for t = 1:ata_model.settings.Tg[g]
                        ata_model.constraints[t1].expected_score.pts =
                            Inputs.expected_score_pts[g]
                        t1 += 1
                    end
                end
            end
            message[2] =
                message[2] * "- Expected score variable and expected_score_pts loaded.\n"

            if size(Inputs.expected_score_min, 1) > 0
                t1 = 1
                for g = 1:ata_model.settings.n_groups
                    for t = 1:ata_model.settings.Tg[g]
                        ata_model.constraints[t1].expected_score.min =
                            Inputs.expected_score_min[g]
                        t1 += 1
                    end
                end
                message[2] = message[2] * "- Minimum expected score constrained.\n"
            end
            if size(Inputs.expected_score_max, 1) > 0
                t1 = 1
                for g = 1:ata_model.settings.n_groups
                    for t = 1:ata_model.settings.Tg[g]
                        ata_model.constraints[t1].expected_score.max =
                            Inputs.expected_score_max[g]
                        t1 += 1
                    end
                end
                message[2] = message[2] * "- Maximum expected score constrained.\n"
            end
            if Inputs.sum_vars != Vector{Vector{String}}(undef, 0)
                t1 = 1
                for g = 1:ata_model.settings.n_groups
                    if !ismissing(sum_vars_min[g])
                        for v = 1:length(Inputs.sum_vars[g])
                            var = ata_model.settings.bank[Symbol(Inputs.sum_vars[g][v])]
                            var[ismissing.(var)] .= zero(Float64)
                            #min
                            if !ismissing(sum_vars_min[g][v])
                                for t = 1:ata_model.settings.Tg[g]
                                    if g > 1
                                        t1 = (g - 1) * Tg[g-1] + t
                                    else
                                        t1 = copy(t)
                                    end
                                    ata_model.constraints[t1].constr_A =
                                        vcat(ata_model.constraints[t1].constr_A, .-var')
                                    ata_model.constraints[t1].constr_b = vcat(
                                        ata_model.constraints[t1].constr_b,
                                        -sum_vars_min[g],
                                    )
                                end
                            end
                            #min
                            if !ismissing(sum_vars_max[g][v])
                                for t = 1:ata_model.settings.Tg[g]
                                    if g > 1
                                        t1 = (g - 1) * Tg[g-1] + t
                                    else
                                        t1 = copy(t)
                                    end
                                    ata_model.constraints[t1].constr_A =
                                        vcat(ata_model.constraints[t1].constr_A, var')
                                    ata_model.constraints[t1].constr_b = vcat(
                                        ata_model.constraints[t1].constr_b,
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
                ata_model.settings.iu.min = Inputs.item_use_min
                message[2] = message[2] * "- Minimum item use constrained.\n"
            end
            if size(Inputs.item_use_max, 1) > 0
                ata_model.settings.iu.max = Inputs.item_use_max
                for v = 1:ata_model.settings.T
                    x_forced0[v][findall(ata_model.settings.iu.max .< 1)] .= false
                end
                message[2] = message[2] * "- Maximum item use constrained.\n"
            end
            if Inputs.obj_type == "MAXIMIN"
                ata_model.obj.cores =
                    [MaximinObjectiveCore() for t = 1:sum(ata_model.settings.Tg)]
                if size(Inputs.obj_points, 1) > 0
                    t1 = 1
                    for g = 1:ata_model.settings.n_groups
                        for t = 1:ata_model.settings.Tg[g]
                            ata_model.obj.cores[t1].points = Inputs.obj_points[g]
                            t1 += 1
                        end
                    end
                    message[2] = message[2] * "- Optimization points loaded.\n"
                else
                    push!(
                        ata_model.output.infos,
                        [
                            "danger",
                            "error: MAXIMIN, CCMAXIMIN and MINIMAX objective types require optimization points. Use obj_points field in input settings.",
                        ],
                    )
                    return nothing
                end
            elseif Inputs.obj_type == "CCMAXIMIN"
                ata_model.obj.cores =
                    [CcMaximinObjectiveCore() for t = 1:sum(ata_model.settings.Tg)]
                if size(Inputs.obj_points, 1) > 0
                    t1 = 1
                    for g = 1:ata_model.settings.n_groups
                        for t = 1:ata_model.settings.Tg[g]
                            ata_model.obj.cores[t1].points = Inputs.obj_points[g]
                            ata_model.obj.cores[t1].alpha = Inputs.obj_aux_float
                            t1 += 1
                        end
                    end
                    message[2] = message[2] * "- Optimization points loaded.\n"
                else
                    push!(
                        ata_model.output.infos,
                        [
                            "danger",
                            "error: MAXIMIN, CCMAXIMIN and MINIMAX objective types require optimization points. Use obj_points field in input settings.",
                        ],
                    )
                    return nothing
                end
            elseif Inputs.obj_type == "MINIMAX"
                ata_model.obj.cores =
                    [MinimaxObjectiveCore() for t = 1:sum(ata_model.settings.Tg)]
                if size(Inputs.obj_points, 1) > 0
                    t1 = 1
                    for g = 1:ata_model.settings.n_groups
                        for t = 1:ata_model.settings.Tg[g]
                            ata_model.obj.cores[t1].points = Inputs.obj_points[g]
                            t1 += 1
                        end
                    end
                    message[2] = message[2] * "- Optimization points loaded.\n"
                else
                    push!(
                        ata_model.output.infos,
                        [
                            "danger",
                            "error: MAXIMIN, CCMAXIMIN and MINIMAX objective types require optimization points. Use obj_points field in input settings.",
                        ],
                    )
                    return nothing
                end
                if (size(Inputs.obj_targets, 1) > 0)
                    t1 = 1
                    for g = 1:ata_model.settings.n_groups
                        if size(Inputs.obj_targets[g], 1) != size(Inputs.obj_points[g], 1)
                            push!(
                                ata_model.output.infos,
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
                        for t = 1:ata_model.settings.Tg[g]
                            ata_model.obj.cores[t1].targets = Inputs.obj_targets[g]
                            t1 += 1
                        end
                    end
                    message[2] = message[2] * "- Targets loaded.\n"
                else
                    push!(
                        ata_model.output.infos,
                        [
                            "danger",
                            "error: MINIMAX objective type requires targets. Use obj_targets field in input settings.",
                        ],
                    )
                    return nothing
                end
            end
            #fictiuos friendSets
            ata_model.settings.fs.counts = ones(ata_model.settings.n_items)
            ata_model.settings.fs.sets = string.(collect(1:ata_model.settings.n_items))
            ata_model.settings.fs.items = [[i] for i = 1:ata_model.settings.n_items]
            ata_model.settings.forced0 = x_forced0
            if Inputs.categories != String[]
                val = Symbol.(Inputs.categories)
                ata_model.output.categories = copy(val)
                write(f, "categories = $val\n\n")
                message[2] = message[2] * "- Categories for output loaded.\n"
            end
        end
        message[2] =
            message[2] * string(
                "ASSEMBLE ",
                ata_model.settings.T,
                " FORMS DIVIDED INTO ",
                ata_model.settings.n_groups,
                " GROUPS.\n",
            )
        #update model
        JLD2.@save "OPT/ata_model.jld2" ata_model
    end
    message[1] = "success"
    push!(ata_model.output.infos, message)
    return ata_model
end