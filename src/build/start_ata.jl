"""
    start_ata(;
        settings::InputSettings = InputSettings(),
        bank::DataFrames.DataFrame = DataFrames.DataFrame(),
        settings_file = "settings_ata.jl",
        bank_file = "bank.csv",
        bank_delim = ";",
    )

# Description

Initialize an empty ATA model and load the test assembly settings contained in the `settings_file`
using item data in file `bank_file` which values are separated by `bank_delim`.
Alternatively, settings object and bank dataframe can be passed using the arguments `settings` and `bank`.

# Arguments

- **`settings::InputSettings`** : Optional. Default:  `InputSettings()` object with custom ATA settings.
- **`bank::DataFrames.DataFrame`** : Optional. Default: `DataFrame()`. A dataframe containing data about the items.
- **`settings_file`** : Optional. Default: "settings_ata.jl". The path of the file containing the ATA settings in the form of an `InputSettings` struct.
- **`bank_file`** : Optional. Default: "bank.csv". The path of the file containing the item pool/bank in the form of custom-separated values.
- **`bank_delim`** : Optional. Default: ";". The custom-separator for the bank_file.

# Output

- An ATA model.

"""
function start_ata(;
    settings::InputSettings = InputSettings(),
    bank::DataFrames.DataFrame = DataFrames.DataFrame(),
    settings_file = "settings_ata.jl",
    bank_file = "bank.csv",
    bank_delim = ";",
)
    message = ["", ""]
    ata_model = NoObjModel()
    Inputs = InputSettings()
    if size(settings.T, 1) > 0
        Inputs = settings
        message[2] = message[2] * "- InputSettings object loaded.\n"
    elseif !isfile(settings_file)
        message[1] = "danger"
        message[2] =
            message[2] * string(
                "Settings file with name ",
                settings_file,
                " does not exist. Provide a valid InputSettings object or a name of an existing file.\n",
            )
        push!(ata_model.output.infos, message)
        return ata_model
    else
        try
            Inputs = read_settings_file(settings_file)
        catch e
            message[1] = "danger"
            message[2] = message[2] * "- Error in reading settings file.\n"
            push!(ata_model.output.infos, message)
            return ata_model
        end
        message[2] = message[2] * "- Settings file loaded.\n"
    end
    try
        infos = ata_model.output.infos
        if Inputs.c
            if Inputs.obj_type == "maximin"
                ata_model = MaximinModel()
            elseif Inputs.obj_type == "minimax"
                ata_model = MinimaxModel()
            elseif Inputs.obj_type == "cc_maximin"
                ata_model = CCMaximinModel()
            elseif Inputs.obj_type == "soyster_maximin"
                ata_model = SoysterMaximinModel()
            elseif Inputs.obj_type == "de_jong_maximin"
                ata_model = DeJongMaximinModel()
            elseif Inputs.obj_type == "custom"
                ata_model = CustomModel()
            else
                error(
                    "Only \"maximin\", \"minimax\", \"cc_maximin\", \"soyster_maximin\", \"de_jong_maximin\", \"custom\" and \"\" objective types are supported.",
                )
            end
            message[2] =
                message[2] *
                string("- ", ata_model.obj.name, " optimization type loaded.\n")
        else
            ata_model = NoObjModel()
        end
        ata_model.output.infos = infos
        #load bank
        if size(bank, 1) > 0
            ata_model.settings.bank = bank
            message[2] = message[2] * "- Item bank dataframe loaded.\n"
        elseif isfile(bank_file)
            try
                ata_model.settings.bank =
                    CSV.read(bank_file, DataFrames.DataFrame, delim = bank_delim)
            catch e
                message[1] = "danger"
                message[2] = message[2] * "- Error in reading the item bank file.\n"
                return nothing
            end
            message[2] = message[2] * "- Item bank file loaded.\n"
        else
            push!(
                ata_model.output.infos,
                [
                    "danger",
                    "Item bank file with name ",
                    bank_file,
                    " does not exist. Provide a valid item bank dataframe or a name of an existing file.",
                ],
            )
            return nothing
        end
        ata_model.settings.n_groups = Inputs.n_groups
        ata_model.settings.n_items = Inputs.n_items
        ata_model.settings.T = Int(sum(Inputs.T))
        ata_model.settings.Tg = Inputs.T
        #initialize constraints
        ata_model.constraints = [Constraint() for t = 1:ata_model.settings.T]
        ata_model.settings.n_groups = size(Inputs.T, 1)
        x_forced0 = Vector{Vector{Bool}}(undef, ata_model.settings.T)
        for t = 1:ata_model.settings.T
            x_forced0[t] = fill(true, ata_model.settings.n_items)
            ata_model.constraints[t].constr_A =
                zeros(Float64, 0, ata_model.settings.n_items)
        end

        if !isdir("opt")
            mkdir("opt")
        end
        open("OPT/Settings.jl", "w") do f
            write(f, "#Settings \n\n")

            ata_model.settings.irt.model = Inputs.irt_model
            ata_model.settings.irt.parameters = DataFrames.DataFrame(
                ata_model.settings.bank[!, Symbol.(Inputs.irt_parameters)],
            )
            ata_model.settings.irt.parametrization = Inputs.irt_parametrization
            ata_model.settings.irt.D = Inputs.irt_D

            if ata_model.settings.irt.model == "1PL"
                DataFrames.rename!(ata_model.settings.irt.parameters, [:b])#nqp values in interval\r\n",
            elseif ata_model.settings.irt.model == "2PL"
                DataFrames.rename!(ata_model.settings.irt.parameters, [:a, :b]) #nqp values in interval\r\n",
            elseif ata_model.settings.irt.model == "3PL"
                DataFrames.rename!(
                    ata_model.settings.irt.parameters,
                    [:a, :b, :c],
                ) #nqp values in interval\r\n",
            else
                push!(
                    ata_model.output.infos,
                    ["danger", "Only 1PL, 2PL and 3PL IRT models are allowed."],
                )
                return nothing
            end
            CSV.write("OPT/irt_parameters.csv", ata_model.settings.irt.parameters)
            message[2] = message[2] * "- IRT item parameters loaded.\n"
                        #test length
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
            #friend sets
            #fictiuos friendSets
            ata_model.settings.fs.counts = ones(ata_model.settings.n_items)
            ata_model.settings.fs.sets = string.(collect(1:ata_model.settings.n_items))
            ata_model.settings.fs.items = [[i] for i = 1:ata_model.settings.n_items]
            ata_model.settings.forced0 = x_forced0
            if Inputs.friend_sets_var != String[]
                ata_model.settings.fs.var = Symbol.(Inputs.friend_sets_var)
                val = Symbol.(Inputs.friend_sets_var)
                _add_friends!(ata_model)
                write(f, "ata_model.settings.fs.var = $val\n\n")
                message[2] = message[2] * "- Variable for Friend sets loaded.\n"
            end
            #enemy sets
            if Inputs.enemy_sets_var != String[]
                ata_model.settings.es.var = Symbol.(Inputs.enemy_sets_var)
                val = Symbol.(Inputs.enemy_sets_var)
                _add_enemies!(ata_model)
                write(f, "enemy_sets_var = $val\n\n")
                message[2] = message[2] * "- Variable for Enemy sets loaded.\n"
            end

            #expected score
            if Inputs.expected_score_var != String[]
                t1 = 1
                for g = 1:ata_model.settings.n_groups
                    for t = 1:ata_model.settings.Tg[g]
                        ata_model.constraints[t1].expected_score.var =
                            Symbol.(Inputs.expected_score_var[g])
                        t1 += 1
                    end
                end
                message[2] =
                message[2] * "- Expected score variable loaded.\n"
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
            else
                for g = 1:ata_model.settings.n_groups
                    for t = 1:ata_model.settings.Tg[g]
                        ata_model.constraints[t1].expected_score.pts =
                            [0.0]
                        t1 += 1
                    end
                end
                message[2] =
                message[2] * "- Expected score points loaded.\n"
            end
            min_exp_score = false
            max_exp_score = false
            if size(Inputs.expected_score_min, 1) > 0
                if any(vcat(Inputs.expected_score_min...) .> 0)
                    min_exp_score = true
                    t1 = 1
                    for g = 1:ata_model.settings.n_groups
                        for t = 1:ata_model.settings.Tg[g]
                            ata_model.constraints[t1].expected_score.min =
                                Inputs.expected_score_min[g]
                            t1 += 1
                        end
                    end
                    message[2] = message[2] * "- Lower bounds for expected score loaded.\n"
                end
            end
            if size(Inputs.expected_score_max, 1) > 0
                if any(vcat(Inputs.expected_score_max...) .< 1.00)
                    max_exp_score = true
                    t1 = 1
                    for g = 1:ata_model.settings.n_groups
                        for t = 1:ata_model.settings.Tg[g]
                            ata_model.constraints[t1].expected_score.max =
                                Inputs.expected_score_max[g]
                            t1 += 1
                        end
                    end
                    message[2] = message[2] * "- Upper bounds for expected score loaded.\n"
                end
            end
            if max_exp_score || min_exp_score
                _add_exp_score!(ata_model)
                message[2] = message[2] * "- Expected scores computed.\n"
            end
            #sum vars
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
                if any(Inputs.item_use_min .> 0)
                    ata_model.settings.iu.min = Inputs.item_use_min
                    ata_model.settings.to_apply[2] = true
                    message[2] = message[2] * "- Minimum item use constrained.\n"
                end
            end
            #item use
            if size(Inputs.item_use_max, 1) > 0
                if any(Inputs.item_use_max .< ata_model.settings.T)
                    ata_model.settings.iu.max = Inputs.item_use_max
                    for v = 1:ata_model.settings.T
                        x_forced0[v][findall(ata_model.settings.iu.max .< 1)] .= false
                    end
                    ata_model.settings.to_apply[1] = true
                    message[2] = message[2] * "- Maximum item use constrained.\n"
                end
            end
            #obj_fun
            if Inputs.obj_type in ["maximin"]
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
                            "- maximin, soyster_maximin, de_jong_maximin, cc_maximin and minimax objective functions require optimization points. Use obj_points field in input settings.\n",
                        ],
                    )
                    return nothing
                end
            elseif Inputs.obj_type == "cc_maximin"
                ata_model.obj.cores =
                    [CCMaximinObjectiveCore() for t = 1:sum(ata_model.settings.Tg)]
                if size(Inputs.obj_points, 1) > 0
                    t1 = 1
                    for g = 1:ata_model.settings.n_groups
                        for t = 1:ata_model.settings.Tg[g]
                            ata_model.obj.cores[t1].points = Inputs.obj_points[g]
                            ata_model.obj.cores[t1].alpha = Inputs.obj_aux_float
                            t1 += 1
                        end
                    end
                    ata_model.obj.R = Inputs.obj_aux_int
                    message[2] = message[2] * "- Optimization points, number of replications, and alpha loaded.\n"
                else
                    push!(
                        ata_model.output.infos,
                        [
                            "danger",
                            "- maximin, cc_maximin, soyster_maximin, de_jong_maximin, and minimax objective functions require optimization points. Use obj_points field in input settings.\n",
                        ],
                    )
                    return nothing
                end
            elseif Inputs.obj_type in ["soyster_maximin", "de_jong_maximin"]
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
                    ata_model.obj.R = Inputs.obj_aux_int
                    message[2] = message[2] * "- Optimization points and number of replications loaded.\n"
                else
                    push!(
                        ata_model.output.infos,
                        [
                            "danger",
                            "- maximin, cc_maximin, soyster_maximin, de_jong_maximin, and minimax objective functions require optimization points. Use obj_points field in input settings.\n",
                        ],
                    )
                    return nothing
                end    
            elseif Inputs.obj_type == "minimax"
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
                            "- maximin, cc_maximin, soyster_maximin, de_jong_maximin, and minimax objective functions require optimization points. Use obj_points field in input settings.\n",
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
                                        "- For group ",
                                        g,
                                        " size of obj_targets (",
                                        size(Inputs.obj_targets[g], 1),
                                        ") is not the same as size of obj_points (",
                                        size(Inputs.obj_points[g], 1),
                                        ").\n",
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
                            "- minimax objective function requires targets. Use obj_targets field in input settings.\n",
                        ],
                    )
                    return nothing
                end
            end

            if Inputs.categories != String[]
                val = Symbol.(Inputs.categories)
                ata_model.output.categories = copy(val)
                write(f, "categories = $val\n\n")
                message[2] = message[2] * "- Categories for output loaded.\n"
            end

            #overlap matrix
            ata_model.settings.ol_max =
                zeros(Float64, ata_model.settings.T, ata_model.settings.T)
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
        message[1] = "success"
        push!(ata_model.output.infos, message)
    catch e
        message[1] = "danger"
        message[2] = message[2] * string("- ", sprint(showerror, e), "\n")
        push!(ata_model.output.infos, message)
    end
    return ata_model
end
