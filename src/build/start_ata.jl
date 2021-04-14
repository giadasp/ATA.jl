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

    ata_model = NoObjModel()

    try
        #LOAD INPUT SETTINGS, if error return ata_model
        Inputs = InputSettings()
        if size(settings.T, 1) > 0
            Inputs = settings
            success!(ata_model, "InputSettings object loaded.")
        elseif !isfile(settings_file)
            error!(
                ata_model,
                string(
                    "Settings file with name ",
                    settings_file,
                    " does not exist. Provide a valid InputSettings object or a name of an existing file.",
                ),
            )
        else
            try
                Inputs = read_settings_file(settings_file)
            catch e
                error!(
                    ata_model,
                    string(
                        "Error in reading the input settings file:\n  ",
                        sprint(showerror, e),
                        ".",
                    ),
                )
            end
            success!(ata_model, "Settings file loaded.")
        end

        #SELECT MODEL TYPE
        infos = ata_model.output.infos
        if Inputs.obj_type != ""
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
            elseif Inputs.obj_type == "robust_maximin"
                ata_model = RobustMaximinModel()
            elseif Inputs.obj_type == "custom"
                ata_model = CustomModel()
            else
                error!(ata_model,
                    "Only \"maximin\", \"minimax\", \"cc_maximin\", \"soyster_maximin\", \"de_jong_maximin\", \"custom\" and \"\" objective types are supported."
                )
            end
            success!(ata_model, string(ata_model.obj.name, " model initialized."))
        end
        ata_model.output.infos = infos

        #START ASSIGNING SETTINGS
        if Inputs.n_groups > 0
            ata_model.settings.n_groups = Inputs.n_groups
        else
            success!(ata_model, "n_groups has been set to 1.")
        end
        if Inputs.n_items > 0
            ata_model.settings.n_items = Inputs.n_items
        else
            error!(ata_model, "n_items must be grater than 0.")
        end
        if Int(sum(Inputs.T)) > 0
            ata_model.settings.T = Int(sum(Inputs.T))
            ata_model.settings.Tg = Inputs.T
        else
            error!(ata_model, "Sum of T must be grater than 0.")
        end

        #INITIALIZE CONSTRAINTS
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

        #ADD BANK
        add_bank!(ata_model, bank = bank, bank_file = bank_file, bank_delim = bank_delim)
        #open("opt/Settings.jl", "w") do f
        #write(f, "#Settings \n")

        #IRT
        if (Inputs.irt_model != "")
            ata_model.settings.irt.model = Inputs.irt_model
            ata_model.settings.irt.parametrization = Inputs.irt_parametrization
            ata_model.settings.irt.D = Inputs.irt_D
            if (size(ata_model.settings.bank, 1) == ata_model.settings.n_items)
                _add_irt!(ata_model)
            end
        end

        #TEST LENGTH
        if Inputs.length_min != Int64[]
            lengthmin = zeros(Int64, ata_model.settings.T)
            lengthweight = ones(Int64, ata_model.settings.T)
            t1 = 1
            for g = 1:ata_model.settings.n_groups, t = 1:ata_model.settings.Tg[g]
                lengthmin[t1] = Int(Inputs.length_min[g])
                lengthweight[t1] = Inputs.length_weight[g]
                t1 += 1
            end
            for t = 1:ata_model.settings.T
                ata_model.constraints[t].constr_A = vcat(
                    ata_model.constraints[t].constr_A,
                    (-lengthweight[t]) .* ones(Float64, ata_model.settings.n_items)',
                )
                ata_model.constraints[t].constr_b =
                    vcat(ata_model.constraints[t].constr_b, -lengthmin[t] * lengthweight[t])
                ata_model.constraints[t].length_min = lengthmin[t]
            end
            #write(f, "length_min = $lengthmin")
            success!(ata_model, "Minimum length of tests constrained.")
        end

        if Inputs.length_max != Int64[]
            lengthmax = zeros(Int64, ata_model.settings.T)
            lengthweight = ones(Int64, ata_model.settings.T)
            t1 = 1
            for g = 1:ata_model.settings.n_groups, t = 1:ata_model.settings.Tg[g]
                lengthmax[t1] = Int(Inputs.length_max[g])
                lengthweight[t1] = Inputs.length_weight[g]
                t1 += 1
            end
            for t = 1:ata_model.settings.T
                ata_model.constraints[t].constr_A = vcat(
                    ata_model.constraints[t].constr_A,
                    (lengthweight[t]) .* ones(ata_model.settings.n_items)',
                )
                ata_model.constraints[t].constr_b =
                    vcat(ata_model.constraints[t].constr_b, lengthmax[t] * lengthweight[t])
                ata_model.constraints[t].length_max = lengthmax[t]
            end
            #write(f, "length_max = $lengthmax")
            success!(ata_model, "Maximum length of tests constrained.")
        end


        #FICTIOUS FRIEND SETS
        ata_model.settings.n_fs = ata_model.settings.n_items
        ata_model.settings.fs.counts = ones(ata_model.settings.n_items)
        ata_model.settings.fs.sets = string.(collect(1:ata_model.settings.n_items))
        ata_model.settings.fs.items = [[i] for i = 1:ata_model.settings.n_items]
        ata_model.settings.forced0 = x_forced0


        if Inputs.friend_sets_var != String[]
            ata_model.settings.fs.var = Symbol.(Inputs.friend_sets_var)
            if size(ata_model.settings.bank, 1) == ata_model.settings.n_items
                val = Symbol.(ata_model.settings.fs.var)
                _add_friends!(ata_model)
            else
                error!(
                    ata_model,
                    string(
                        "Variable for Friend sets NOT loaded because the provided item bank has length not equal to ",
                        ata_model.settings.n_items,
                        ".",
                    ),
                )
            end
        end
        #enemy sets
        if Inputs.enemy_sets_var != String[]
            ata_model.settings.es.var = Symbol.(Inputs.enemy_sets_var)
            if size(ata_model.settings.bank, 1) == ata_model.settings.n_items
                val = Symbol.(Inputs.enemy_sets_var)
                _add_enemies!(ata_model)
            else
                error!(
                    ata_model,
                    string(
                        "Variable for Enemy sets NOT loaded because the provided item bank has length not equal to ",
                        ata_model.settings.n_items,
                        ".",
                    ),
                )
            end
        end

        #EXPECTED SCORE
        if Inputs.expected_score_var != String[]
            t1 = 1
            for g = 1:ata_model.settings.n_groups, t = 1:ata_model.settings.Tg[g]
                ata_model.constraints[t1].expected_score.var =
                    Symbol.(Inputs.expected_score_var[g])
                t1 += 1
            end
            success!(ata_model, "Expected score variable loaded.")
        end
        if Inputs.expected_score_pts != Float64[]
            t1 = 1
            for g = 1:ata_model.settings.n_groups, t = 1:ata_model.settings.Tg[g]
                ata_model.constraints[t1].expected_score.pts = Inputs.expected_score_pts[g]
                t1 += 1
            end
        else
            t1 = 1
            for g = 1:ata_model.settings.n_groups, t = 1:ata_model.settings.Tg[g]
                ata_model.constraints[t1].expected_score.pts = [0.0]
                t1 += 1
            end
            success!(ata_model, "Expected score points loaded.")
        end
        min_exp_score = false
        max_exp_score = false
        if size(Inputs.expected_score_min, 1) > 0
            if any(vcat(Inputs.expected_score_min...) .> 0)
                min_exp_score = true
                t1 = 1
                for g = 1:ata_model.settings.n_groups, t = 1:ata_model.settings.Tg[g]
                    ata_model.constraints[t1].expected_score.min =
                        Inputs.expected_score_min[g]
                    t1 += 1
                end
                success!(ata_model, "Lower bounds for expected score loaded.")
            end
        end
        if size(Inputs.expected_score_max, 1) > 0
            if any(vcat(Inputs.expected_score_max...) .< 1.00)
                max_exp_score = true
                t1 = 1
                for g = 1:ata_model.settings.n_groups, t = 1:ata_model.settings.Tg[g]
                    ata_model.constraints[t1].expected_score.max =
                        Inputs.expected_score_max[g]
                    t1 += 1
                end
                success!(ata_model, "Upper bounds for expected score loaded.")
            end
        end
        if (max_exp_score || min_exp_score)
            if size(ata_model.settings.bank, 1) == ata_model.settings.n_items
                _add_exp_score!(ata_model)
            else
                error!(
                    ata_model,
                    string(
                        "Expected scores (ICFs) based on IRT parameters NOT computed because the provided item bank has length not equal to ",
                        ata_model.settings.n_items,
                        ".",
                    ),
                )
            end
        end

        #ITEM USE
        if size(Inputs.item_use_min, 1) > 0
            if any(Inputs.item_use_min .> 0)
                ata_model.settings.iu.min = Inputs.item_use_min
                ata_model.settings.to_apply[2] = true
                success!(ata_model, "Minimum item use constrained.")
            end
        end
        if size(Inputs.item_use_max, 1) > 0
            if any(Inputs.item_use_max .< ata_model.settings.T)
                ata_model.settings.iu.max = Inputs.item_use_max
                for v = 1:ata_model.settings.T
                    x_forced0[v][findall(ata_model.settings.iu.max .< 1)] .= false
                end
                ata_model.settings.to_apply[1] = true
                success!(ata_model, "Maximum item use constrained.")
            end
        end

        #OBJECTIVE FUNCTION
        if Inputs.obj_type in ["maximin"]
            ata_model.obj.cores =
                [MaximinObjectiveCore() for t = 1:sum(ata_model.settings.Tg)]
            if size(Inputs.obj_points, 1) > 0
                t1 = 1
                for g = 1:ata_model.settings.n_groups, t = 1:ata_model.settings.Tg[g]
                    ata_model.obj.cores[t1].points = Inputs.obj_points[g]
                    t1 += 1
                end
                success!(ata_model, "Optimization points loaded.")
            else
                error!(
                    ata_model,
                    "maximin, cc_maximin, soyster_maximin, de_jong_maximin, robust_maximin, and minimax objective functions require optimization points. Use obj_points field in input settings.",
                )
                return nothing
            end
        elseif Inputs.obj_type == "cc_maximin"
            ata_model.obj.cores =
                [CCMaximinObjectiveCore() for t = 1:sum(ata_model.settings.Tg)]
            if size(Inputs.obj_points, 1) > 0
                t1 = 1
                for g = 1:ata_model.settings.n_groups, t = 1:ata_model.settings.Tg[g]
                    ata_model.obj.cores[t1].points = Inputs.obj_points[g]
                    ata_model.obj.cores[t1].alpha = Inputs.obj_aux_float
                    t1 += 1
                end
                ata_model.obj.R = Inputs.obj_aux_int
                success!(
                    ata_model,
                    string(
                        "Optimization points, number of samples (",
                        ata_model.obj.R,
                        "), and alpha (",
                        Inputs.obj_aux_float,
                        ") loaded.",
                    ),
                )
            else
                error!(
                    ata_model,
                    "maximin, cc_maximin, soyster_maximin, de_jong_maximin, robust_maximin, and minimax objective functions require optimization points. Use obj_points field in input settings.",
                )
                return nothing
            end
        elseif Inputs.obj_type in ["soyster_maximin", "de_jong_maximin"]
            ata_model.obj.cores =
                [MaximinObjectiveCore() for t = 1:sum(ata_model.settings.Tg)]
            if size(Inputs.obj_points, 1) > 0
                t1 = 1
                for g = 1:ata_model.settings.n_groups, t = 1:ata_model.settings.Tg[g]
                    ata_model.obj.cores[t1].points = Inputs.obj_points[g]
                    t1 += 1
                end
                ata_model.obj.R = Inputs.obj_aux_int
                success!(
                    ata_model,
                    "Optimization points and number of replications loaded.",
                )
            else
                error!(
                    ata_model,
                    "maximin, cc_maximin, soyster_maximin, de_jong_maximin, robust_maximin, and minimax objective functions require optimization points. Use obj_points field in input settings.",
                )
                return nothing
            end
        elseif Inputs.obj_type in ["robust_maximin"]
            ata_model.obj.cores =
                [RobustMaximinObjectiveCore() for t = 1:sum(ata_model.settings.Tg)]
            if size(Inputs.obj_points, 1) > 0
                t1 = 1
                for g = 1:ata_model.settings.n_groups, t = 1:ata_model.settings.Tg[g]
                    ata_model.obj.cores[t1].points = Inputs.obj_points[g]
                    t1 += 1
                end
                ata_model.obj.R = Inputs.obj_aux_int
                ata_model.obj.Gamma = Inputs.obj_aux_float
                success!(
                    ata_model,
                    string(
                        "Optimization points, number of samples (",
                        ata_model.obj.R,
                        "), and Gamma (",
                        ata_model.obj.Gamma,
                        ") loaded.",
                    ),
                )
            else
                error!(
                    ata_model,
                    "maximin, cc_maximin, soyster_maximin, de_jong_maximin, robust_maximin, and minimax objective functions require optimization points. Use obj_points field in input settings.",
                )
                return nothing
            end
        elseif Inputs.obj_type == "minimax"
            ata_model.obj.cores =
                [MinimaxObjectiveCore() for t = 1:sum(ata_model.settings.Tg)]
            if size(Inputs.obj_points, 1) > 0
                t1 = 1
                for g = 1:ata_model.settings.n_groups, t = 1:ata_model.settings.Tg[g]
                    ata_model.obj.cores[t1].points = Inputs.obj_points[g]
                    t1 += 1
                end
                success!(ata_model, "Optimization points loaded.")
            else
                error!(
                    ata_model,
                    "maximin, cc_maximin, soyster_maximin, de_jong_maximin, robust_maximin, and minimax objective functions require optimization points. Use obj_points field in input settings.",
                )
                return nothing
            end
            if (size(Inputs.obj_targets, 1) > 0)
                t1 = 1
                for g = 1:ata_model.settings.n_groups
                    if size(Inputs.obj_targets[g], 1) != size(Inputs.obj_points[g], 1)
                        error!(
                            ata_model,
                            string(
                                "For group ",
                                g,
                                " size of obj_targets (",
                                size(Inputs.obj_targets[g], 1),
                                ") is not the same as size of obj_points (",
                                size(Inputs.obj_points[g], 1),
                                ").",
                            ),
                        )
                        return nothing
                    end
                    for t = 1:ata_model.settings.Tg[g]
                        ata_model.obj.cores[t1].targets = Inputs.obj_targets[g]
                        t1 += 1
                    end
                end
                success!(ata_model, "Targets loaded.")
            else
                error!(
                    ata_model,
                    "minimax objective function requires targets. Use obj_targets field in input settings.",
                )
                return nothing
            end
        end

        if Inputs.categories != String[]
            val = Symbol.(Inputs.categories)
            ata_model.output.categories = copy(val)
            #write(f, "categories = $val\n")
            success!(ata_model, "Categories for output loaded.")
        end

        #overlap matrix
        ata_model.settings.ol_max =
            zeros(Float64, ata_model.settings.T, ata_model.settings.T)
        #end
        success!(
            ata_model,
            string(
                "ASSEMBLE ",
                ata_model.settings.T,
                " FORMS DIVIDED INTO ",
                ata_model.settings.n_groups,
                " GROUPS.",
            ),
        )
        #update model
        JLD2.@save "opt/ata_model.jld2" ata_model
    catch e
        error!(ata_model, string(sprint(showerror, e)))
    end
    return ata_model
end
