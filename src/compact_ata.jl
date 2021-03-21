function compact_ata(;
    settings::InputSettings = InputSettings(),
    bank::DataFrames.DataFrame = DataFrames.DataFrame(),
    settings_file = "",
    bank_file = "",
    bank_delim = ";",
    add_friends = true,
    add_enemies = true,
    add_constraints = true,
    constraints::DataFrames.DataFrame = DataFrames.DataFrame(),
    constraints_file = "",
    constraints_delim = ";",
    add_overlap = true,
    overlap::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),
    overlap_file = "",
    overlap_delim = ";",
    add_exp_score = true,
    group_by_friends = true,
    add_obj_fun = true,
    solver = "siman", #"jumpATA",
    print_it = true,
    print_folder = "RESULTS",
    plot_it = true,
    plot_folder = "PLOTS",
    kwargs...,
)
    ata_model = start_ata(;
        settings = settings,
        bank = bank,
        settings_file = settings_file,
        bank_file = bank_file,
        bank_delim = bank_delim,
    )
    print_last_info(ata_model)
    if add_friends
        add_friends!(ata_model)
        print_last_info(ata_model)
    end
    if add_enemies
        add_enemies!(ata_model)
        print_last_info(ata_model)
    end
    if add_constraints
        add_constraints!(
            ata_model;
            constraints = constraints,
            constraints_file = constraints_file,
            constraints_delim = constraints_delim,
        )
        print_last_info(ata_model)
    end
    if add_overlap
        add_overlap!(
            ata_model;
            overlap = overlap,
            overlap_file = overlap_file,
            overlap_delim = overlap_delim
        )
        print_last_info(ata_model)
    end
    if add_exp_score
        add_exp_score!(ata_model)
        print_last_info(ata_model)
    end
    if group_by_friends
        group_by_friends!(ata_model)
        print_last_info(ata_model)
    end
    if add_obj_fun
        add_obj_fun!(ata_model)
        print_last_info(ata_model)
    end
    # 9. assemble
    assemble!(ata_model; solver = solver, kwargs...)
    if print_it
        print_results(
            ata_model;
            group_by_fs = group_by_friends,
            results_folder = print_folder,
        )
    end
    if plot_it
        plot_results(
            ata_model;
            group_by_fs = group_by_friends,
            results_folder = plot_folder,
        )
    end
    return ata_model::AbstractModel
end
