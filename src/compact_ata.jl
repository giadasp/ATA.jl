function compact_ata(;
    settings::InputSettings = InputSettings(),
    bank::DataFrames.DataFrame = DataFrames.DataFrame(),
    settings_file = "",
    bank_file = "",
    bank_delim = ";",
    add_constraints = true,
    constraints::DataFrames.DataFrame = DataFrames.DataFrame(),
    constraints_file = "",
    constraints_delim = ";",
    add_overlap = true,
    overlap::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),
    overlap_file = "",
    overlap_delim = ";",
    add_obj_fun = true,
    solver = "siman", #"jump",
    return_results = true,
    results_folder = "results",
    return_plots = true,
    plots_folder = "plots",
    sim_pool = DataFrames.DataFrame(),
    psychometrics = false,
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
            overlap_delim = overlap_delim,
        )
        print_last_info(ata_model)
    end
    if add_obj_fun
        add_obj_fun!(ata_model; psychometrics = psychometrics)
        print_last_info(ata_model)
    end
    # 9. assemble
    assemble!(ata_model; solver = solver, kwargs...)
    if return_results
        print_results(ata_model; results_folder = results_folder, sim_pool = sim_pool)
    end
    if return_plots
        plot_results(ata_model; plots_folder = plots_folder)
    end
    return ata_model::AbstractModel
end
