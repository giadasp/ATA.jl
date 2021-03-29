# Print and Plot

```@meta
CurrentModule = ATA
```

## Print Results on Text File

```@docs
print_results(
        ata_model::NoObjModel;
        group_by_fs = false,
        results_folder = "results",
        sim_pool::DataFrames.DataFrame = DataFrame(),
)
print_results(
        ata_model::MaximinModel;
        group_by_fs = false,
        results_folder = "results",
        sim_pool::DataFrames.DataFrame = DataFrame(),
)
print_results(
        ata_model::MinimaxModel;
        group_by_fs = false,
        results_folder = "results",
        sim_pool::DataFrames.DataFrame = DataFrame(),
)
print_results(
    ata_model::CCMaximinModel;
    group_by_fs = false,
    results_folder = "results",
    sim_pool::DataFrames.DataFrame = DataFrame(),
)
```

## Plot TIFs and ICFs

```@docs
plot_results(
    ata_model::AbstractModel;
    group_by_fs = false,
    results_folder = "plots",
)
```