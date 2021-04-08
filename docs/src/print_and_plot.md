# Print and Plot

```@meta
CurrentModule = ATA
```

## Print Results on Text File

```@docs
print_results(
        ata_model::NoObjModel;
        results_folder = "results",
        sim_pool::DataFrames.DataFrame = DataFrame(),
)
print_results(
        ata_model::MaximinModel;
        results_folder = "results",
        sim_pool::DataFrames.DataFrame = DataFrame(),
)
print_results(
        ata_model::MinimaxModel;
        results_folder = "results",
        sim_pool::DataFrames.DataFrame = DataFrame(),
)
print_results(
    ata_model::CCMaximinModel;
    results_folder = "results",
    sim_pool::DataFrames.DataFrame = DataFrame(),
)
```

## Plot TIFs and ICFs

```@docs
plot_results(
    ata_model::AbstractModel;
    plots_folder = "plots",
)
```