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
    ata_model::Union{MaximinModel, SoysterMaximinModel, DeJongMaximinModel, RobustMaximinModel, MinimaxModel};
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

To plot the TIFs and ICTs you must install (add) and load (using) the Julia package [`ATAPlot`](https://github.com/giadasp/ATAPlot.jl).
It requires to have the Julia package `PGFPlotsX.jl` and Miktex installed. 
A new version of `ATAPlot.jl` which use the package `Plots` and the backend `GR` is under development. 

```@docs
plot_results(
    ata_model::Union{MaximinModel, SoysterMaximinModel, DeJongMaximinModel, RobustMaximinModel, CCMaximinModel, MinimaxModel};
    plots_folder = "plots",
)
```