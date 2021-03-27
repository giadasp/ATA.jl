# Build the ATA model

```@meta
CurrentModule = ATA
```

## Start ATA Model

```@docs
start_ata(;
    settings::InputSettings = InputSettings(),
    bank::DataFrames.DataFrame = DataFrames.DataFrame(),
    settings_file = "settings_ata.jl",
    bank_file = "bank.csv",
    bank_delim = ";",
)
```

## Add Constraints

```@docs
add_constraints!(
    ata_model::AbstractModel;
    constraints::DataFrames.DataFrame = DataFrames.DataFrame(),
    constraints_file = "constraints.csv",
    constraints_delim = ";",
)
```

## Add Overlap Matrix

```@docs
add_overlap!(
    ata_model::AbstractModel;
    overlap::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),
    overlap_file = "overlap_matrix.csv",
    overlap_delim = ";",
)
```

## Group by Friend Sets

```@docs
group_by_friends!(ata_model::AbstractModel)
```

## Add Objective Function

```@docs
add_obj_fun!(ata_model::Union{MaximinModel,MinimaxModel}; kwargs...)
add_obj_fun!(
    ata_model::SoysterMaximinModel;
    psychometrics = false,
    items::Vector{Psychometrics.Item} = Psychometrics.Item[],
    items_file = "",
    kwargs...
)
add_obj_fun!(
    ata_model::DeJongMaximinModel;
    psychometrics = false,
    items::Vector{Psychometrics.Item} = Psychometrics.Item[],
    items_file = "",
    kwargs...
)
add_obj_fun!(
    ata_model::CcMaximinModel;
    psychometrics = false,
    items_file = "",
    items::Vector{Psychometrics.Item} = Psychometrics.Item[],
    kwargs...
)
```

## Load Design

A design can be loaded as a starting solution for the ATA or to print (or plot) the results.

```@docs
load_design!(design::Matrix{Any}, ata_model::AbstractModel)
```


