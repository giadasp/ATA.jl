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

