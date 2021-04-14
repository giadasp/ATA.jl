

"""
    add_bank!(
        ata_model::AbstractModel;
        settings::InputSettings = InputSettings(),
        bank::DataFrames.DataFrame = DataFrames.DataFrame(),
        settings_file = "settings_ata.jl",
        bank_file = "bank.csv",
        bank_delim = ";",
    )

# Description

Attach an item bank to an initialized ATA model.
The name of a `bank_delim` separated value file can be passed with the argument `bank_file`.
Alternatively, a bank dataframe can be passed using the argument `bank`.

# Arguments

- **`ata_model::AbstractModel`** : Required.
- **`bank::DataFrames.DataFrame`** : Optional. Default: `DataFrame()`. A dataframe containing data about the items.
- **`bank_file`** : Optional. Default: "bank.csv". The path of the file containing the item pool/bank in the form of custom-separated values.
- **`bank_delim`** : Optional. Default: ";". The custom-separator for the bank_file.

# Output

- An ATA model.
"""
function add_bank!(
    ata_model::AbstractModel;
    settings::InputSettings = InputSettings(),
    bank::DataFrames.DataFrame = DataFrames.DataFrame(),
    settings_file = "settings_ata.jl",
    bank_file = "bank.csv",
    bank_delim = ";",
)
    try
        #load bank
        bank_loaded = DataFrame()
        if size(bank, 1) > 0
            bank_loaded = bank
        elseif isfile(bank_file)
            try
                bank_loaded = CSV.read(bank_file, DataFrames.DataFrame, delim = bank_delim)
            catch e
                error!(
                    ata_model,
                    string(
                        "Error in reading the item bank file:\n  ",
                        sprint(showerror, e),
                        ".",
                    ),
                )
            end
        else
            error!(
                ata_model,
                string(
                    "Item bank file with name ",
                    bank_file,
                    " does not exist.\n  Provide a valid item bank dataframe or a name of an existing file.",
                ),
            )
        end
        if size(bank_loaded, 1) == ata_model.settings.n_items
            success!(ata_model, "Item bank file loaded.")
            ata_model.settings.bank = bank_loaded
            _add_irt!(ata_model)
            if ata_model.settings.fs.var != Symbol[]
                _add_friends!(ata_model)
            end
            if ata_model.settings.es.var != Symbol[]
                _add_enemies!(ata_model)
            end
            if any(
                vcat(
                    [
                        ata_model.constraints[t].expected_score.max for
                        t = 1:ata_model.settings.T
                    ]...,
                ) .< 1.00,
            ) || any(
                vcat(
                    [
                        ata_model.constraints[t].expected_score.min for
                        t = 1:ata_model.settings.T
                    ]...,
                ) .> 0.0,
            )
                _add_exp_score!(ata_model)
            end
        else
            error!(
                ata_model,
                string(
                    "Error in loading the item bank dataframe:\n number of rows is not equal to n_items.",
                ),
            )
        end
    catch e
        error!(ata_model, string(sprint(showerror, e)))
    end
    
    return nothing
end
