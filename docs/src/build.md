# Build the ATA model

```@docs
start_ATA()
load_settings!(::ATA.Model; settings_file = "SettingsATA.jl", bank_file = "bank.csv", bank_delim = ";")
add_friends!(::ATA.Model)
add_enemies!(::ATA.Model)
add_constraints!(::ATA.Model; constraints_file = "Constraints.csv", constraints_delim = ";")
add_overlap!(::ATA.Model; overlap_file = "OverlapMatrix.csv", overlap_delim=";")
add_exp_score!(::ATA.Model)
add_obj_fun!(::ATA.Model)
group_by_friends!(::ATA.Model)
load_design!(::Matrix{Any}, ::ATA.Model)
```