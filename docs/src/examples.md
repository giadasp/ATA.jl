# Examples

## Example with JuMP (0.21.3) and Cbc


```jldoctest
#cd("folder in which the package is saved")
# using Pkg
# Pkg.activate(".")  # required
# Pkg.instantiate()
# cd("where your input files are")
julia> using ATA
julia> using JuMP
# add Cbc by running import Pkg; Pkg.add("Cbc")
julia> using Cbc
```

1. Resetting the ATA process (Needed)

```jldoctest
julia> ATAmodel = start_ATA()
```

2. Add file with custom settings (Needed)
```jldoctest
julia> @info load_settings!(ATAmodel; settings_file="settingsATA.jl", bank_file="data/bank.csv", bank_delim=";")[2]
```

3. Add friend set variables (Optional)
```jldoctest
julia> @info add_friends!(ATAmodel)[2]
```

4. Add enemy set variables (Optional)
```jldoctest
julia> @info add_enemies!(ATAmodel)[2]
```

5. Add categorical constraints (Optional)
```jldoctest
julia> @info add_constraints!(ATAmodel; constraints_file="constraints.csv", constraints_delim=";")[2]
```

6. Add overlap maxima (Optional)
```jldoctest
julia> @info add_overlap!(ATAmodel; overlap_file="Overlap Matrix.csv", overlap_delim=";")[2]
```

7. Add expected score constraints (Optional)
```jldoctest
julia> @info add_exp_score!(ATAmodel)[2]
```

8. Add overlap maxima (Optional, Needed if add_friends!(model) hase been run)
```jldoctest
julia> @info group_by_friends!(ATAmodel)[2]
```

9. Add objective function (Optional)
```jldoctest
julia> @info add_obj_fun!(ATAmodel)[2] 
```

### Assembly settings

Set the solver, "siman" for simulated annealing, "jumpATA" for MILP solver.
```jldoctest
julia> solver = "jumpATA"
```

MILP (Not suggested for large scale ATA)
Select the solver, Cbc as open-source is a good option:
```jldoctest
julia> optimizer_constructor = "Cbc"
```

Optimizer attributes:
```jldoctest
julia> optimizer_attributes = [("seconds", 100), ("logLevel", 1)]
```

10. assemble
```jldoctest
julia> assemble!(ATAmodel;
    solver=solver,
    optimizer_attributes=optimizer_constructor,
    optimizer_constructor=optimizer_attributes
    )
```


All the settings and outputs from optimization are in ATAmodel object.
See the struct in ATA.jl to understand how to retrieve all the information.
A summary of the resulting tests is available in <results_folder>/Results.txt

```jldoctest
julia> print_results!(ATAmodel;
group_by_fs=true,
plots_out=true,
results_folder="RESULTS")
```

