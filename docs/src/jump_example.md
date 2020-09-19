

## Example with JuMP (0.21.3) and Cbc

0\. Add and load required packages.
```julia
#cd("folder in which the package is saved")
using Pkg
Pkg.activate(".")  # required
Pkg.instantiate()
using ATA
cd("examples")
Pkg.add("JuMP@0.21.3")
using JuMP
Pkg.add("Cbc")
using Cbc
```

1\. Resetting the ATA process (Needed)
  
```julia
ATAmodel = start_ATA()
```

2\. Add file with custom settings (Needed)

```julia
@info load_settings!(ATAmodel; settings_file="SettingsATA.jl", bank_file="data/Bank.csv", bank_delim=";")[2]
```

3\. Add friend set variables (Optional)

```julia
@info add_friends!(ATAmodel)[2]
```

4\. Add enemy set variables (Optional)

```julia
@info add_enemies!(ATAmodel)[2]
```

5\. Add categorical constraints (Optional)

```julia
@info add_constraints!(ATAmodel; constraints_file="Constraints.csv", constraints_delim=";")[2]
```

6\. Add overlap maxima (Optional)

```julia
@info add_overlap!(ATAmodel; overlap_file="OverlapMatrix.csv", overlap_delim=";")[2]
```

7\. Add expected score constraints (Optional)

```julia
@info add_exp_score!(ATAmodel)[2]
```

8\. Group items by friend sets (Optional, Needed if `add_friends!(model)` has been run)

```julia
@info group_by_friends!(ATAmodel)[2]
```

9\. Add objective function (Optional)

```julia
@info add_obj_fun!(ATAmodel)[2] 
```

### Assembly settings

Set the solver, "siman" for simulated annealing, "jumpATA" for MILP solver.
```julia
solver = "jumpATA"
```

MILP (Not suggested for large scale ATA)
Select the solver, Cbc as open-source is a good option:
```julia
optimizer_constructor = "Cbc"
```

Optimizer attributes:
```julia
optimizer_attributes = [("seconds", 100), ("logLevel", 1)]
```

10\. Assemble

```julia
  assemble!(ATAmodel;
      solver=solver,
      optimizer_attributes=optimizer_constructor,
      optimizer_constructor=optimizer_attributes
  )
```

11\. Print and Plot

All the settings and outputs from optimization are in ATAmodel object.
See the struct in ATA.jl to understand how to retrieve all the information.
A summary of the resulting tests is available in <results_folder>/Results.txt.
To print the plots the packages Plots and PGFPlotsX are required together with an installed implementation of LaTeX such as MikTein.

```julia

Pkg.add("Plots")
using Plots
Pkg.add("PGFPlotsX")
using PGFPlotsX

print_results!(ATAmodel;
group_by_fs=true,
plots_out=true,
results_folder="RESULTS")
```
