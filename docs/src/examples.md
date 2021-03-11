# Examples

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

1\. Resetting the ATA process and load settings (Needed)

```julia
ATAmodel = start_ATA(settings_file="examples/settingsATA maximin.jl", bank_file="examples/data/Bank.csv", bank_delim=";")
print_last_info(ATAmodel)
```

2\. Add friend set variables (Optional)

```julia
add_friends!(ATAmodel)
print_last_info(ATAmodel)
```

3\. Add enemy set variables (Optional)

```julia
add_enemies!(ATAmodel)
print_last_info(ATAmodel)
```

4\. Add constraints (Optional)

```julia
add_constraints!(
    ATAmodel;
    constraints_file = "examples/Constraints.csv",
    constraints_delim = ";",
)
print_last_info(ATAmodel)
```

5\. Add overlap maxima (Optional)

```julia
add_overlap!(ATAmodel; overlap_file = "examples/OverlapMatrix.csv", overlap_delim = ";")
print_last_info(ATAmodel)
```

6\. Add expected score constraints (Optional)

```julia
add_exp_score!(ATAmodel)
print_last_info(ATAmodel)
```

7\. Group items by friend sets (Optional, Needed if `add_friends!(model)` has been run)

```julia
group_by_friends!(ATAmodel)
print_last_info(ATAmodel)
```

8\. Add objective function (Optional)

```julia
add_obj_fun!(ATAmodel)
print_last_info(ATAmodel)
```

9\. Assemble

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

```julia
  assemble!(ATAmodel;
      solver=solver,
      optimizer_attributes = optimizer_attributes,
      optimizer_constructor = optimizer_constructor
  )
```

10\. Print and Plot

All the settings and outputs are in the ATAmodel object.
See the struct in ATA.jl to understand how to retrieve all the information.
A summary of the resulting tests is available in <results_folder>/Results.txt
If siman is chosen, the optimality and feasibility of the best neighbourhood
is reported in <results_folder>/ResultsATA.jl

```julia
print_results(ATAmodel; group_by_fs = true, results_folder = "RESULTS")
```

```julia
using ATAPlot
plot_results(ATAmodel; group_by_fs = true, results_folder = "PLOTS")
```
