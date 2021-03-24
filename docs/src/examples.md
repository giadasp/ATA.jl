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
ata_model = start_ata(settings_file="examples/settings_ata_maximin.jl", bank_file="examples/data/bank.csv", bank_delim=";")
print_last_info(ata_model)
```

2\. Add friend set variables (Optional)

```julia
add_friends!(ata_model)
print_last_info(ata_model)
```

3\. Add enemy set variables (Optional)

```julia
add_enemies!(ata_model)
print_last_info(ata_model)
```

4\. Add constraints (Optional)

```julia
add_constraints!(
    ata_model;
    constraints_file = "examples/constraints.csv",
    constraints_delim = ";",
)
print_last_info(ata_model)
```

5\. Add overlap maxima (Optional)

```julia
add_overlap!(ata_model; overlap_file = "examples/overlap_matrix.csv", overlap_delim = ";")
print_last_info(ata_model)
```

6\. Add expected score constraints (Optional)

```julia
add_exp_score!(ata_model)
print_last_info(ata_model)
```

7\. Group items by friend sets (Optional, Needed if `add_friends!(model)` has been run)

```julia
group_by_friends!(ata_model)
print_last_info(ata_model)
```

8\. Add objective function (Optional)

```julia
add_obj_fun!(ata_model)
print_last_info(ata_model)
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
  assemble!(ata_model;
      solver=solver,
      optimizer_attributes = optimizer_attributes,
      optimizer_constructor = optimizer_constructor
  )
```

10\. Print and Plot

All the settings and outputs are in the ata_model object.
See the struct in ATA.jl to understand how to retrieve all the information.
A summary of the resulting tests is available in <results_folder>/results.txt
If siman is chosen, the optimality and feasibility of the best neighbourhood
is reported in <results_folder>/results_ata.jl

```julia
print_results(ata_model; group_by_fs = true, results_folder = "results")
```

```julia
using ATAPlot
plot_results(ata_model; group_by_fs = true, results_folder = "plots")
```
