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

10\. Assemble

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

## Example with siman solver (Simulated Annealing Heuristic)

Suitable for large-scale ATA instances and custom objective functions optimization.

0\. Add and load required packages.
```julia
#cd("folder in which the package is saved")
using Pkg
Pkg.activate(".")  # required
Pkg.instantiate()
using ATA
cd("examples")
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

10\. Assemble.

Set the solver, "siman" for simulated annealing, "jumpATA" for MILP solver.
```julia
solver = "siman"
```

SIMAN (Suggested for Large scale ATA):

```julia
start_temp = 0.0001
```
Default: `0.1`. Values:  `[0, Inf]`. 
Starting temperature, set to minimum for short journeys (if 0 worse solutions will never be accepted).

```julia
geom_temp = 0.1
```
Default: `0.1`. Values:  `[0, Inf)`.
Decreasing geometric factor.

```julia
n_item_sample = Inf
```
Default: `1`. Values: `[1, Inf]`. 
Number of items to alter. Set to minimum for a shallow analysis, set to maximum for a deep analysis of the neighbourhoods.

```julia
n_test_sample = ATAmodel.settings.T
```
Default: `1`. Values: `[1, Inf]`. 
Number of tests to alter. Set to minimum for a shallow analysis, set to maximum for a deep analysis of the neighbourhoods.

```julia
opt_feas = 0.9
```
Default: `0.0`. Values: `[0.0, Inf)`. 
Optimality/feasibility balancer, if `0.0` only feasibility of solution is analysed. Viceversa, if `1.0`, only optimality is considered (uncontrained model). All the other values are accepted but produce uninterpretable results.
`
```julia
n_fill = 1
```
Default: `1`. Values: `[0, Inf)`.
Number of fill-up phases, usually 1 is sufficient, if start_temp is high it can be higher. 
If a `starting_design` is supplied, it should be set to `0`.

```julia
verbosity = 2
```
Default: `2`. Values: `1` (minimal), `2` (detailed).
Verbosity level. In the console '+' stands for improvement, '_' for accepting worse solution.
The dots are the fill-up improvement steps.

Termination criteria: 

```julia
max_time = 10.0
```
Default: `1000.0`. Values: `[0, Inf)`.
Time limit in seconds.

```julia
max_conv = 5
```
Default: `2`. Values: `[1, Inf)`. 
Maximum convergence, stop when, after max_conv rounds no improvements have been found. 
Set to minimum for shallow analysis, increase it for deep analysis of neighbourhoods.

```julia
feas_nh = 1
```
Default: `0`. Values: `[1, Inf)`. 
Maximum number of Feasibility neighbourhoods to explore, set to the minimum if the model is small or not highly constrained.

```julia
opt_nh = Inf
```
Default: `5`. Values: `[1, Inf)`. 
Maximum number of Optimality neighbourhoods to explore, set to the minimum if the model is highly constrained.


```julia
assemble!(ATAmodel;
    solver=solver,
    max_time=max_time,
    start_temp=start_temp,
    geom_temp=geom_temp,
    results_folder = results_folder,
    n_item_sample=n_item_sample,
    n_test_sample=n_test_sample,
    verbosity=verbosity,
    max_conv=max_conv,
    opt_feas=opt_feas,
    n_fill=n_fill,
    feas_nh=feas_nh,
    opt_nh=opt_nh
    )
```

For printing resulting tests and plots see step 11. in [first example](#Example-with-JuMP-(0.21.3)-and-Cbc)


## Custom objective function (only siman)

To add a custom objective function modify step 8 in previous examples.
It works only with the siman solver.

8\. Add objective function (Optional)

```julia
@info add_obj_fun!(ATAmodel)[2] 

ATAmodel.obj.type = "custom"

ATAmodel.obj.fun = function (x::Matrix{Float64}, obj_args::NamedTuple)
	IIF = obj_args.IIF
	T = size(IIF, 1)
	TIF = zeros(Float64, T)
	for t = 1:T
        K, I = size(IIF[t])
        # ungroup items
        xₜ = FS_to_items(x[:, t], obj_args.FS_items)
		if K > 1
			TIF[t] = Inf
			for k = 1:K
				TIF[t] = min(TIF, LinearAlgebra.dot(IIF[1, :], xₜ))
			end
		else
			TIF[t] = LinearAlgebra.BLAS.gemv('N', IIF, xₜ)[1]
		end
    end
    min_TIF = minimum(TIF)
    TIF = [min_TIF for t = 1:T] 
    # Must return a vector of length T.
    # The resulting objective function is the minimum of all values in this vector.
	return TIF::Vector{Float64}
end

ATAmodel.obj.args = (IIF = FileIO.load("data/IIF.jld2", "IIF"), FS_items = ATAmodel.settings.FS.items)
```

## Parallel neighbouhoods evaluation (only siman)

By using multiple cores, the Simulated Annealing solver evaluates several neighbourhoods at a time reducing the time needed for the exploration of the solution space.  
`Julia` starts with 1 core by default, to run `Julia` with `<ncores>`, type `julia -p <ncores>` in the command prompt and follow the [example](#Example-with-siman-solver-(Simulated-Annealing-Heuristic)) as usual.

