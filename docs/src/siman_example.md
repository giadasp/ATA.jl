## Example with siman solver (Simulated Annealing heuristic)

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

### Assembly settings

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
Default: 2. Values: `1` (minimal), `2` (detailed).
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

10\. Assemble.

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
