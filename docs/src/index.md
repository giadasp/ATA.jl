```@meta
CurrentModule = ATA
DocTestSetup = quote
    using ATA
end
```
# ATA.jl: Automated Test Assembly Made Easy

Version 0.10.0

A package for automated test assembly (ATA) written in Julia.
Simulated Annealing algorithm is available for large scale ATA models.
Otherwise, any MILP solver compatible with JuMP can be used.
Interfaced with Dash or pure Julia code.

## Documentation

Documentation on exported functions available at [link](https://giadasp.github.io/ATA.jl/docs) (in progress).

## Objectives

1. no objective;
2. MAXIMIN TIF (minimum between multiple ability points supported);
3. Chance-constrained MAXIMIN TIF;
4. MINIMAX TIF (maximum between multiple ability points supported);
5. custom objective function (only for siman solver)

## Constraints

- parallel tests (one group of tests);
- non-parallel tests (more than one group of tests);
- minimum and/or maximum length;
- minimum and/or maximum number of items with certain features (categorical constraints);
- maximum and/or minimum sum of numerical variables (quantitative constraints);
- maximum and/or minimum expected scores at multiple ability point; 
- maximum and/or minimum expected score based on the IRT paradigm or by a given column in the pool dataframe;
- minimum and/or maximum item use;
- maximum and/or minimum mean of numerical variables (quantitative constraints); (will be implemented soon)
- maximum overlap between tests; (increases dramatically the size and complexity of MILP ATA models `:boom:`, we suggest to avoid this constraint if the ATA model is already very large or, alternatively, to use the siman solver `:smirk_cat:`)
- group by units (friend sets);
- items exclusivity (enemy sets);

## Solvers

- for objectives 1, 2, and 4: JuMP MILP solvers (see the list at https://jump.dev/JuMP.jl/v0.21.1/installation/). Follow the installation manuals of the solver you want to use. We suggest CBC or GLPK (default) as open-source solvers. Also commercial solvers are supported, CPLEX and Gurobi are examples.
- for objectives 1, 2, 3, 4 and 5 (all): Simulated Annealing solver (siman), pure Julia ATA solver.

## Report

Summarizing features of the assembled tests and plots are available.
Add [`ATAPlot`](https://github.com/giadasp/ATAPlot.jl) to plot the Test Information Functions (TIFs) and Item Characteristic Functions (ICTs).

## How to

If you want to play with this package:

Install Julia-1.6.0-rc1 at `https://julialang.org/downloads/`

run `Julia.exe` and type:

```
] add https://github.com/giadasp/ATA.jl
```

Load the package by

```
using ATA
```

Play with the test files in folder "examples".

For a quick and compact ATA instance, run the code in "example_compact.jl".

If you want to dig in all the ATA building, assembly and output steps, run the code in the other files in the folder.
The comments explain which arguments must be customized in order to solve your ATA model. 
If you do not modify the arguments you just solve a toy ATA model with 3 non parallel groups of tests assembled starting from a 366 items bank.

For an even more easier ATA experience, look at "example_dash_app.jl".
It requires the installation of the package [`ATADash`](https://github.com/giadasp/ATADash.jl)

As a remind, to plot the TIFs and ICTs you must install (add) and load (using) the package [`ATAPlot`](https://github.com/giadasp/ATAPlot.jl).


## Documentation Contents

```@contents
Pages = ["build.md", "opt.md", "compact.md", "print_and_plot.md", "structs.md", "utils.md", "examples.md"]
Depth = 3
```

