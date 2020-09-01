```@meta
CurrentModule = ATA
DocTestSetup = quote
    using ATA
end
```
# ATA.jl: Automated Test Assembly Made Easy

A package for automated test assembly (ATA) written in Julia.
Simulated Annealing algorithm is available for large scale ATA models.
Otherwise, any MILP solver compatible with JuMP can be used.
Interfaced with Dash or pure Julia code.

## Objectives
1. no objective;
2. MAXIMIN TIF (minimum between multiple ability points supported);
3. Chance-constrained MAXIMIN TIF.


## Constraints
- parallel tests (one group of tests);
- non-parallel tests (more than one group of tests);
- minimum and/or maximum length;
- minimum and/or maximum number of items with certain features (categorical constraints);
- maximum and/or minimum sum of numerical variables (quantitative constraints);
- maximum and/or minimum expected scores at multiple ability point;
- minimum and/or maximum item use;
- maximum and/or minimum mean of numerical variables (quantitative constraints); (beta)
- maximum and/or minimum expected score based on the IRT paradigm or by a given column in the pool dataframe;
- maximum overlap between tests;
- group by units (friend sets);
- items exclusivity (enemy sets);

## Solvers
- for objectives 1 and 2: JuMP MILP solvers (see the list at https://jump.dev/JuMP.jl/v0.21.1/installation/). Follow the installation manuals of the solver you want to use. For open-source we suggest CBC or GLPK (default). For commercial we suggest Gurobi or CPLEX.
- for objectives 1, 2 and 3: Simulated Annealing solver, pure Julia ATA solver.

## Report
Summarizing features of the assembled tests and Plots of the ICFs and TIFs are available.

## How to

If you want to play with this package:

Install Julia-1.5.1 at `https://julialang.org/downloads/`

run `Julia.exe` and type:

```
] add https://github.com/giadasp/ATA.jl
```

Load the package by

```
using ATA
```

Play with the test files in folder "examples".

Set all the specifications by modifying the files in "examples" folder and run "example.jl" or "example with custom objective function.jl" following the manual in the comments.

Distributed analysis of the neighbourhoods is available, look at "examples/parallelExample.jl".


