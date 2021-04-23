# ATA.jl: Automated Test Assembly Made Easy

Version 0.16.0

A package for automated test assembly (ATA) written in Julia.

ATA.jl is an open-source tool for building optimal test forms.
It is completely free (if you choose an open-source solver) and, 
on an average device (e.g. 
[Colab notebook](https://github.com/giadasp/ATA.jl/blob/master/examples/ATA.jl%20automated%20test%20assembly%20with%20Julia.ipynb)) is capable of producing a large number of test forms which can meet several
content and/or psychometric specifications.

For likely unsolvable extremely large scale ATA models, a solver based on the
simulated annealing heuristics is available (siman solver).
Otherwise, any mixed-integer linear programming (MILP) solver compatible with `JuMP.jl` can be used.

Look at the paper [^Spaccapanico2020] for some examples of ATA constraints and 
feasibility and size issues that can be encountered while optimizing ATA models.

Interfaced by `Dash` (via [`ATADash.jl`](https://github.com/giadasp/ATADash.jl)),
or running pure Julia code ([examples folder in this repo](https://github.com/giadasp/ATA.jl/blob/master/examples))
locally or on a hosted machine ([Colab notebook](https://github.com/giadasp/ATA.jl/blob/master/examples/ATA.jl%20automated%20test%20assembly%20with%20Julia.ipynb)).

## Documentation

Documentation on exported functions available at [link](https://giadasp.github.io/ATA.jl/docs) (in progress).

## Objectives

1. no objective;
2. MAXIMIN TIF (minimum between multiple ability points supported);
3. Chance-constrained MAXIMIN TIF [^Spaccapanico2021];
4. 3-standard deviations corrected MAXIMIN TIF [^Soyster1973];
5. 1-standard deviations corrected MAXIMIN TIF [^DeJong2009];
6. robust MAXIMIN TIF [^Veldkamp2013];
7. MINIMAX TIF (maximum between multiple ability points supported);
8. custom objective function.

## Constraints

- parallel tests (one group of tests);
- non-parallel tests (more than one group of tests);
- minimum and/or maximum length;
- minimum and/or maximum number of items with certain features (categorical constraints);
- maximum and/or minimum sum of numerical variables (quantitative constraints);
- maximum and/or minimum expected scores at multiple ability point; 
- maximum and/or minimum expected score based on the IRT paradigm or by a given column in the pool dataframe;
- minimum and/or maximum item use;
- maximum and/or minimum mean of numerical variables (quantitative constraints); (will be implemented soon, in the meanwhile look at the Colab notebook for a workaround)
- maximum overlap between tests; (increases dramatically the size and complexity of MILP ATA models `:boom:`, we suggest to avoid this constraint if the ATA model is already very large or, alternatively, to use the siman solver `:smirk_cat:`)
- group by units (friend sets);
- items exclusivity (enemy sets);

## Solvers

- for objectives 1, 2, 4, 5, 6, and 7: JuMP MILP solvers (see the list at https://jump.dev/JuMP.jl/v0.21.1/installation/). Follow the installation manuals of the solver you want to use. We suggest CBC or GLPK (default) as open-source solvers. Also commercial solvers are supported, CPLEX and Gurobi are examples.
- for objectives 1, 2, 3, 4, 5, 7, and 8 : Simulated Annealing solver (siman), pure Julia ATA solver.

## Report

Summarizing features of the assembled tests and plots are available.
Add [`ATAPlot`](https://github.com/giadasp/ATAPlot.jl) to plot the Test Information Functions (TIFs) and Item Characteristic Functions (ICTs).

## How to

If you want to play with this package:

Install Julia-1.6.0 at `https://julialang.org/downloads/`

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
It requires the installation of the package [`ATADash.jl`](https://github.com/giadasp/ATADash.jl)

As a remind, to plot the TIFs and ICTs you must install (add) and load (using) the package [`ATAPlot.jl`](https://github.com/giadasp/ATAPlot.jl).
It requires to have `PGFPlotsX.jl` and Miktex installed. 
A new version of `ATAPlot.jl` which use the package `Plots` and the backend `GR` is under development. 

## Documentation Contents

```@contents
Pages = ["build.md", "opt.md", "compact.md", "print_and_plot.md", "structs.md", "utils.md", "examples.md"]
Depth = 3
```

## Bug reports

Please report any issues via the Github [issue tracker]. All types of issues 
are welcome and encouraged; this includes bug reports, documentation typos,
feature requests, etc.

[issue tracker]: https://github.com/giadasp/ATA.jl/issues

## Citing ATA.jl

If you find `ATA.jl` useful in your work, we kindly request that you cite the
following github repository (it hasn't a DOI yet):

```bibtex
@misc{ATAjl,
    author       = {Spaccapanico Proietti, Giada},
    title        = {{ATA.jl: Automated Test Assembly Made Easy}},
    month        = apr,
    year         = 2021,
    version      = {0.16.0},
    url          = {https://github.com/giadasp/ATA.jl}
    }
```

[^Spaccapanico2020]: Spaccapanico Proietti, G., Matteucci, M., Mignani, S. (2020). Automated Test Assembly for Large-Scale Standardized Assessments: Practical Issues and Possible Solutions. Psych.; 2(4):315-337. https://doi.org/10.3390/psych2040024

[^Spaccapanico2021]: Spaccapanico Proietti, G., Matteucci, M., Mignani, S., Veldkamp, B. P. (2020). Chance-Constrained Automated Test Assembly. http://amsacta.unibo.it/6401/

[^Soyster1973]: Soyster, A.L. (1973) Technical Note—Convex Programming with Set-Inclusive Constraints and Applications to Inexact Linear Programming. Operations Research 21 (5) 1154-1157 https://doi.org/10.1287/opre.21.5.1154

[^DeJong2009]: De Jong, M. G., Steenkamp, J.-B. E. M., Veldkamp, B. P. (2009). A model for the construction of country-specific yet internationally comparable short-form marketing scales. Marketing Science, 28, 674-689. doi:10.1287/mksc.1080.0439