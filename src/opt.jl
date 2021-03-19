include("simAn.jl")

"""
    assemble!(ata_model::AbstractModel; solver="jumpATA", starting_design=Matrix{Float64}(undef, 0, 0), results_folder="RESULTS", start_temp=0.1, geom_temp=0.1, n_item_sample=1, n_test_sample=1, opt_feas=0.0, n_fill=1, max_time=1000.00, max_conv=2, feas_nh=0, opt_nh=5, verbosity=2, optimizer_constructor="GLPK", optimizer_attributes=[("tm_lim", 1000)])

# Description

Assemble the tests.

# Arguments

- **`ata_model::AbstractModel`** : Required. The model built with ATA fuctions. 
- **`solver`** : Optional. Default: `"jumpATA"`. Values: `"jumpATA"`, `"siman"`. The solving interface to be used (JuMP or internal solver based on Simulated Annealing).
- **`starting_design`** : Optional. Default: `Matrix{Float64}(undef, 0, 0)`. The starting design matrix. Must be a `Matrix{Float64}`.
- **`results_folder`** : Optional. Default: `"RESULTS"`. The folder in which the output is stored.

## siman arguments

  - **`start_temp`** : Optional. Default: `0.1`. Values: `[0, Inf]`. Starting temperature, set to minimum for short journeys (if 0 worse solutions will never be accepted).
  - **`geom_temp`** : Optional. Default: `0.1`. Values: `[0, Inf)`. Decreasing geometric factor.
  - **`n_item_sample`** : Optional. Default: `1`. Values: `[1, Inf]`. Number of items to alter. Set to minimum for a shallow analysis, set to maximum for a deep analysis of the neighbourhoods.
  - **`n_test_sample`** : Optional. Default: `1`. Values: `[1, Inf]`. Number of tests to alter. Set to minimum for a shallow analysis, set to maximum for a deep analysis of the neighbourhoods.
  - **`opt_feas`** : Optional. Default: `0.0`. Values: `[0, Inf)`. Optimality/feasibility balancer, if 0 only feasibility of solution is analysed. Viceversa, if 1, only optimality is considered (uncontrained model). All the other values are accepted but produce uninterpretable results.
  - **`n_fill`** : Optional. Default: `1`. Values: `[0, Inf)`. Number of fill-up phases, usually 1 is sufficient, if start_temp is high it can be higher. If a starting_design is supplied, it can be set to 0.
  - **`verbosity`** : Optional. Default: `2`. Values: `1` (minimal), `2` (detailed). Verbosity level. In the console '+' stands for improvement, '_' for accepting worse solution. The dots are the fill-up improvement steps.
    
    * Termination criteria

      - **`max_time`** : Optional. Default: `1000.0`. Values: `[0, Inf)`. Time limit in seconds.
      - **`max_conv`** : Optional. Default: `2`. Values: `[1, Inf)`. Maximum convergence, stop when, after max_conv rounds no improvements have been found. Set to minimum for shallow analysis, increase it for deep analysis of neighbourhoods.
      - **`feas_nh`** : Optional. Default: `0`. Values: `[1, Inf)`. Maximum number of Feasibility neighbourhoods to explore, set to the minimum if the model is small or not highly constrained.
      - **`opt_nh`** : Optional. Default: `5`. Values: `[1, Inf)`. Maximum number of Optimality neighbourhoods to explore, set to the minimum if the model is highly constrained.

## jumpATA arguments

  - **`optimizer_constructor`** : Optional. Default: `"GLPK"`. Values: `"GLPK"`, `"Knitro"`, `"Gurobi"`, `"Cbc"`, `"CPLEX"`, `"Xpress"`, `"SCIP"`, `"Juniper"`. JuMP solver selection. Remember to load the required package before assemble!.
  - **`optimizer_attributes`** : Optional. Default: `[("tm_lim", 1000)]`. Values: An array of pairs `(attribute, value)`. Attributes and related values for the JuMP solver.

## other keyword arguments
  - **`kwargs...`** : Optional. 
"""
function assemble!(
    ata_model::AbstractModel;
    solver = "jumpATA",
    starting_design = Matrix{Float64}(undef, 0, 0),
    results_folder = "RESULTS",
    start_temp = 0.1,
    geom_temp = 0.1,
    n_item_sample = 1,
    n_test_sample = 1,
    opt_feas = 0.0,
    n_fill = 1,
    max_time = 1000.0,
    max_conv = 2,
    feas_nh = 0,
    opt_nh = 5,
    verbosity = 2,
    optimizer_constructor = "GLPK",
    optimizer_attributes = [("tm_lim", 1000)],
    kwargs...
)
    if solver == "siman"
        siman!(
            ata_model;
            results_folder = results_folder,
            starting_design = starting_design,
            start_temp = start_temp,
            geom_temp = geom_temp,
            max_time = max_time,
            n_item_sample = n_item_sample,
            n_test_sample = n_test_sample,
            max_conv = max_conv,
            verbosity = verbosity,
            opt_feas = opt_feas,
            n_fill = n_fill,
            feas_nh = feas_nh,
            opt_nh = opt_nh,
            kwargs...
        )
    elseif solver == "jumpATA"
        jumpATA!(
            ata_model;
            results_folder = results_folder,
            starting_design = starting_design,
            optimizer_constructor = optimizer_constructor,
            optimizer_attributes = optimizer_attributes,
            kwargs...
        )
    else
      push!(ata_model.output.infos, ["danger", "only \"siman\" and \"jumpATA\" are supported.\n"])
    end
    return nothing
end
