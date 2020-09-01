################################################################################
###                    Parallel siman algorithm                           ######
################################################################################

#1. First, save the package folder in you pc.
#2. Then, modify all the input files (included this) in the test folder as required.
#3. When the ATA model is ready run Julia.exe and write:
#cd("where your input files are")
# include("LocalTest.jl")

#[Beta: for running ATA with multiple cores do:
#add julia.exe to the path and run Julia with the command prompt:
#julia -p numberOfCores]

using Distributed #this is not needed if Julia has been run with numberOfCores>1
@everywhere   using DataFrames.DataFrames
@everywhere   using Distributed
@everywhere   cd("folder in which the package is saved")
@everywhere   using Pkg
@everywhere   Pkg.activate(".")  # required
@everywhere   Pkg.instantiate()
@everywhere   using ATA

#@everywhere cd("where your input files are")
#Input files needed: settingsATA.jl, bank.csv, CategoricalConstraints.csv, OverlapMatrix.csv, BSpar.jld2 (only for Chance-Contrained (CC))
#settingsATA.jl : overall features of the tests, such as length, expected score, item use, ecc...
#CategoricalConstraints.csv : categorical constraints
#OverlapMatrix.csv : maximum overlap allowed between test forms, it is a TxT matrix. If the assembly is non parallel (i.e. there is more than one group of parallel tests, n_groups>1)
#it is a n_groups x n_groups matrix.
#BSpar.jld2 : it is a nPar way array (Array{Float64,nPar}) where nPar is the number of IRT parameters. Each sub array is a n_items x R matrix (Matrix{Float64}(.,n_items,R)).

#for resetting the ATA process (Needed)
ATAmodel = start_ATA()

#Each of the following commands returns a string vector, the second element is a message describing the result.
#1. Add file with custom settings (Needed)
println(load_settings!(ATAmodel; settings_file = "settingsATA.jl", bank_file = "bank.csv", bank_delim = ";")[2])

#2. Add friend set variables (Optional)
println(add_friends!(ATAmodel)[2])

#3. Add enemy set variables (Optional)
println(add_enemies!(ATAmodel)[2])

#4. Add categorical constraints (Optional)
println(add_constraints!(ATAmodel; constraints_file = "constraints.csv", constraints_delim=";")[2])

#5. Add overlap maxima (Optional)
println(add_overlap!(ATAmodel; overlap_file = "Overlap Matrix.csv", overlap_delim=";")[2])

#6. Add expected score constraints (Optional)
println(add_exp_score!(ATAmodel)[2])

#7. Add overlap maxima (Optional, Needed if add_friends!(model) hase been run)
println(group_by_friends!(ATAmodel)[2])

#8. Add objective function (Optional)
println(add_obj_fun!(ATAmodel)[2])

#Assembly settings

# Set the solver, "siman" for simulated annealing, "jumpATA" for MILP solver.
solver = "siman"

# SIMAN (Suggested for Large scale ATA):

start_temp = 0.0001
# Default: `0.1`. Values:  `[0, Inf]`. 
# Starting temperature, set to minimum for short journeys (if 0 worse solutions will never be accepted).

geom_temp = 0.1
# Default: `0.1`. Values:  `[0, Inf)`.
# Decreasing geometric factor.

n_item_sample = Inf
# Default: 1. Values: `[1, Inf]`. 
# Number of items to alter. Set to minimum for a shallow analysis, 
# set to maximum for a deep analysis of the neighbourhoods.

n_test_sample = ATAmodel.settings.T
# Default: 1. Values: `[1, Inf]`. 
# Number of tests to alter. Set to minimum for a shallow analysis, set to maximum for a deep analysis of the neighbourhoods.

opt_feas = 0.9
# Default: 0.0. Values: `[0, Inf)`. 
# Optimality/feasibility balancer, if 0 only feasibility of solution is analysed. Viceversa, if 1, only optimality is considered (uncontrained model). All the other values are accepted but produce uninterpretable results.

n_fill = 1
# Default: 1. Values: `[0, Inf)`.
# Number of fill-up phases, usually 1 is sufficient, if start_temp is high it can be higher. 
# If a starting_design is supplied, it should be set to 0.
 
verbosity = 2
# Default: 2. Values: `1` (minimal), `2` (detailed).
# Verbosity level. In the console '+' stands for improvement, '_' for accepting worse solution.
# The dots are the fill-up improvement steps.

#! Termination criteria: 

max_time = 10.0
# Default: `1000.0`. Values: `[0, Inf)`.
# Time limit in seconds.

max_conv = 5
# Default: `2`. Values: `[1, Inf)`. 
# Maximum convergence, stop when, after max_conv rounds no improvements have been found. 
# Set to minimum for shallow analysis, increase it for deep analysis of neighbourhoods.

feas_nh = 1
# Default: `0`. Values: `[1, Inf)`. 
# Maximum number of Feasibility neighbourhoods to explore, set to the minimum if the model is small or not highly constrained.

opt_nh = Inf
# Default: `5`. Values: `[1, Inf)`. 
# Maximum number of Optimality neighbourhoods to explore, set to the minimum if the model is highly constrained.

#9. assemble
assemble!(ATAmodel;
    solver = solver,
    max_time = max_time,
    start_temp = start_temp,
    geom_temp = geom_temp,
    n_item_sample = n_item_sample,
    n_test_sample = n_test_sample,
    verbosity = verbosity,
    max_conv = max_conv,
    opt_feas = opt_feas,
    n_fill = n_fill,
    feas_nh = feas_nh,
    opt_nh = opt_nh,
    optimizer_attributes = optimizer_constructor,
    optimizer_constructor =optimizer_attributes
    )

# All the settings and outputs from optimization are in ATAmodel object.
# See the struct in ATA.jl to understand how to retrieve all the information.
# A summary of the resulting tests is available in results_folder/Results.txt
# If siman is chosen, the optimality and feasibility of the best neighbourhood
# is reported in "RESULTS/ResultsATA.jl"
print_results(ATAmodel;
group_by_fs = true,
plots_out = false,
results_folder = "RESULTS")


#to stop all the processes do:
#ctrl+C
#interrupt()
#addprocs(numberOfProcs)
#@everywhere cd("where your input files are")
#include("LocalTest.jl")
#Julia will compile the package again.
