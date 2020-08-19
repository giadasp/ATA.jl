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
ATAmodel = StartATA()

#Each of the following commands returns a string vector, the second element is a message describing the result.
#1. Add file with custom settings (Needed)
println(LoadSettings!(ATAmodel; settings_file = "settingsATA.jl", bank_file = "bank.csv", bank_delim = ";")[2])

#2. Add friend set variables (Optional)
println(AddFriendSets!(ATAmodel)[2])

#3. Add enemy set variables (Optional)
println(AddEnemySets!(ATAmodel)[2])

#4. Add categorical constraints (Optional)
println(AddConstr!(ATAmodel; constraints_file = "constraints.csv", constraints_delim=";")[2])

#5. Add overlap maxima (Optional)
println(AddOverlaps!(ATAmodel; overlap_file = "Overlap Matrix.csv", overlap_delim=";")[2])

#6. Add expected score constraints (Optional)
println(AddExpScore!(ATAmodel)[2])

#7. Add overlap maxima (Optional, Needed if AddFriendSets!(model) hase been run)
println(GroupByFriendSet!(ATAmodel)[2])

#8. Add objective function (Optional)
println(AddObjFun!(ATAmodel)[2])

#Assembly settings

#Set the solver, "siman" for simulated annealing, "jumpATA" for MILP solver.
solver = "siman"

#SIMAN (Suggested for Large scale ATA):

#Stopping criteria:
#Time limit in seconds
max_time = 100
#Minimum convergence, stop when, after max_conv rounds no improvements have been found.
#Set to minimum for shallow analysis, increase it for deep analysis.
conv_max = 5 # 2 <= conv_max <= Inf

#Starting temperature, set to minimum for short journeys (if 0 worse solutions will never be accepted)
start_temp = 0.0001 #0 <= start_temp <= Inf
#Deacresing geometric factor
geom_temp  = 0.1  #0 <= start_temp <= Inf

#Number of fill up phases, usually 1 is sufficient, if start_temp is high can be setted high.
# if a starting_design is supplied, can be setted to 0.
n_fill = 1
#Set deep analysis of neighbourhoods, set both to 1 for a shallow analysis.
n_item_sample = maximum([ATAmodel.Cosntraints[t].length_max for t=1:ATAmodel.settings.T])
# 1 <= n_item_sample <= Inf
n_test_sample = ATAmodel.settings.T
# 1 <= n_test_sample <= Inf

#Number of Feasibility neighbourhoods to explore, set to the minimum if the model is small or not highly constrained
feas_nh = 1 # 0 <= feas_nh <= Inf
#Number of Optimality neighbourhoods to explore, set to the minimum if the model is highly constrained
opt_nh = 200 # 0 <= feas_nh <= Inf
#Optimality/feasibility balancer, if 0 only feasibility of solution is analysed.
#Viceversa, if 1, only optimality is considered (uncontrained model)
opt_feas = 0.9 # 0 <= opt_feas <= 1

#Verbosity, in the console '+' stands for improvement, '_' for accepting worse solution.
verbosity = 2 # 1 = minimal, 2 = detailed.

#MILP (Look at test.jl for use MILP solvers)

#9. Assemble
Assemble!(ATAmodel;
    solver = solver,
    max_time = max_time,
    start_temp = start_temp,
    geom_temp = geom_temp,
    n_item_sample = n_item_sample,
    n_test_sample = n_test_sample,
    verbosity = verbosity,
    conv_max = conv_max,
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
PrintResults(ATAmodel;
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
