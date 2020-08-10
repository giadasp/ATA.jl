cd("where your input files are")
using ATAjl

#For using the Dash APP:

RunATA!()
#and navigate with the browser to localhost:8080
# Before running the app, if you want to use a MILP solver, remember to load it
# (ex: using Cbc; RunATA!()).

#If you prefere to use Julia code:

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
println(AddConstr!(ATAmodel; constraints_file = "Constraints.csv", constraints_delim=";")[2])

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
n_item_sample = maximum([ATAmodel.Cosntraints[t].length_max for t=1:ATAmodel.Settings.T])
# 1 <= n_item_sample <= Inf
n_test_sample = ATAmodel.Settings.T
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

#MILP (Not suggested for large scale ATA)
#Select the solver, Cbc as open-source is a good option.
#add Cbc by running import Pkg; Pkg.add("Cbc")
using Cbc
optimizer_constructor = Cbc.Optimizer
#Optimizer attributes
optimizer_attributes = [("seconds",100), ("logLevel",1)]

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
