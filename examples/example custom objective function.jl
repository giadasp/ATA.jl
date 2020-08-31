# cd("folder in which the package is saved")
# using Pkg
# Pkg.activate(".")  # required
# Pkg.instantiate()
# cd("where your input files are")
using ATA

# for custom obj_fun:
using CSV
using FileIO

# For using the Dash APP:
# run_app!()
# Navigate with the browser to localhost:8080
# Before running the app, if you want to use a MILP solver, remember to load it
# (ex: using Cbc; run_app!()).

# If you prefere to use Julia code:

# for resetting the ATA process (Needed)
ATAmodel = start_ATA()

# Each of the following commands returns a string vector, the second element is a message describing the result.
# 1. Add file with custom settings (Needed)
@info load_settings!(ATAmodel; settings_file="settingsATA.jl", bank_file="data/bank.csv", bank_delim=";")[2]

# 2. Add friend set variables (Optional)
@info add_friends!(ATAmodel)[2]

# 3. Add enemy set variables (Optional)
@info add_enemies!(ATAmodel)[2]

# 4. Add categorical constraints (Optional)
@info add_constraints!(ATAmodel; constraints_file="constraints.csv", constraints_delim=";")[2]

# 5. Add overlap maxima (Optional)
@info add_overlap!(ATAmodel; overlap_file="Overlap Matrix.csv", overlap_delim=";")[2]

# 6. Add expected score constraints (Optional)
@info add_exp_score!(ATAmodel)[2]

# 7. Add overlap maxima (Optional, Needed if add_friends!(model) hase been run)
@info group_by_friends!(ATAmodel)[2]

# 8. Add objective function (Optional)
@info add_obj_fun!(ATAmodel)[2] 

# custom objective type, function and arguments
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

# Assembly settings

# SIMAN (Suggested for Large scale ATA):
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


# 9. assemble
assemble!(ATAmodel;
    solver=solver,
    max_time=max_time,
    start_temp=start_temp,
    geom_temp=geom_temp,
    n_item_sample=n_item_sample,
    n_test_sample=n_test_sample,
    verbosity=verbosity,
    max_conv=max_conv,
    opt_feas=opt_feas,
    n_fill=n_fill,
    feas_nh=feas_nh,
    opt_nh=opt_nh# ,
    # optimizer_attributes = optimizer_constructor,
    # optimizer_constructor =optimizer_attributes
    )

# All the settings and outputs from optimization are in ATAmodel object.
# See the struct in ATA.jl to understand how to retrieve all the information.
# A summary of the resulting tests is available in results_folder/Results.txt
# If siman is chosen, the optimality and feasibility of the best neighbourhood
# is reported in "RESULTS/ResultsATA.jl"

using Plots
using PGFPlotsX

# requires lualatex installed 
ATAmodel.obj.type = "MAXIMIN"
print_results!(ATAmodel;
group_by_fs=true,
plots_out=false,
results_folder="RESULTS")
