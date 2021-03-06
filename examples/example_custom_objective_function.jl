using ATA
# For load data to input in custom obj_fun:
using CSV
using FileIO
using LinearAlgebra
cd("where your input files are")

# 1. Start ATA and add file with custom settings (Needed)
ata_model = start_ata(;
    settings_file = "settings_ata_custom.jl",
    bank_file = "data/bank.csv",
    bank_delim = ";",
)

# 2. Add categorical constraints (Optional)
add_constraints!(ata_model; constraints_file = "constraints.csv", constraints_delim = ";")

# 3. Add overlap maxima (Optional)
add_overlap!(ata_model; overlap_file = "overlap_matrix.csv", overlap_delim = ";")

# custom objective type, function and arguments
ata_model.obj.name = "custom"

ata_model.obj.fun = function (x::Matrix{Float64}, obj_args::NamedTuple)
    IIF = obj_args.IIF
    T = obj_args.T
    TIF = zeros(Float64, T)
    for t = 1:T
        K, I = size(IIF[t])
        # ungroup items
        xₜ = fs_to_items(x[:, t], I, obj_args.fs_items)
        if K > 1
            TIF[t] = Inf
            for k = 1:K
                TIF[t] = min(TIF[t], LinearAlgebra.dot(IIF[t][k, :], xₜ))
            end
        else
            TIF[t] = LinearAlgebra.dot(IIF[t][1, :], xₜ)[1]
        end
    end
    min_TIF = minimum(TIF)
    TIF = [min_TIF for t = 1:T]
    # Must return a vector of length T.
    # The resulting objective function is the minimum of all values in this vector.
    return TIF::Vector{Float64}
end

ata_model.obj.args = (
    T = ata_model.settings.T,
    IIF = FileIO.load("data/IIF.jld2", "IIF"),
    fs_items = ata_model.settings.fs.items,
)

print_infos(ata_model)
# Assembly settings

# SIMAN (Suggested for Large scale ATA):
# Set the solver, "siman" for simulated annealing, "jump" for MILP solver.
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

n_test_sample = ata_model.settings.T
# Default: 1. Values: `[1, Inf]`. 
# Number of tests to alter. Set to minimum for a shallow analysis, set to maximum for a deep analysis of the neighbourhoods.

opt_feas = 0.9
# Default: 0.0. Values: `[0, Inf)`. 
# Optimality/feasibility balancer, if 0 only feasibility of solution is analysed. Viceversa, if 1, only optimality is considered (uncontrained model). All the other values are accepted but produce uninterpretable results.

n_fill = 1
# Default: 1. Values: `[0, Inf)`.
# Number of fill-up phases, usually 1 is sufficient, if start_temp is high it can be higher. 
# If a starting_design is supplied, it should be set to 0.

verbosity = 1
# Default: 2. Values: `1` (minimal), `2` (detailed).
# Verbosity level. In the console '+' stands for improvement, '_' for accepting worse solution.
# The dots are the fill-up improvement steps.

#! Termination criteria: 

max_time = 100.0
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


# 4. assemble
assemble!(
    ata_model;
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
    opt_nh = opt_nh,# ,
    # optimizer_attributes = optimizer_constructor,
    # optimizer_constructor =optimizer_attributes
)

# All the settings and outputs from optimization are in ata_model object.
# See the struct in ATA.jl to understand how to retrieve all the information.
# A summary of the resulting tests is available in results_folder/results.txt
# If siman is chosen, the optimality and feasibility of the best neighbourhood
# is reported in "results/results_ata.jl"

ata_model.obj.name = "maximin"
print_results(ata_model; results_folder = "results")


#]add https://github.com/giadasp/ATAPlot.jl
using ATAPlot
plot_results(ata_model; plots_folder = "plots")
