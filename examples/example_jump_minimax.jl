using ATA
using JuMP
# add Cbc by running import Pkg; Pkg.add("Cbc")
using Cbc
# cd("where your input files are")

# 1. Start ATA and add file with custom settings (Needed)
ata_model = start_ata(
    settings_file = "settings_ata minimax.jl",
    bank_file = "data/bank.csv",
    bank_delim = ";",
)
print_last_info(ata_model)

# 2. Add categorical constraints (Optional)
add_constraints!(ata_model; constraints_file = "constraints.csv", constraints_delim = ";")
print_last_info(ata_model)

# 3. Add overlap maxima (Optional, Needed if friend sets needs to be taken into account.)
group_by_friends!(ata_model)
print_last_info(ata_model)

# 4. Add objective function (Optional)
add_obj_fun!(ata_model)
print_last_info(ata_model)

# Assembly settings

# Set the solver, "siman" for simulated annealing, "jumpATA" for MILP solver.
solver = "jumpATA"

# MILP (Not suggested for large scale ATA)
# Select the solver, Cbc as open-source is a good option.
optimizer_constructor = "Cbc"
# #Optimizer attributes
optimizer_attributes = [("seconds", 500), ("logLevel", 1)]

# 5. assemble
assemble!(
    ata_model;
    solver = solver,
    optimizer_attributes = optimizer_attributes,
    optimizer_constructor = optimizer_constructor,
)

# All the settings and outputs from optimization are in ata_model object.
# See the struct in ATA.jl to understand how to retrieve all the information.
# A summary of the resulting tests is available in results_folder/results.txt
# If siman is chosen, the optimality and feasibility of the best neighbourhood
# is reported in "results/results_ata.jl"

print_results(ata_model; group_by_fs = true, results_folder = "results")

#]add https://github.com/giadasp/ATAPlot.jl
using ATAPlot
plot_results(ata_model; group_by_fs = true, results_folder = "plots")
