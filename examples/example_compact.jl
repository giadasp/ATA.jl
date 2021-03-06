using Distributed  #this is not needed if Julia has been run with <numberOfCores> >1
@everywhere using ATA
@everywhere cd("examples")
using ATAPlot
ata_model = compact_ata(;
    settings_file = "settings_ata maximin.jl",
    bank_file = "data/bank.csv",
    constraints_file = "constraints.csv",
    overlap_file = "overlap_matrix.csv",
    solver = "siman",
    max_time = 200,
)
