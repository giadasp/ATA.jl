__precompile__()
module ATA

import Printf
import CSV
import DelimitedFiles
import DataFrames.DataFrames
import Distributions.Distributions
import Distributed
import Random
import LinearAlgebra
import JLD2
import StatsBase
using Requires

include("structs.jl")
include("utils.jl")
include("build.jl")
include("opt.jl")
include("print.jl")

function __init__()
    @require Cbc = "9961bab8-2fa3-5c5a-9d89-47fab24efd76" begin
        function add_Cbc!(m)
            JuMP.set_optimizer(m, Cbc.Optimizer)
        end
    end
    @require Gurobi = "2e9cd046-0924-5485-92f1-d5272153d98b" begin
        function add_Gurobi!(m)
            JuMP.set_optimizer(m, Gurobi.Optimizer)
        end
    end
    @require CPLEX = "a076750e-1247-5638-91d2-ce28b192dca0" begin
        function add_CPLEX!(m)
            JuMP.set_optimizer(m, CPLEX.Optimizer)
        end
    end
    @require GLPK = "60bf3e95-4087-53dc-ae20-288a0d20c6a6" begin
        function add_GLPK!(m)
            JuMP.set_optimizer(m, GLPK.Optimizer)
        end
    end
    @require Xpress = "9e70acf3-d6c9-5be6-b5bd-4e2c73e3e054" begin
        function add_Xpress!(m)
            JuMP.set_optimizer(m, Gurobi.Optimizer)
        end
    end
    @require KNITRO = "67920dd8-b58e-52a8-8622-53c4cffbe346" begin
        function add_KNITRO!(m)
            JuMP.set_optimizer(m, KNITRO.Optimizer)
        end
    end
    @require Juniper = "2ddba703-00a4-53a7-87a5-e8b9971dde84" begin
        function add_Juniper!(m)
            JuMP.set_optimizer(m, Juniper.Optimizer)
        end
    end
    @require SCIP = "bdc5fa1a-e2a0-4ec3-bc77-c59bd234be9e" begin
        function add_SCIP!(m)
            JuMP.set_optimizer(m, SCIP.Optimizer)
        end
    end
    @require ATAPlot = "372623d9-1dd6-4d20-a513-3f20705132c0" begin
        include("plot.jl")
    end
    @require ATADash = "236b7dbe-4167-40cb-a459-bcf8ce4b2cbd" begin
        function run_app() 
            ATADash.ATA_app() 
        end
    end
    @require JuMP = "4076af6c-e467-56ae-b986-b466b2749572" include("jump.jl")
end

export mean,
    start_ATA,
    load_settings!,
    add_friends!,
    add_enemies!,
    add_constraints!,
    add_overlap!,
    add_exp_score!,
    group_by_friends!,
    add_obj_fun!,
    assemble!,
    print_results,
    plot_results,
    load_design!,
    item_info,
    item_char,
    resp_gen,
    student_likelihood,
    FS_to_items,#_mycopy,  fast_sort!, _myqle#, optimize
    run_app
ATA

end
