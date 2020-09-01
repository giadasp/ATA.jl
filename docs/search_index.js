var documenterSearchIndex = {"docs":
[{"location":"opt/#Assemble-tests","page":"Assemble tests","title":"Assemble tests","text":"","category":"section"},{"location":"opt/","page":"Assemble tests","title":"Assemble tests","text":"Modules = [ATA]\nPages   = [\"opt.jl\"]","category":"page"},{"location":"opt/#ATA.assemble!-Tuple{ATA.Model}","page":"Assemble tests","title":"ATA.assemble!","text":"assemble!(ATAmodel::Model; solver=\"jumpATA\", starting_design=Matrix{Float64}(undef, 0, 0), results_folder=\"RESULTS\", start_temp=0.1, geom_temp=0.1, n_item_sample=1, n_test_sample=1, opt_feas=0.0, n_fill=1, max_time=1000.00, max_conv=2, feas_nh=0, opt_nh=5, verbosity=2, optimizer_constructor=\"GLPK\", optimizer_attributes=[(\"tm_lim\", 1000)])\n\nDescription\n\nAssemble the tests.\n\nArguments\n\nATAmodel::Model : Required. The model built with ATA fuctions (loadsettings!, addfriends, etc...). \nsolver : Optional. Default: \"jumpATA\". Values: \"jumpATA\", \"siman\". The solving interface to be used (JuMP or internal solver based on Simulated Annealing)\nstarting_design : Optional. Default: Matrix{Float64}(undef, 0, 0). The starting design matrix. Must be a Matrix{Float64}.\nresults_folder : Optional. Default: \"RESULTS\". The folder in which the output is stored.\n\nsiman arguments:\n\n- **`start_temp`** : Optional. Default: `0.1`. Values: `[0, Inf]`. Starting temperature, set to minimum for short journeys (if 0 worse solutions will never be accepted).\n- **`geom_temp`** : Optional. Default: `0.1`. Values: `[0, Inf)`. Decreasing geometric factor.\n- **`n_item_sample`** : Optional. Default: `1`. Values: `[1, Inf]`. Number of items to alter. Set to minimum for a shallow analysis, set to maximum for a deep analysis of the neighbourhoods.\n- **`n_test_sample`** : Optional. Default: `1`. Values: `[1, Inf]`. Number of tests to alter. Set to minimum for a shallow analysis, set to maximum for a deep analysis of the neighbourhoods.\n- **`opt_feas`** : Optional. Default: `0.0`. Values: `[0, Inf)`. Optimality/feasibility balancer, if 0 only feasibility of solution is analysed. Viceversa, if 1, only optimality is considered (uncontrained model). All the other values are accepted but produce uninterpretable results.\n- **`n_fill`** : Optional. Default: `1`. Values: `[0, Inf)`. Number of fill-up phases, usually 1 is sufficient, if start_temp is high it can be higher. If a starting_design is supplied, it can be set to 0.\n- **`verbosity`** : Optional. Default: `2`. Values: `1` (minimal), `2` (detailed). Verbosity level. In the console '+' stands for improvement, '_' for accepting worse solution. The dots are the fill-up improvement steps.\n* Termination criteria: \n- **`max_time`** : Optional. Default: `1000.0`. Values: `[0, Inf)`. Time limit in seconds.\n- **`max_conv`** : Optional. Default: `2`. Values: `[1, Inf)`. Maximum convergence, stop when, after max_conv rounds no improvements have been found. Set to minimum for shallow analysis, increase it for deep analysis of neighbourhoods.\n- **`feas_nh`** : Optional. Default: `0`. Values: `[1, Inf)`. Maximum number of Feasibility neighbourhoods to explore, set to the minimum if the model is small or not highly constrained.\n- **`opt_nh`** : Optional. Default: `5`. Values: `[1, Inf)`. Maximum number of Optimality neighbourhoods to explore, set to the minimum if the model is highly constrained.\n\njumpATA arguments:\n\n- **`optimizer_constructor`** : Optional. Default: `\"GLPK\"`. Values: `\"GLPK\"`, `\"Knitro\"`, `\"Gurobi\"`, `\"Cbc\"`, `\"CPLEX\"`, `\"Xpress\"`, `\"SCIP\"`, `\"Juniper\"`. JuMP solver selection. Remember to load the required package before assemble!.\n- **`optimizer_attributes`** : Optional. Default: `[(\"tm_lim\", 1000)]`. Values: An array of pairs `(attribute, value)`. Attributes and related values for the JuMP solver.\n\nExamples\n\njulia> nothing\n\n\n\n\n\n","category":"method"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"CurrentModule = ATA\nDocTestSetup = quote\n    using ATA\nend","category":"page"},{"location":"#ATA.jl:-Automated-Test-Assembly-Made-Easy","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"A package for automated test assembly (ATA) written in Julia. Simulated Annealing algorithm is available for large scale ATA models. Otherwise, any MILP solver compatible with JuMP can be used. Interfaced with Dash or pure Julia code.","category":"page"},{"location":"#Objectives","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Objectives","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"no objective;\nMAXIMIN TIF (minimum between multiple ability points supported);\nChance-constrained MAXIMIN TIF.","category":"page"},{"location":"#Constraints","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Constraints","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"parallel tests (one group of tests);\nnon-parallel tests (more than one group of tests);\nminimum and/or maximum length;\nminimum and/or maximum number of items with certain features (categorical constraints);\nmaximum and/or minimum sum of numerical variables (quantitative constraints);\nmaximum and/or minimum expected scores at multiple ability point;\nminimum and/or maximum item use;\nmaximum and/or minimum mean of numerical variables (quantitative constraints); (beta)\nmaximum and/or minimum expected score based on the IRT paradigm or by a given column in the pool dataframe;\nmaximum overlap between tests;\ngroup by units (friend sets);\nitems exclusivity (enemy sets);","category":"page"},{"location":"#Solvers","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Solvers","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"for objectives 1 and 2: JuMP MILP solvers (see the list at https://jump.dev/JuMP.jl/v0.21.1/installation/). Follow the installation manuals of the solver you want to use. For open-source we suggest CBC or GLPK (default). For commercial we suggest Gurobi or CPLEX.\nfor objectives 1, 2 and 3: Simulated Annealing solver, pure Julia ATA solver.","category":"page"},{"location":"#Report","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Report","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Summarizing features of the assembled tests and Plots of the ICFs and TIFs are available.","category":"page"},{"location":"#How-to","page":"ATA.jl: Automated Test Assembly Made Easy","title":"How to","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"If you want to play with this package:","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Install Julia-1.5.1 at https://julialang.org/downloads/","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"run Julia.exe and type:","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"] add https://github.com/giadasp/ATA.jl","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Load the package by","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"using ATA","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Play with the test files in folder \"examples\".","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Set all the specifications by modifying the files in \"examples\" folder and run \"example.jl\" or \"example with custom objective function.jl\" following the manual in the comments.","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Distributed analysis of the neighbourhoods is available, look at \"examples/parallelExample.jl\".","category":"page"},{"location":"utils/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"utils/","page":"Utilities","title":"Utilities","text":"Modules = [ATA]\nPages   = [\"utils.jl\"]","category":"page"},{"location":"utils/#ATA.FS_to_items-Tuple{Array{Float64,1},Int64,Array{Array{Int64,1},1}}","page":"Utilities","title":"ATA.FS_to_items","text":"FS_to_items(xₜ, n_items, FS_items)\n\nDescription\n\nUngroup items grouped by friend sets. \n\nArguments\n\nxₜ::Vector{Float64} : Required. grouped items 0-1 vector.\nn_items::Int64 : Required. Number if items.\nFS_items::Vector{Vector{Int64}} : Required. vector of items included in the friend sets.\n\nOutput\n\nx_Iₜ::Vector{Float64} Returns a 0-1 vector at item level.\n\n\n\n\n\n","category":"method"},{"location":"utils/#ATA.item_char-Tuple{DataFrames.DataFrame,Array{Float64,1}}","page":"Utilities","title":"ATA.item_char","text":"item_char(pars, theta; model = \"2PL\", parametrization = \"at-b\", D = 1, derivatives = false)\n\nDescription\n\nCompute the item characteristic function (probability of correct answer).\n\nArguments\n\npars::DataFrames.DataFrame : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).\ntheta::Vector{Float64} : Required. Vector of ability points.\nmodel : Optional. Default: \"2PL\". Values:  \"1PL\", \"2PL\", \"3PL\". IRT model.\nparametrization : Optional. Default: \"at-b\". Values:  \"at-b\", \"at-ab\", \"at+b\", \"at+ab\". IRT model parametrization. Ex: at-b is Da(theta-b).\nD : Optional. Default: 1. \nderivatives : Optional. Default: false. If it's true compute the derivatives. Ow a zeros(0,0) matrix is returned. \n\nOutput\n\np::Matrix{Float64}: Matrix of probabilites. \npder::Matrix{Float64}: Matrix of derivatives.\n\n\n\n\n\n","category":"method"},{"location":"utils/#ATA.item_info-Tuple{DataFrames.DataFrame,Any}","page":"Utilities","title":"ATA.item_info","text":"item_info(pars, theta; model = \"2PL\", parametrization = \"at-b\", D = 1, derivatives = false)\n\nDescription\n\nCompute the item Fisher information function.\n\nArguments\n\npars::DataFrames.DataFrame : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).\ntheta::Vector{Float64} : Required. Vector of ability points.\nmodel : Optional. Default: \"2PL\". Values:  \"1PL\", \"2PL\", \"3PL\". IRT model.\nparametrization : Optional. Default: \"at-b\". Values:  \"at-b\", \"at-ab\", \"at+b\", \"at+ab\". IRT model parametrization. Ex: at-b is Da(theta-b).\nD : Optional. Default: 1. \n\nOutput\n\ni::Matrix{Float64}: Matrix of the item information. \n\n\n\n\n\n","category":"method"},{"location":"utils/#ATA.resp_gen-Tuple{DataFrames.DataFrame,Array{Float64,1},Array{Float64,2}}","page":"Utilities","title":"ATA.resp_gen","text":"resp_gen(f, pars, theta, design)\n\nDescription\n\nGenerate dichotomous responses from a dataframe of item parameters and a vector of ability points (theta).\n\nArguments\n\npars::DataFrames.DataFrame : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).\ntheta::Vector{Float64} : Required. Vector of ability points.\ndesign::Matrix{Float64} : Optional. Default: \"at-b\". Values:  \"at-b\", \"at-ab\", \"at+b\", \"at+ab\". IRT model parametrization. Ex: at-b is Da(theta-b).\n\nOutput\n\ni::Matrix{Float64}: Matrix of the item information. \n\n\n\n\n\n","category":"method"},{"location":"build/#Build-the-ATA-model","page":"Build the ATA model","title":"Build the ATA model","text":"","category":"section"},{"location":"build/","page":"Build the ATA model","title":"Build the ATA model","text":"Modules = [ATA]\nPages   = [\"build.jl\"]","category":"page"}]
}
