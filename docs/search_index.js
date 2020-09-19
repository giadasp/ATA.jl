var documenterSearchIndex = {"docs":
[{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/#Example-with-JuMP-(0.21.3)-and-Cbc","page":"Examples","title":"Example with JuMP (0.21.3) and Cbc","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"0. Add and load required packages.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"#cd(\"folder in which the package is saved\")\nusing Pkg\nPkg.activate(\".\")  # required\nPkg.instantiate()\nusing ATA\ncd(\"examples\")\nPkg.add(\"JuMP@0.21.3\")\nusing JuMP\nPkg.add(\"Cbc\")\nusing Cbc","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"1. Resetting the ATA process (Needed)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"ATAmodel = start_ATA()","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"2. Add file with custom settings (Needed)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"@info load_settings!(ATAmodel; settings_file=\"settingsATA.jl\", bank_file=\"data/bank.csv\", bank_delim=\";\")[2]","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"3. Add friend set variables (Optional)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"@info add_friends!(ATAmodel)[2]","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"4. Add enemy set variables (Optional)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"@info add_enemies!(ATAmodel)[2]","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"5. Add categorical constraints (Optional)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"@info add_constraints!(ATAmodel; constraints_file=\"constraints.csv\", constraints_delim=\";\")[2]","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"6. Add overlap maxima (Optional)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"@info add_overlap!(ATAmodel; overlap_file=\"Overlap Matrix.csv\", overlap_delim=\";\")[2]","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"7. Add expected score constraints (Optional)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"@info add_exp_score!(ATAmodel)[2]","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"8. Add overlap maxima (Optional, Needed if add_friends!(model) hase been run)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"@info group_by_friends!(ATAmodel)[2]","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"9. Add objective function (Optional)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"@info add_obj_fun!(ATAmodel)[2] ","category":"page"},{"location":"examples/#Assembly-settings","page":"Examples","title":"Assembly settings","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Set the solver, \"siman\" for simulated annealing, \"jumpATA\" for MILP solver.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"solver = \"jumpATA\"","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"MILP (Not suggested for large scale ATA) Select the solver, Cbc as open-source is a good option:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"optimizer_constructor = \"Cbc\"","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Optimizer attributes:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"optimizer_attributes = [(\"seconds\", 100), (\"logLevel\", 1)]","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"10. assemble","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"  assemble!(ATAmodel;\n      solver=solver,\n      optimizer_attributes=optimizer_constructor,\n      optimizer_constructor=optimizer_attributes\n  )","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"All the settings and outputs from optimization are in ATAmodel object. See the struct in ATA.jl to understand how to retrieve all the information. A summary of the resulting tests is available in <results_folder>/Results.txt. To print the plots the packages Plots and PGFPlotsX are required together with an installed implementation of LaTeX such as MikTein.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"\nPkg.add(\"Plots\")\nusing Plots\nPkg.add(\"PGFPlotsX\")\nusing PGFPlotsX\n\nprint_results!(ATAmodel;\ngroup_by_fs=true,\nplots_out=true,\nresults_folder=\"RESULTS\")","category":"page"},{"location":"opt/#Assemble-tests","page":"Assemble tests","title":"Assemble tests","text":"","category":"section"},{"location":"opt/","page":"Assemble tests","title":"Assemble tests","text":"Modules = [ATA]\nPages   = [\"opt.jl\"]","category":"page"},{"location":"opt/#ATA.assemble!-Tuple{ATA.Model}","page":"Assemble tests","title":"ATA.assemble!","text":"assemble!(ATAmodel::Model; solver=\"jumpATA\", starting_design=Matrix{Float64}(undef, 0, 0), results_folder=\"RESULTS\", start_temp=0.1, geom_temp=0.1, n_item_sample=1, n_test_sample=1, opt_feas=0.0, n_fill=1, max_time=1000.00, max_conv=2, feas_nh=0, opt_nh=5, verbosity=2, optimizer_constructor=\"GLPK\", optimizer_attributes=[(\"tm_lim\", 1000)])\n\nDescription\n\nAssemble the tests.\n\nArguments\n\nATAmodel::Model : Required. The model built with ATA fuctions (loadsettings!, addfriends, etc...). \nsolver : Optional. Default: \"jumpATA\". Values: \"jumpATA\", \"siman\". The solving interface to be used (JuMP or internal solver based on Simulated Annealing).\nstarting_design : Optional. Default: Matrix{Float64}(undef, 0, 0). The starting design matrix. Must be a Matrix{Float64}.\nresults_folder : Optional. Default: \"RESULTS\". The folder in which the output is stored.\n\nsiman arguments\n\nstart_temp : Optional. Default: 0.1. Values: [0, Inf]. Starting temperature, set to minimum for short journeys (if 0 worse solutions will never be accepted).\ngeom_temp : Optional. Default: 0.1. Values: [0, Inf). Decreasing geometric factor.\nn_item_sample : Optional. Default: 1. Values: [1, Inf]. Number of items to alter. Set to minimum for a shallow analysis, set to maximum for a deep analysis of the neighbourhoods.\nn_test_sample : Optional. Default: 1. Values: [1, Inf]. Number of tests to alter. Set to minimum for a shallow analysis, set to maximum for a deep analysis of the neighbourhoods.\nopt_feas : Optional. Default: 0.0. Values: [0, Inf). Optimality/feasibility balancer, if 0 only feasibility of solution is analysed. Viceversa, if 1, only optimality is considered (uncontrained model). All the other values are accepted but produce uninterpretable results.\nn_fill : Optional. Default: 1. Values: [0, Inf). Number of fill-up phases, usually 1 is sufficient, if starttemp is high it can be higher. If a startingdesign is supplied, it can be set to 0.\nverbosity : Optional. Default: 2. Values: 1 (minimal), 2 (detailed). Verbosity level. In the console '+' stands for improvement, '_' for accepting worse solution. The dots are the fill-up improvement steps.\nTermination criteria\nmax_time : Optional. Default: 1000.0. Values: [0, Inf). Time limit in seconds.\nmax_conv : Optional. Default: 2. Values: [1, Inf). Maximum convergence, stop when, after max_conv rounds no improvements have been found. Set to minimum for shallow analysis, increase it for deep analysis of neighbourhoods.\nfeas_nh : Optional. Default: 0. Values: [1, Inf). Maximum number of Feasibility neighbourhoods to explore, set to the minimum if the model is small or not highly constrained.\nopt_nh : Optional. Default: 5. Values: [1, Inf). Maximum number of Optimality neighbourhoods to explore, set to the minimum if the model is highly constrained.\n\njumpATA arguments\n\noptimizer_constructor : Optional. Default: \"GLPK\". Values: \"GLPK\", \"Knitro\", \"Gurobi\", \"Cbc\", \"CPLEX\", \"Xpress\", \"SCIP\", \"Juniper\". JuMP solver selection. Remember to load the required package before assemble!.\noptimizer_attributes : Optional. Default: [(\"tm_lim\", 1000)]. Values: An array of pairs (attribute, value). Attributes and related values for the JuMP solver.\n\n\n\n\n\n","category":"method"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"CurrentModule = ATA\nDocTestSetup = quote\n    using ATA\nend","category":"page"},{"location":"#ATA.jl:-Automated-Test-Assembly-Made-Easy","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"A package for automated test assembly (ATA) written in Julia. Simulated Annealing algorithm is available for large scale ATA models. Otherwise, any MILP solver compatible with JuMP can be used. Interfaced with Dash or pure Julia code.","category":"page"},{"location":"#Objectives","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Objectives","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"no objective;\nMAXIMIN TIF (minimum between multiple ability points supported);\nChance-constrained MAXIMIN TIF.","category":"page"},{"location":"#Constraints","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Constraints","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"parallel tests (one group of tests);\nnon-parallel tests (more than one group of tests);\nminimum and/or maximum length;\nminimum and/or maximum number of items with certain features (categorical constraints);\nmaximum and/or minimum sum of numerical variables (quantitative constraints);\nmaximum and/or minimum expected scores at multiple ability point;\nminimum and/or maximum item use;\nmaximum and/or minimum mean of numerical variables (quantitative constraints); (beta)\nmaximum and/or minimum expected score based on the IRT paradigm or by a given column in the pool dataframe;\nmaximum overlap between tests;\ngroup by units (friend sets);\nitems exclusivity (enemy sets);","category":"page"},{"location":"#Solvers","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Solvers","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"for objectives 1 and 2: JuMP MILP solvers (see the list at https://jump.dev/JuMP.jl/v0.21.1/installation/). Follow the installation manuals of the solver you want to use. For open-source we suggest CBC or GLPK (default). For commercial we suggest Gurobi or CPLEX.\nfor objectives 1, 2 and 3: Simulated Annealing solver, pure Julia ATA solver.","category":"page"},{"location":"#Report","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Report","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Summarizing features of the assembled tests and Plots of the ICFs and TIFs are available.","category":"page"},{"location":"#How-to","page":"ATA.jl: Automated Test Assembly Made Easy","title":"How to","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"If you want to play with this package:","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Install Julia-1.5.1 at https://julialang.org/downloads/","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"run Julia.exe and type:","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"] add https://github.com/giadasp/ATA.jl","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Load the package by","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"using ATA","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Play with the test files in folder \"examples\".","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Set all the specifications by modifying the files in \"examples\" folder and run \"example.jl\" or \"example custom objective function.jl\" following the manual in the comments.","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Distributed analysis of the neighbourhoods is available, look at \"examples/example parallel.jl\".","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"A Dash app is available for an even more easier experience. Look at \"examples/example Dash app.jl\"","category":"page"},{"location":"utils/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"utils/","page":"Utilities","title":"Utilities","text":"Modules = [ATA]\nPages   = [\"utils.jl\"]","category":"page"},{"location":"utils/#ATA.FS_to_items-Tuple{Array{Float64,1},Int64,Array{Array{Int64,1},1}}","page":"Utilities","title":"ATA.FS_to_items","text":"FS_to_items(xₜ, n_items, FS_items)\n\nDescription\n\nUngroup items grouped by friend sets. \n\nArguments\n\nxₜ::Vector{Float64} : Required. grouped items 0-1 vector.\nn_items::Int64 : Required. Number if items.\nFS_items::Vector{Vector{Int64}} : Required. vector of items included in the friend sets.\n\nOutput\n\nx_Iₜ::Vector{Float64} Returns a 0-1 vector at item level.\n\n\n\n\n\n","category":"method"},{"location":"utils/#ATA.item_char-Tuple{DataFrames.DataFrame,Array{Float64,1}}","page":"Utilities","title":"ATA.item_char","text":"item_char(pars, theta; model = \"2PL\", parametrization = \"at-b\", D = 1, derivatives = false)\n\nDescription\n\nCompute the item characteristic function (probability of correct answer).\n\nArguments\n\npars::DataFrames.DataFrame : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).\ntheta::Vector{Float64} : Required. Vector of ability points.\nmodel : Optional. Default: \"2PL\". Values:  \"1PL\", \"2PL\", \"3PL\". IRT model.\nparametrization : Optional. Default: \"at-b\". Values:  \"at-b\", \"at-ab\", \"at+b\", \"at+ab\". IRT model parametrization. Ex: at-b is Da(theta-b).\nD : Optional. Default: 1. \nderivatives : Optional. Default: false. If it's true compute the derivatives. Ow a zeros(0,0) matrix is returned. \n\nOutput\n\np::Matrix{Float64}: Matrix of probabilites. \npder::Matrix{Float64}: Matrix of derivatives.\n\n\n\n\n\n","category":"method"},{"location":"utils/#ATA.item_info-Tuple{DataFrames.DataFrame,Any}","page":"Utilities","title":"ATA.item_info","text":"item_info(pars, theta; model = \"2PL\", parametrization = \"at-b\", D = 1, derivatives = false)\n\nDescription\n\nCompute the item Fisher information function.\n\nArguments\n\npars::DataFrames.DataFrame : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).\ntheta::Vector{Float64} : Required. Vector of ability points.\nmodel : Optional. Default: \"2PL\". Values:  \"1PL\", \"2PL\", \"3PL\". IRT model.\nparametrization : Optional. Default: \"at-b\". Values:  \"at-b\", \"at-ab\", \"at+b\", \"at+ab\". IRT model parametrization. Ex: at-b is Da(theta-b).\nD : Optional. Default: 1. \n\nOutput\n\ni::Matrix{Float64}: Matrix of the item information. \n\n\n\n\n\n","category":"method"},{"location":"utils/#ATA.resp_gen-Tuple{DataFrames.DataFrame,Array{Float64,1},Array{Float64,2}}","page":"Utilities","title":"ATA.resp_gen","text":"resp_gen(f, pars, theta, design)\n\nDescription\n\nGenerate dichotomous responses from a dataframe of item parameters and a vector of ability points (theta).\n\nArguments\n\npars::DataFrames.DataFrame : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).\ntheta::Vector{Float64} : Required. Vector of ability points.\ndesign::Matrix{Float64} : Optional. Default: \"at-b\". Values:  \"at-b\", \"at-ab\", \"at+b\", \"at+ab\". IRT model parametrization. Ex: at-b is Da(theta-b).\n\nOutput\n\ni::Matrix{Float64}: Matrix of the item information. \n\n\n\n\n\n","category":"method"},{"location":"build/#Build-the-ATA-model","page":"Build the ATA model","title":"Build the ATA model","text":"","category":"section"},{"location":"build/","page":"Build the ATA model","title":"Build the ATA model","text":"start_ATA()\nload_settings!(::ATA.Model; settings_file = \"SettingsATA.jl\", bank_file = \"bank.csv\", bank_delim = \";\")\nadd_friends!(::ATA.Model)\nadd_enemies!(::ATA.Model)\nadd_constraints!(::ATA.Model; constraints_file = \"Constraints.csv\", constraints_delim = \";\")\nadd_overlap!(::ATA.Model; overlap_file = \"OverlapMatrix.csv\", overlap_delim=\";\")\nadd_exp_score!(::ATA.Model)\nadd_obj_fun!(::ATA.Model)\ngroup_by_friends!(::ATA.Model)\nload_design!(::Matrix{Any}, ::ATA.Model)","category":"page"},{"location":"build/#ATA.start_ATA-Tuple{}","page":"Build the ATA model","title":"ATA.start_ATA","text":"start_ATA()\n\nDescription\n\nStart the test assembly instance.\n\nOutput\n\nIt returns an empty ATA.Model object.\n\n\n\n\n\n","category":"method"},{"location":"build/#ATA.load_settings!-Tuple{ATA.Model}","page":"Build the ATA model","title":"ATA.load_settings!","text":"load_settings!(ATAmodel::Model; settings_file = \"SettingsATA.jl\", bank_file = \"bank.csv\", bank_delim = \";\")\n\nDescription\n\nLoad the test assembly settings contained in the settings_file. It requires the start_ATA step.\n\nArguments\n\nATAmodel::Model : Required. The model built with start_ATA() function.\nsettings_file : Optional. Default: \"SettingsATA.jl\". The path of the file containing the ATA settings in the form of an InputSettings struct.\nbank_file : Optional. Default: \"bank.csv\". The path of the file containing the item pool/bank in the form of custom-separated values.\nbank_delim : Optional. Default: \";\". The custom-separator for the bank_file.\n\n\n\n\n\n","category":"method"},{"location":"build/#ATA.add_friends!-Tuple{ATA.Model}","page":"Build the ATA model","title":"ATA.add_friends!","text":"add_friends!(ATAmodel::Model)\n\nDescription\n\nAdd friend sets to the ATA.Model as specified in the settings_file loaded by load_settings! function. It requires the load_settings! step.\n\nArguments\n\nATAmodel::Model : Required. The model built with start_ATA() and with settings loaded by load_settings! function.\n\n\n\n\n\n","category":"method"},{"location":"build/#ATA.add_enemies!-Tuple{ATA.Model}","page":"Build the ATA model","title":"ATA.add_enemies!","text":"add_enemies!(ATAmodel::Model)\n\nDescription\n\nAdd enemy sets to the ATA.Model as specified in the settings_file loaded by load_settings! function. It requires the load_settings! step.\n\nArguments\n\nATAmodel::Model : Required. The model built with start_ATA() and with settings loaded by load_settings! function.\n\n\n\n\n\n","category":"method"},{"location":"build/#ATA.add_constraints!-Tuple{ATA.Model}","page":"Build the ATA model","title":"ATA.add_constraints!","text":"add_constraints!(ATAmodel::Model; constraints_file = \"Constraints.csv\", constraints_delim = \";\")\n\nDescription\n\nAdd categorical and sum constraints to the ATA.Model as specified in the constraints_file. It requires the load_settings! step.\n\nArguments\n\nATAmodel::Model : Required. The model built with start_ATA() and with settings loaded by load_settings! function.\nconstraints_file : Optional. Default: \"Constraints.csv\". The path of the file containing the categorical and sum constraints in the form of custom-separated values.\nconstraints_delim : Optional. Default: \";\". The custom-separator for the constraints_file.\n\n\n\n\n\n","category":"method"},{"location":"build/#ATA.add_overlap!-Tuple{ATA.Model}","page":"Build the ATA model","title":"ATA.add_overlap!","text":"add_overlap!(ATAmodel::Model; overlap_file = \"OverlapMatrix.csv\", overlap_delim=\";\")\n\nDescription\n\nAdd maximum overlap constraints to the ATA.Model as specified by the overlap_file which contains a TxT matrix defining the maximum number of shared items between each pair of tests. It requires the load_settings! step.\n\nArguments\n\nATAmodel::Model : Required. The model built with start_ATA() and with settings loaded by load_settings! function.\noverlap_file : Optional. Default: \"OverlapMatrix.csv\". The path of the file containing the maximum overlap between each pair of tests in the form of a matrix with custom-separated values.\noverlap_delim : Optional. Default: \";\". The custom-separator for the overlap_file.\n\n\n\n\n\n","category":"method"},{"location":"build/#ATA.add_exp_score!-Tuple{ATA.Model}","page":"Build the ATA model","title":"ATA.add_exp_score!","text":"add_exp_score!(ATAmodel::Model)\n\nDescription\n\nAdd the expected score constraints as specified in the settings_file. If the names of the expected score columns of the item bank are not provided in the settings_file they are computed as the item characteristic functions following the IRT model always specified in the settings_file. It requires the load_settings! step.\n\nArguments\n\nATAmodel::Model : Required. The model built with start_ATA() and with settings loaded by load_settings! function.\n\n\n\n\n\n","category":"method"},{"location":"build/#ATA.add_obj_fun!-Tuple{ATA.Model}","page":"Build the ATA model","title":"ATA.add_obj_fun!","text":"add_obj_fun!(ATAmodel::Model)\n\nDescription\n\nAdd the objective function as specified in the settings_file. It requires the load_settings! step.  \n\nArguments\n\nATAmodel::Model : Required. The model built with start_ATA() and with settings loaded by load_settings! function.\n\n\n\n\n\n","category":"method"},{"location":"build/#ATA.group_by_friends!-Tuple{ATA.Model}","page":"Build the ATA model","title":"ATA.group_by_friends!","text":"group_by_friends!(ATAmodel::Model)\n\nDescription\n\nGroup the items by friend sets once the friend sets have been added to the ATA.Model by add_friends! function and as last operation of the assembly.\n\nArguments\n\nATAmodel::Model : Required. The model built with start_ATA(), settings loaded by load_settings! function and friend sets added by add_friends! function.\n\n\n\n\n\n","category":"method"},{"location":"build/#ATA.load_design!-Tuple{Array{Any,2},ATA.Model}","page":"Build the ATA model","title":"ATA.load_design!","text":"load_design!(design::Matrix{Any}, ATAmodel::Model)\n\nDescription\n\nLoad a 0-1 IxT (or nFSxT if the items are grouped by friend sets) design matrix into the ATA model. Useful for loading a custom starting design before the assemble! step or to print/plot the features of the tests produced by a custom design before running print_results!\n\nArguments\n\nATAmodel::Model : Required. The model built with start_ATA(), settings loaded by load_settings!.\n\n\n\n\n\n","category":"method"}]
}
