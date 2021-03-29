var documenterSearchIndex = {"docs":
[{"location":"print_and_plot/#Print-and-Plot","page":"Print and Plot","title":"Print and Plot","text":"","category":"section"},{"location":"print_and_plot/","page":"Print and Plot","title":"Print and Plot","text":"CurrentModule = ATA","category":"page"},{"location":"print_and_plot/#Print-Results-on-Text-File","page":"Print and Plot","title":"Print Results on Text File","text":"","category":"section"},{"location":"print_and_plot/","page":"Print and Plot","title":"Print and Plot","text":"print_results(\n        ata_model::NoObjModel;\n        group_by_fs = false,\n        results_folder = \"results\",\n        sim_pool::DataFrames.DataFrame = DataFrame(),\n)\nprint_results(\n        ata_model::MaximinModel;\n        group_by_fs = false,\n        results_folder = \"results\",\n        sim_pool::DataFrames.DataFrame = DataFrame(),\n)\nprint_results(\n        ata_model::MinimaxModel;\n        group_by_fs = false,\n        results_folder = \"results\",\n        sim_pool::DataFrames.DataFrame = DataFrame(),\n)\nprint_results(\n    ata_model::CCMaximinModel;\n    group_by_fs = false,\n    results_folder = \"results\",\n    sim_pool::DataFrames.DataFrame = DataFrame(),\n)","category":"page"},{"location":"print_and_plot/#ATA.print_results-Tuple{ATA.NoObjModel}","page":"Print and Plot","title":"ATA.print_results","text":"print_results(\n    ata_model::NoObjModel;\n    group_by_fs = false,\n    results_folder = \"results\",\n    sim_pool::DataFrames.DataFrame = DataFrame(),\n)\n\nDescription\n\nPrint the features of the assembled tests.\n\nArguments\n\nata_model::NoObjModel : Required. The model built with ATA fuctions, ata_model.design matrix must be IxT or nfsxT if the items are grouped by friend sets. \ngroup_by_fs : Optional. Default: false. Set to true if items have been grouped by friend sets by group_by_friends!.\nresults_folder : Optional. Default: \"results\". The folder in which the output is stored.\nsim_pool::DataFrames.DataFrame : Optional. Default: DataFrame(). The pool with true item paramaters. For simulation studies.\n\n\n\n\n\n","category":"method"},{"location":"print_and_plot/#ATA.print_results-Tuple{ATA.MaximinModel}","page":"Print and Plot","title":"ATA.print_results","text":"print_results(\n    ata_model::Union{MaximinModel, SoysterMaximinModel, DeJongMaximinModel};\n    group_by_fs = false,\n    results_folder = \"results\",\n    sim_pool::DataFrames.DataFrame = DataFrame(),\n)\n\nDescription\n\nPrint the features of the assembled tests.\n\nArguments\n\nata_model::Union{MaximinModel, SoysterMaximinModel, DeJongMaximinModel} : Required. The model built with ATA fuctions, ata_model.design matrix must be IxT or nfsxT if the items are grouped by friend sets. \ngroup_by_fs : Optional. Default: false. Set to true if items have been grouped by friend sets by group_by_friends!.\nresults_folder : Optional. Default: \"results\". The folder in which the output is stored.\nsim_pool::DataFrames.DataFrame : Optional. Default: DataFrame(). The pool with true item paramaters. For simulation studies.\n\n\n\n\n\n","category":"method"},{"location":"print_and_plot/#ATA.print_results-Tuple{ATA.CCMaximinModel}","page":"Print and Plot","title":"ATA.print_results","text":"print_results(\n    ata_model::CCMaximinModel;\n    group_by_fs = false,\n    results_folder = \"results\",\n    sim_pool::DataFrames.DataFrame = DataFrame(),\n)\n\nDescription\n\nPrint the features of the assembled tests.\n\nArguments\n\nata_model::CCMaximinModel : Required. The model built with ATA fuctions, ata_model.design matrix must be IxT or nfsxT if the items are grouped by friend sets. \ngroup_by_fs : Optional. Default: false. Set to true if items have been grouped by friend sets by group_by_friends!.\nresults_folder : Optional. Default: \"results\". The folder in which the output is stored.\nsim_pool::DataFrames.DataFrame : Optional. Default: DataFrame(). The pool with true item paramaters. For simulation studies.\n\n\n\n\n\n","category":"method"},{"location":"print_and_plot/#Plot-TIFs-and-ICFs","page":"Print and Plot","title":"Plot TIFs and ICFs","text":"","category":"section"},{"location":"print_and_plot/","page":"Print and Plot","title":"Print and Plot","text":"plot_results(\n    ata_model::AbstractModel;\n    group_by_fs = false,\n    results_folder = \"plots\",\n)","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"The folder \"examples\" contains several examples of ATA models solved with ATA.jl.","category":"page"},{"location":"structs/#Structs","page":"Structs","title":"Structs","text":"","category":"section"},{"location":"structs/","page":"Structs","title":"Structs","text":"Modules = [ATA]\nPages   = [\"structs/structs.jl\"]","category":"page"},{"location":"opt/#Assemble-tests","page":"Assemble tests","title":"Assemble tests","text":"","category":"section"},{"location":"opt/","page":"Assemble tests","title":"Assemble tests","text":"Modules = [ATA]\nPages   = [\"opt.jl\"]","category":"page"},{"location":"opt/#ATA.assemble!-Tuple{ATA.AbstractModel}","page":"Assemble tests","title":"ATA.assemble!","text":"assemble!(ata_model::AbstractModel; solver=\"jumpATA\", starting_design=Matrix{Float64}(undef, 0, 0), results_folder=\"results\", start_temp=0.1, geom_temp=0.1, n_item_sample=1, n_test_sample=1, opt_feas=0.0, n_fill=1, max_time=1000.00, max_conv=2, feas_nh=0, opt_nh=5, verbosity=2, optimizer_constructor=\"GLPK\", optimizer_attributes=[(\"tm_lim\", 1000)])\n\nDescription\n\nAssemble the tests.\n\nArguments\n\nata_model::AbstractModel : Required. The model built with ATA fuctions. \nsolver : Optional. Default: \"jumpATA\". Values: \"jumpATA\", \"siman\". The solving interface to be used (JuMP or internal solver based on Simulated Annealing).\nstarting_design : Optional. Default: Matrix{Float64}(undef, 0, 0). The starting design matrix. Must be a Matrix{Float64}.\nresults_folder : Optional. Default: \"results\". The folder in which the output is stored.\n\nsiman arguments\n\nstart_temp : Optional. Default: 0.1. Values: [0, Inf]. Starting temperature, set to minimum for short journeys (if 0 worse solutions will never be accepted).\ngeom_temp : Optional. Default: 0.1. Values: [0, Inf). Decreasing geometric factor.\nn_item_sample : Optional. Default: 1. Values: [1, Inf]. Number of items to alter. Set to minimum for a shallow analysis, set to maximum for a deep analysis of the neighbourhoods.\nn_test_sample : Optional. Default: 1. Values: [1, Inf]. Number of tests to alter. Set to minimum for a shallow analysis, set to maximum for a deep analysis of the neighbourhoods.\nopt_feas : Optional. Default: 0.0. Values: [0, Inf). Optimality/feasibility balancer, if 0 only feasibility of solution is analysed. Viceversa, if 1, only optimality is considered (uncontrained model). All the other values are accepted but produce uninterpretable results.\nn_fill : Optional. Default: 1. Values: [0, Inf). Number of fill-up phases, usually 1 is sufficient, if starttemp is high it can be higher. If a startingdesign is supplied, it can be set to 0.\nverbosity : Optional. Default: 2. Values: 1 (minimal), 2 (detailed). Verbosity level. In the console '+' stands for improvement, '_' for accepting worse solution. The dots are the fill-up improvement steps.\nTermination criteria\nmax_time : Optional. Default: 1000.0. Values: [0, Inf). Time limit in seconds.\nmax_conv : Optional. Default: 2. Values: [1, Inf). Maximum convergence, stop when, after max_conv rounds no improvements have been found. Set to minimum for shallow analysis, increase it for deep analysis of neighbourhoods.\nfeas_nh : Optional. Default: 0. Values: [1, Inf). Maximum number of Feasibility neighbourhoods to explore, set to the minimum if the model is small or not highly constrained.\nopt_nh : Optional. Default: 5. Values: [1, Inf). Maximum number of Optimality neighbourhoods to explore, set to the minimum if the model is highly constrained.\n\njumpATA arguments\n\noptimizer_constructor : Optional. Default: \"GLPK\". Values: \"GLPK\", \"Knitro\", \"Gurobi\", \"Cbc\", \"CPLEX\", \"Xpress\", \"SCIP\", \"Juniper\". JuMP solver selection. Remember to load the required package before assemble!.\noptimizer_attributes : Optional. Default: [(\"tm_lim\", 1000)]. Values: An array of pairs (attribute, value). Attributes and related values for the JuMP solver.\n\nother keyword arguments\n\nkwargs... : Optional. \n\n\n\n\n\n","category":"method"},{"location":"compact/#Compact-ATA","page":"Compact ATA","title":"Compact ATA","text":"","category":"section"},{"location":"compact/","page":"Compact ATA","title":"Compact ATA","text":"Modules = [ATA]\nPages   = [\"compact_ata.jl\"]","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"CurrentModule = ATA\nDocTestSetup = quote\n    using ATA\nend","category":"page"},{"location":"#ATA.jl:-Automated-Test-Assembly-Made-Easy","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Version 0.10.0","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"A package for automated test assembly (ATA) written in Julia. Simulated Annealing algorithm is available for large scale ATA models. Otherwise, any MILP solver compatible with JuMP can be used. Interfaced with Dash or pure Julia code.","category":"page"},{"location":"#Documentation","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Documentation","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Documentation on exported functions available at link (in progress).","category":"page"},{"location":"#Objectives","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Objectives","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"no objective;\nMAXIMIN TIF (minimum between multiple ability points supported);\nChance-constrained MAXIMIN TIF;\nMINIMAX TIF (maximum between multiple ability points supported);\ncustom objective function (only for siman solver)","category":"page"},{"location":"#Constraints","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Constraints","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"parallel tests (one group of tests);\nnon-parallel tests (more than one group of tests);\nminimum and/or maximum length;\nminimum and/or maximum number of items with certain features (categorical constraints);\nmaximum and/or minimum sum of numerical variables (quantitative constraints);\nmaximum and/or minimum expected scores at multiple ability point; \nmaximum and/or minimum expected score based on the IRT paradigm or by a given column in the pool dataframe;\nminimum and/or maximum item use;\nmaximum and/or minimum mean of numerical variables (quantitative constraints); (will be implemented soon)\nmaximum overlap between tests; (increases dramatically the size and complexity of MILP ATA models :boom:, we suggest to avoid this constraint if the ATA model is already very large or, alternatively, to use the siman solver :smirk_cat:)\ngroup by units (friend sets);\nitems exclusivity (enemy sets);","category":"page"},{"location":"#Solvers","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Solvers","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"for objectives 1, 2, and 4: JuMP MILP solvers (see the list at https://jump.dev/JuMP.jl/v0.21.1/installation/). Follow the installation manuals of the solver you want to use. We suggest CBC or GLPK (default) as open-source solvers. Also commercial solvers are supported, CPLEX and Gurobi are examples.\nfor objectives 1, 2, 3, 4 and 5 (all): Simulated Annealing solver (siman), pure Julia ATA solver.","category":"page"},{"location":"#Report","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Report","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Summarizing features of the assembled tests and plots are available. Add ATAPlot to plot the Test Information Functions (TIFs) and Item Characteristic Functions (ICTs).","category":"page"},{"location":"#How-to","page":"ATA.jl: Automated Test Assembly Made Easy","title":"How to","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"If you want to play with this package:","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Install Julia-1.6.0-rc1 at https://julialang.org/downloads/","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"run Julia.exe and type:","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"] add https://github.com/giadasp/ATA.jl","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Load the package by","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"using ATA","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Play with the test files in folder \"examples\".","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"For a quick and compact ATA instance, run the code in \"example_compact.jl\".","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"If you want to dig in all the ATA building, assembly and output steps, run the code in the other files in the folder. The comments explain which arguments must be customized in order to solve your ATA model.  If you do not modify the arguments you just solve a toy ATA model with 3 non parallel groups of tests assembled starting from a 366 items bank.","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"For an even more easier ATA experience, look at \"exampledashapp.jl\". It requires the installation of the package ATADash","category":"page"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"As a remind, to plot the TIFs and ICTs you must install (add) and load (using) the package ATAPlot.","category":"page"},{"location":"#Documentation-Contents","page":"ATA.jl: Automated Test Assembly Made Easy","title":"Documentation Contents","text":"","category":"section"},{"location":"","page":"ATA.jl: Automated Test Assembly Made Easy","title":"ATA.jl: Automated Test Assembly Made Easy","text":"Pages = [\"build.md\", \"opt.md\", \"compact.md\", \"print_and_plot.md\", \"structs.md\", \"utils.md\", \"examples.md\"]\nDepth = 3","category":"page"},{"location":"utils/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"utils/","page":"Utilities","title":"Utilities","text":"Modules = [ATA]\nPages   = [\"utils.jl\"]","category":"page"},{"location":"utils/#ATA.fs_to_items-Tuple{Vector{Float64}, Int64, Vector{Vector{Int64}}}","page":"Utilities","title":"ATA.fs_to_items","text":"fs_to_items(xₜ, n_items, fs_items)\n\nDescription\n\nUngroup items grouped by friend sets. \n\nArguments\n\nxₜ::Vector{Float64} : Required. grouped items 0-1 vector.\nn_items::Int64 : Required. Number if items.\nfs_items::Vector{Vector{Int64}} : Required. vector of items included in the friend sets.\n\nOutput\n\nx_Iₜ::Vector{Float64} Returns a 0-1 vector at item level.\n\n\n\n\n\n","category":"method"},{"location":"utils/#ATA.item_char-Tuple{DataFrames.DataFrame, Float64}","page":"Utilities","title":"ATA.item_char","text":"item_char(pars::DataFrames.DataFrame, theta::Float64; model = \"2PL\", parametrization = \"at-b\", D = 1, derivatives = false)\n\nDescription\n\nCompute the item characteristic function (probability of correct answer).\n\nArguments\n\npars::DataFrames.DataFrame : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).\ntheta::Float64 : Required. Ability point.\nmodel : Optional. Default: \"2PL\". Values:  \"1PL\", \"2PL\", \"3PL\". IRT model.\nparametrization : Optional. Default: \"at-b\". Values:  \"at-b\", \"at-ab\", \"at+b\", \"at+ab\". IRT model parametrization. Ex: at-b is Da(theta-b).\nD : Optional. Default: 1. \nderivatives : Optional. Default: false. If it's true compute the derivatives. Ow a zeros(0,0) matrix is returned. \n\nOutput\n\np::Matrix{Float64}: Matrix of probabilites. \npder::Matrix{Float64}: Matrix of derivatives.\n\n\n\n\n\n","category":"method"},{"location":"utils/#ATA.item_char-Tuple{DataFrames.DataFrame, Vector{Float64}}","page":"Utilities","title":"ATA.item_char","text":"item_char(pars::DataFrames.DataFrame, theta::Vector{Float64}; model = \"2PL\", parametrization = \"at-b\", D = 1, derivatives = false)\n\nDescription\n\nCompute the item characteristic function (probability of correct answer).\n\nArguments\n\npars::DataFrames.DataFrame : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).\ntheta::Vector{Float64} : Required. Vector of ability points.\nmodel : Optional. Default: \"2PL\". Values:  \"1PL\", \"2PL\", \"3PL\". IRT model.\nparametrization : Optional. Default: \"at-b\". Values:  \"at-b\", \"at-ab\", \"at+b\", \"at+ab\". IRT model parametrization. Ex: at-b is Da(theta-b).\nD : Optional. Default: 1. \nderivatives : Optional. Default: false. If it's true compute the derivatives. Ow a zeros(0,0) matrix is returned. \n\nOutput\n\np::Matrix{Float64}: Matrix of probabilites. \npder::Matrix{Float64}: Matrix of derivatives.\n\n\n\n\n\n","category":"method"},{"location":"utils/#ATA.item_info-Tuple{DataFrames.DataFrame, Float64}","page":"Utilities","title":"ATA.item_info","text":"item_info(\n    pars::DataFrames.DataFrame,\n    theta::Float64;\n    model = \"2PL\", #1PL, 2PL, 3PL, grm\n    parametrization = \"at-b\", #\"at-b, at-ab, at+b, at+ab\"\n    D = 1,\n)\n\nDescription\n\nCompute the item Fisher information function.\n\nArguments\n\npars::DataFrames.DataFrame : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).\ntheta::Float64 : Required. Vector of ability points.\nmodel : Optional. Default: \"2PL\". Values:  \"1PL\", \"2PL\", \"3PL\". IRT model.\nparametrization : Optional. Default: \"at-b\". Values:  \"at-b\", \"at-ab\", \"at+b\", \"at+ab\". IRT model parametrization. Ex: at-b is Da(theta-b).\nD : Optional. Default: 1. \n\nOutput\n\ni::Matrix{Float64}: Matrix of the item information. \n\n\n\n\n\n","category":"method"},{"location":"utils/#ATA.item_info-Tuple{DataFrames.DataFrame, Vector{Float64}}","page":"Utilities","title":"ATA.item_info","text":"item_info(\n    pars::DataFrames.DataFrame,\n    theta::Vector{Float64};\n    model = \"2PL\", #1PL, 2PL, 3PL, grm\n    parametrization = \"at-b\", #\"at-b, at-ab, at+b, at+ab\"\n    D = 1,\n)\n\nDescription\n\nCompute the item Fisher information function.\n\nArguments\n\npars::DataFrames.DataFrame : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).\ntheta::Vector{Float64} : Required. Vector of ability points.\nmodel : Optional. Default: \"2PL\". Values:  \"1PL\", \"2PL\", \"3PL\". IRT model.\nparametrization : Optional. Default: \"at-b\". Values:  \"at-b\", \"at-ab\", \"at+b\", \"at+ab\". IRT model parametrization. Ex: at-b is Da(theta-b).\nD : Optional. Default: 1. \n\nOutput\n\ni::Matrix{Float64}: Matrix of the item information. \n\n\n\n\n\n","category":"method"},{"location":"build/#Build-the-ATA-model","page":"Build the ATA model","title":"Build the ATA model","text":"","category":"section"},{"location":"build/","page":"Build the ATA model","title":"Build the ATA model","text":"CurrentModule = ATA","category":"page"},{"location":"build/#Start-ATA-Model","page":"Build the ATA model","title":"Start ATA Model","text":"","category":"section"},{"location":"build/","page":"Build the ATA model","title":"Build the ATA model","text":"start_ata(;\n    settings::InputSettings = InputSettings(),\n    bank::DataFrames.DataFrame = DataFrames.DataFrame(),\n    settings_file = \"settings_ata.jl\",\n    bank_file = \"bank.csv\",\n    bank_delim = \";\",\n)","category":"page"},{"location":"build/#ATA.start_ata-Tuple{}","page":"Build the ATA model","title":"ATA.start_ata","text":"start_ata(;\n    settings::InputSettings = InputSettings(),\n    bank::DataFrames.DataFrame = DataFrames.DataFrame(),\n    settings_file = \"settings_ata.jl\",\n    bank_file = \"bank.csv\",\n    bank_delim = \";\",\n)\n\nDescription\n\nInitialize an empty ATA model and load the test assembly settings contained in the settings_file using item data in file bank_file which values are separated by bank_delim. Alternatively, settings object and bank dataframe can be passed using the arguments settings and bank.\n\nArguments\n\nsettings::InputSettings : Optional. Default:  InputSettings() object with custom ATA settings.\nbank::DataFrames.DataFrame : Optional. Default: DataFrame(). A dataframe containing data about the items.\nsettings_file : Optional. Default: \"settings_ata.jl\". The path of the file containing the ATA settings in the form of an InputSettings struct.\nbank_file : Optional. Default: \"bank.csv\". The path of the file containing the item pool/bank in the form of custom-separated values.\nbank_delim : Optional. Default: \";\". The custom-separator for the bank_file.\n\nOutput\n\nAn ATA model.\n\n\n\n\n\n","category":"method"},{"location":"build/#Add-Constraints","page":"Build the ATA model","title":"Add Constraints","text":"","category":"section"},{"location":"build/","page":"Build the ATA model","title":"Build the ATA model","text":"add_constraints!(\n    ata_model::AbstractModel;\n    constraints::DataFrames.DataFrame = DataFrames.DataFrame(),\n    constraints_file = \"constraints.csv\",\n    constraints_delim = \";\",\n)","category":"page"},{"location":"build/#ATA.add_constraints!-Tuple{ATA.AbstractModel}","page":"Build the ATA model","title":"ATA.add_constraints!","text":"add_constraints!(ata_model::AbstractModel; constraints_file = \"constraints.csv\", constraints_delim = \";\")\n\nDescription\n\nAdd categorical and sum constraints to the ATA.AbstractModel as specified in the constraints_file. Alternatively, the constraint dataframe can be passed using the argument constraints. It requires the start_ata step.\n\nArguments\n\nata_model::AbstractModel : Required. The model built with start_ata() and with settings loaded by start_ata function.\n`constraints::DataFrames.DataFrame = DataFrames.DataFrame()`. Optional. Default: DataFrames.DataFrame(). The dataframe containing the categorical and quantitative constraints.\nconstraints_file : Optional. Default: \"constraints.csv\". The path of the file containing the categorical and sum constraints in the form of custom-separated values.\nconstraints_delim : Optional. Default: \";\". The custom-separator for the constraints_file.\n\n\n\n\n\n","category":"method"},{"location":"build/#Add-Overlap-Matrix","page":"Build the ATA model","title":"Add Overlap Matrix","text":"","category":"section"},{"location":"build/","page":"Build the ATA model","title":"Build the ATA model","text":"add_overlap!(\n    ata_model::AbstractModel;\n    overlap::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),\n    overlap_file = \"overlap_matrix.csv\",\n    overlap_delim = \";\",\n)","category":"page"},{"location":"build/#ATA.add_overlap!-Tuple{ATA.AbstractModel}","page":"Build the ATA model","title":"ATA.add_overlap!","text":"add_overlap!(ata_model::AbstractModel; overlap_file = \"overlap_matrix.csv\", overlap_delim=\";\")\n\nDescription\n\nAdd maximum overlap constraints to the ATA.AbstractModel as specified by the overlap_file which contains a TxT matrix defining the maximum number of shared items between each pair of tests. Alternatively, the overlap matrix can be passed using the argument overlap. It requires the start_ata step.\n\nArguments\n\nata_model::AbstractModel : Required. The model built with start_ata() and with settings loaded by start_ata function.\noverlap::Matrix{Float64} = Matrix{Float64}(undef, 0, 0): Optional. \noverlap_file : Optional. Default: \"overlap_matrix.csv\". The path of the file containing the maximum overlap between each pair of tests in the form of a matrix with custom-separated values.\noverlap_delim : Optional. Default: \";\". The custom-separator for the overlap_file.\n\n\n\n\n\n","category":"method"},{"location":"build/#Group-by-Friend-Sets","page":"Build the ATA model","title":"Group by Friend Sets","text":"","category":"section"},{"location":"build/","page":"Build the ATA model","title":"Build the ATA model","text":"group_by_friends!(ata_model::AbstractModel)","category":"page"},{"location":"build/#ATA.group_by_friends!-Tuple{ATA.AbstractModel}","page":"Build the ATA model","title":"ATA.group_by_friends!","text":"group_by_friends!(ata_model::AbstractModel)\n\nDescription\n\nGroup the items by friend sets once the friend sets have been added to the ATA.AbstractModel by add_friends! function and as last operation of the assembly.\n\nArguments\n\nata_model::AbstractModel : Required. The model built with start_ata(), settings loaded by start_ata function and friend sets added by add_friends! function.\n\n\n\n\n\n","category":"method"},{"location":"build/#Add-Objective-Function","page":"Build the ATA model","title":"Add Objective Function","text":"","category":"section"},{"location":"build/","page":"Build the ATA model","title":"Build the ATA model","text":"add_obj_fun!(ata_model::Union{MaximinModel, MinimaxModel}; kwargs...)\nadd_obj_fun!(\n    ata_model::SoysterMaximinModel;\n    psychometrics = false,\n    items::Vector{Psychometrics.Item} = Psychometrics.Item[],\n    items_file = \"\",\n    kwargs...\n)\nadd_obj_fun!(\n    ata_model::DeJongMaximinModel;\n    psychometrics = false,\n    items::Vector{Psychometrics.Item} = Psychometrics.Item[],\n    items_file = \"\",\n    kwargs...\n)\nadd_obj_fun!(\n    ata_model::CCMaximinModel;\n    psychometrics = false,\n    items_file = \"\",\n    items::Vector{Psychometrics.Item} = Psychometrics.Item[],\n    kwargs...\n)","category":"page"},{"location":"build/#ATA.add_obj_fun!-Tuple{Union{ATA.MaximinModel, ATA.MinimaxModel}}","page":"Build the ATA model","title":"ATA.add_obj_fun!","text":"add_obj_fun!(ata_model::Union{MaximinModel,MinimaxModel}; kwargs...)\n\nDescription\n\nAdd the objective function as specified in the settings_file. It requires the start_ata build step.   Computes the IIFs at predefined ability points.\n\nArguments\n\nata_model::Union{MaximinModel, MinimaxModel} : Required. The model built with start_ata() and with settings loaded by start_ata function.\n\n\n\n\n\n","category":"method"},{"location":"build/#ATA.add_obj_fun!-Tuple{ATA.SoysterMaximinModel}","page":"Build the ATA model","title":"ATA.add_obj_fun!","text":"add_obj_fun!(\n    ata_model::SoysterMaximinModel;\n    psychometrics = false,\n    items::Vector{Psychometrics.Item} = Psychometrics.Item[],\n    items_file = \"\",\n    kwargs...\n)\n\nDescription\n\nAdd the objective function as specified in the settings_file. It requires the start_ata build step.   Computes the IIFs at predefined ability points using R sampled item parameters and set the minimum across all the samples as the IIFs.\n\nArguments\n\nata_model::SoysterMaximinModel : Required. The model built with start_ata() and with settings loaded by start_ata function.\n\n\n\n\n\n","category":"method"},{"location":"build/#ATA.add_obj_fun!-Tuple{ATA.DeJongMaximinModel}","page":"Build the ATA model","title":"ATA.add_obj_fun!","text":"add_obj_fun!(\n    ata_model::DeJongMaximinModel;\n    psychometrics = false,\n    items::Vector{Psychometrics.Item} = Psychometrics.Item[],\n    items_file = \"\",\n    kwargs...\n)\n\nDescription\n\nAdd the objective function as specified in the settings_file. It requires the start_ata build step.   Computes the IIFs at predefined ability points using R sampled item parameters, and set the mean minus the standard deviation across all the samples as the IIFs.\n\nArguments\n\nata_model::DeJongMaximinModel : Required. The model built with start_ata() and with settings loaded by start_ata function.\n\n\n\n\n\n","category":"method"},{"location":"build/#ATA.add_obj_fun!-Tuple{ATA.CCMaximinModel}","page":"Build the ATA model","title":"ATA.add_obj_fun!","text":"add_obj_fun!(\n    ata_model::CCMaximinModel;\n    psychometrics = false,\n    items_file = \"\",\n    items::Vector{Psychometrics.Item} = Psychometrics.Item[],\n    kwargs...\n)\n\nDescription\n\nAdd the objective function as specified in the settings_file. It requires the start_ata build step.   Computes the IIFs at predefined ability points using R sampled item parameters.\n\nArguments\n\nata_model::CCMaximinModel : Required. The model built with start_ata() and with settings loaded by start_ata function.\n\n\n\n\n\n","category":"method"},{"location":"build/#Load-Design","page":"Build the ATA model","title":"Load Design","text":"","category":"section"},{"location":"build/","page":"Build the ATA model","title":"Build the ATA model","text":"A design can be loaded as a starting solution for the ATA or to print (or plot) the results.","category":"page"},{"location":"build/","page":"Build the ATA model","title":"Build the ATA model","text":"load_design!(design::Matrix{Any}, ata_model::AbstractModel)","category":"page"}]
}
