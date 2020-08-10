#cd("where your input files are")
#using ATAjl

ATAmodel = StartATA() #for resetting the ATA process (Needed)
println(LoadSettings!(ATAmodel; settings_file = "settingsATA.jl", bank_file = "bank.csv", bank_delim = ";")[2]) #add file with custom settings (Needed)
println(AddFriendSets!(ATAmodel)[2]) #add friend set variables (Optional)
println(AddEnemySets!(ATAmodel)[2]) #add enemy set variables (Optional)
println(AddConstr!(ATAmodel; constraints_file = "Constraints.csv", constraints_delim=";")[2]) #add categorical constraints (Optional)
println(AddOverlaps!(ATAmodel; overlap_file = "Overlap Matrix.csv", overlap_delim=";")[2]) #add overlap maxima (Optional)
println(AddExpScore!(ATAmodel)[2]) #add expected score constraints (Optional)
println(GroupByFriendSet!(ATAmodel)[2]) #add overlap maxima (Optional)
println(AddObjFun!(ATAmodel)[2]) #add objective function (Needed)
#Assemble!

#, optimizer_attributes = [("seconds",1000),("logLevel",3)], optimizer_constructor = Cbc.Optimizer ,
Assemble!(ATAmodel; solver="siman", start_temp=0.0001, geom_temp=0.1, max_time = 1000.0, n_item_sample = 40, n_test_sample = 36, verbosity = 2, conv_max=5, opt_feas = 0.9, n_fill = 1, feas_nh = 1, opt_nh = 200, optimizer_attributes = [("CPX_PARAM_TILIM",200)], optimizer_constructor = GLPK.Optimizer  )
#All the settings and outputs from optimization are in ATAmodel. See the struct in ATA.jl to understand how to retrieve all the information.
PrintResults(ATAmodel; GroupByFS = true, BSfolder = "BS") #BSfolder is for CCATA, save all the plots in RESULTS folder
