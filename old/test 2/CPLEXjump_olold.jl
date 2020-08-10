optimizer_constructor = CPLEX.Optimizer
optimizer_attributes =  [("CPX_PARAM_TILIM",500)]
results_folder= "RESULTS"
starting_design=Matrix{Float64}(undef,0,0)
	message = ""
	if !(results_folder in readdir())
		mkdir(results_folder)
	else
		message *= string("You have already a folder with this name, files in ", results_folder," will be overwritten.\n")
	end

	nItems = ATAmodel.Settings.nItems
	if ATAmodel.Settings.OptType == "MAXIMIN"
		if isfile("OPT/IIF.jld2")
			JLD2.@load "OPT/IIF.jld2" IIF
			message *= "- Assembling tests with MAXIMIN..."
		else
			message *= "No IIF.jld2 file in OPT folder, Run AddObjFun!() first!"
			return message
		end
	elseif ATAmodel.Settings.OptType == "CC"
		message *= "You must use the Simulated Annealing algorithm to assemble tests with CC objective function."
		return message
	elseif ATAmodel.Settings.OptType == ""
		IIF = []
		message *= "Assembling tests with NO objective function..."
	end

	if ATAmodel.Settings.nFS==0
		ATAmodel.Settings.nFS = ATAmodel.Settings.nItems
	end

	#OVERLAP new
	Pairs_t = combinations(ATAmodel.Settings.T)
	nPairs_t = size(Pairs_t,1)

	#overlap old
	opMatrix = ATAmodel.Settings.olMax
	Pairs = combinations(ATAmodel.Settings.T)
	nPairs = size(Pairs,1)
	ol_max = Array{Int64,1}(undef,nPairs)
	fInd = [Pairs[pair][1] for pair in 1:nPairs]
	fIndFirst = [Pairs[pair][2] for pair in 1:nPairs]
	for pair in 1:nPairs
		ol_max[pair] = opMatrix[fInd[pair],fIndFirst[pair]]
	end


#Freind sets groups
if ATAmodel.Settings.nFS!=ATAmodel.Settings.nItems
	#group IIFs
	IIF_new = Vector{Matrix{Float64}}(undef,ATAmodel.Settings.T)
	ICF_new = Vector{Matrix{Float64}}(undef,ATAmodel.Settings.T)
	for t = 1:ATAmodel.Settings.T
		if size(IIF,1)>0
			IIF_new[t] = Matrix{Float64}(undef, size(IIF[t],2), ATAmodel.Settings.nFS)
		end
		if size(ATAmodel.Constraints[t].ExS.Val,1) > 0
			ICF_new[t] = Matrix{Float64}(undef, size(ATAmodel.Constraints[t].ExS.Val,2), ATAmodel.Settings.nFS)
		end
	end
	for fs = 1:ATAmodel.Settings.nFS
		for t = 1:ATAmodel.Settings.T
			if size(IIF,1)>0
				IIF_new[t] = sum(IIF[t][:, ATAmodel.Settings.FS.Items[fs]], dims = 2)
			end
			if size(ATAmodel.Constraints[t].ExS.Val,1) > 0
				ICF_new[t] = sum(ATAmodel.Constraints[t].ExS.Val[:, ATAmodel.Settings.FS.Items[fs]], dims = 2)
			end
		end
	end
end

################################################################################
#                                count constraints
################################################################################
	global c = 1
	#length
	for t=1:ATAmodel.Settings.T
		global c += ATAmodel.Constraints[t].length_max .> 0 ? 1 : 0
		global c += ATAmodel.Constraints[t].length_min .> 0 ? 1 : 0
	end

	# item use
	global c += size(ATAmodel.IU.Max, 1)
	global c += sum(ATAmodel.IU.Min .> 0)

	#overlap new
	#global c += ATAmodel.Settings.T
	#overlap old
	global c += nPairs

	#expected score
	for t=1:ATAmodel.Settings.T
		if size(ICF_new[t],1) > 0
			for k=1:size(ICF_new[t],1)
				 if ATAmodel.Constraints[t].ExS.Min[k] > 0
					 global c+=1
				 end
				 if ATAmodel.Constraints[t].ExS.Max[k] < 1
					 global c+=1
				 end
		 	end
		end
	end

	#Constraints
	for t=1:ATAmodel.Settings.T
		if size(ATAmodel.Constraints[t].catConstrA,1)>0
			 global c += size(ATAmodel.Constraints[t].catConstrA,1)
		end
	end
	ncons=copy(c)-1
################################################################################
#                                JuMP model
################################################################################

	m = JuMP.Model()
	JuMP.set_optimizer(m, optimizer_constructor)
	for (name, value) in optimizer_attributes
		JuMP.set_optimizer_attribute(m, name, value)
	end

	#starting design check
	if size(starting_design,1)>0
		if ATAmodel.Settings.nFS != ATAmodel.Settings.nItems
			if (size(starting_design,1) != ATAmodel.Settings.nFS || size(starting_design,2) != ATAmodel.Settings.T)
				message *= "- Starting design must be of size: (nFS x T).\n"
				return message
			end
		else
			if (size(starting_design,1) != ATAmodel.Settings.nItems || size(starting_design,2) != ATAmodel.Settings.T)
				message *= "- Starting design must be of size: (nItems x T).\n"
				return message
			end
		end
		if (any(starting_design != 0 && starting_design != 1))
			message *= "- Starting design must contain only 1 or 0.\n"
			return message
		end
	else
		starting_design = zeros(ATAmodel.Settings.nFS,ATAmodel.Settings.T)
	end

	#decision variables
	JuMP.@variables(m, begin
       x[i = 1:ATAmodel.Settings.nFS , t = 1:ATAmodel.Settings.T] <= Int(ATAmodel.Settings.forced0[t][i]), Bin, (start = starting_design[i, t])
    end)

	#Overlap Vars
	JuMP.@variable(m, y[i = 1:ATAmodel.Settings.nFS, p = 1:nPairs], Bin)

	#Constraints vars
	#JuMP.@variable(m, z[global c = 1:ncons] >= 0)
	global c = 1

	#length
	for t=1:ATAmodel.Settings.T
		if ATAmodel.Constraints[t].length_max>0
			JuMP.@constraint(m,  sum(x[i,t]*ATAmodel.Settings.FS.Counts[i] for i=1:ATAmodel.Settings.nFS)  <= ATAmodel.Constraints[t].length_max)
			global c+=1
		end
		if ATAmodel.Constraints[t].length_min>0
			JuMP.@constraint(m,  sum(x[i,t]*ATAmodel.Settings.FS.Counts[i] for i=1:ATAmodel.Settings.nFS) >= ATAmodel.Constraints[t].length_min)
			 global c+=1
		end
	end

	#overlap classic
	for p=1:nPairs
		JuMP.@constraint(m, sum(y[i,p]*ATAmodel.Settings.FS.Counts[i] for i=1:ATAmodel.Settings.nFS) <=ol_max[p])
		 global c+=1
		JuMP.@constraint(m,  [i = 1:ATAmodel.Settings.nFS], 2*y[i,p] <= x[i,fInd[p]]+x[i,fIndFirst[p]])
		JuMP.@constraint(m,  [i = 1:ATAmodel.Settings.nFS], y[i,p] >= x[i,fInd[p]]+x[i,fIndFirst[p]]-1)
	end

	# Item Use
	for i in 1:ATAmodel.Settings.nFS
		if ATAmodel.IU.Min[i] .> 0
			JuMP.@constraint(m,   sum(x[i,t] for t=1:ATAmodel.Settings.T) >= ATAmodel.IU.Min[i] )
			 global c+=1
		end
		JuMP.@constraint(m,  sum(x[i,t] for t=1:ATAmodel.Settings.T)  <= ATAmodel.IU.Max[i])
		 global c+=1
	end

	#expected score
	for t=1:ATAmodel.Settings.T
		if size(ICF_new[t],1) > 0
			for k=1:size(ICF_new[t],1)
				 if ATAmodel.Constraints[t].ExS.Min[k] > 0
					 JuMP.@constraint(m, sum(x[i,t]*ATAmodel.Constraints[t].ExS.Val[k,i] for i=1:ATAmodel.Settings.nFS) >= ATAmodel.Constraints[t].ExS.Min[k]*ATAmodel.Constraints[t].length_min)
					 global c+=1
				 end
				 if ATAmodel.Constraints[t].ExS.Max[k] < 1
					 JuMP.@constraint(m, sum(x[i,t]*ATAmodel.Constraints[t].ExS.Val[k,i] for i=1:ATAmodel.Settings.nFS) <= ATAmodel.Constraints[t].ExS.Max[k]*ATAmodel.Constraints[t].length_max)
					 global c+=1
				 end
			end
		end
	end

	# #Constraints
	for t in 1:ATAmodel.Settings.T
		if size(ATAmodel.Constraints[t].catConstrA,1)>0
			for constr=1:(size(ATAmodel.Constraints[t].catConstrA,1))
				JuMP.@constraint(m, sum(x[i,t]*ATAmodel.Constraints[t].catConstrA[constr,i] for i=1:ATAmodel.Settings.nFS) <= ATAmodel.Constraints[t].catConstrb[constr])
				global c+=1
			end
		end
	end


	ncons=copy(c)-1

	if ATAmodel.Settings.OptType == "MAXIMIN"
		#Objective bound
		JuMP.@variable(m, w >= 0)
		for t = 1:ATAmodel.Settings.T
			for k = 1:size(ATAmodel.Obj.OptPts[t],1)
				JuMP.@constraint(m, sum(round(IIF[t][k,i];digits = 4)*(x[i,t]) for i = 1:ATAmodel.Settings.nFS) >= w)
			end
		end
		JuMP.@objective(m, Min, (- 0.1 * w))# + (0.9 * sum(0.0 for global c = 1:ncons)))
	elseif ATAmodel.Settings.OptType == ""
		#JuMP.@objective(m, Min, (sum(0.0 for c=1:ncons)))
	end

	JuMP.optimize!(m)
	message *= string("The model has termination status:", JuMP.termination_status(m))
	println(string("The model has termination status:", JuMP.termination_status(m)))
	design = abs.(JuMP.value.(x))
	writedlm(string(results_folder,"/design.csv"),design)
	ATAmodel.Output.Design = design
	JLD2.@save string(results_folder,"/ATAmodel.jld2") ATAmodel
	return message
