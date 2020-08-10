optimizer_constructor = CPLEX.Optimizer
optimizer_attributes =  [("CPX_PARAM_TILIM",500),("CPX_PARAM_PREIND",1)]
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

	#OVERLAP
	Pairs_t = combinations(ATAmodel.Settings.T)
	nPairs_t = size(Pairs_t,1)
	# ol_max = Array{Int64,1}(undef,nPairs_t)
	# fInd = [Pairs_t[pair][1] for pair in 1:nPairs_t]
	# fIndFirst = [Pairs_t[pair][2] for pair in 1:nPairs_t]
	# # for pair in 1:nPairs_t
	# 	ol_max[pair] = minimum([ATAmodel.Settings.olMax[fInd[pair],fIndFirst[pair]],minimum(ATAmodel.Settings.olMax[fInd[pair],:]),minimum(ATAmodel.Settings.olMax[fIndFirst[pair],:])])
	# end

	# #no overlap between x and y
	# Pairs_i = combinations(ATAmodel.Settings.nFS)
	# nPairs_i = size(Pairs_i,1)

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

	#overlap
	global c += ATAmodel.Settings.T


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


	#Overlap Vars new
	JuMP.@variable(m, x[i = 1:ATAmodel.Settings.nFS, t = 1:(ATAmodel.Settings.T)] <= Int(ATAmodel.Settings.forced0[t][i]), Bin)#, start=starting_design[i,t])
	JuMP.@variable(m, y[i = 1:ATAmodel.Settings.nFS, t = 1:(ATAmodel.Settings.T)] <= Int(ATAmodel.Settings.forced0[t][i]), Bin)#, start=zeros(ATAmodel.Settings.nFS,(ATAmodel.Settings.T-1))[i,t])
	#Constraints vars
	#JuMP.@variable(m, z[global c = 1:ncons] >= 0)
	global c = 1
	JuMP.@variable(m, z >= 0)

	#length
	for t=1:ATAmodel.Settings.T
		if ATAmodel.Constraints[t].length_max>0
			#if t==ATAmodel.Settings.T
				#JuMP.@constraint(m,  sum(x[i,t]*Int(ATAmodel.Settings.FS.Counts[i]) for i=1:ATAmodel.Settings.nFS)  <= ATAmodel.Constraints[t].length_max)
			#else
				JuMP.@constraint(m,  sum(x[i,t]*Int(ATAmodel.Settings.FS.Counts[i]) for i=1:ATAmodel.Settings.nFS) + sum(y[i,t]*Int(ATAmodel.Settings.FS.Counts[i]) for i=1:(ATAmodel.Settings.nFS)) <= ATAmodel.Constraints[t].length_max)
			#end
		end
		if ATAmodel.Constraints[t].length_min>0
			#if t==ATAmodel.Settings.T
			#	JuMP.@constraint(m,  sum(x[i,t]*Int(ATAmodel.Settings.FS.Counts[i]) for i=1:ATAmodel.Settings.nFS) >=  ATAmodel.Constraints[t].length_min)
			#else
				JuMP.@constraint(m,  sum(x[i,t]*Int(ATAmodel.Settings.FS.Counts[i]) for i=1:ATAmodel.Settings.nFS) + sum(y[i,t]*Int(ATAmodel.Settings.FS.Counts[i]) for i=1:(ATAmodel.Settings.nFS))>=  ATAmodel.Constraints[t].length_min )
			#end
			 global c+=1
		end
	end


	#Item Use
	for i in 1:ATAmodel.Settings.nFS
		if ATAmodel.IU.Min[i] .> 0
			JuMP.@constraint(m, sum(x[i,t] for t=1:ATAmodel.Settings.T) + sum(y[i,t] for t=1:(ATAmodel.Settings.T)) >=  ATAmodel.IU.Min[i] )
			 global c+=1
		end
		JuMP.@constraint(m, sum(x[i,t] for t=1:ATAmodel.Settings.T) + sum(y[i,t] for t=1:(ATAmodel.Settings.T)) <= ATAmodel.IU.Max[i])
		 global c+=1
	end

	#overlap new
	for t=1:ATAmodel.Settings.T
		#if t<ATAmodel.Settings.T
			JuMP.@constraint(m,  sum(y[i,t]*Int(ATAmodel.Settings.FS.Counts[i]) for i=1:ATAmodel.Settings.nFS) <= z)#Int(maximum(ATAmodel.Settings.olMax[t,:])))
			for i=1:ATAmodel.Settings.nFS
			 	JuMP.@constraint(m,  x[i,t] + y[i,t] <= 1)
			end
		#else
			# for i=1:ATAmodel.Settings.nFS
			# 	if Int(ATAmodel.Settings.forced0[t][i])<1
			#  		JuMP.upper_bound(x[i,t], Int(ATAmodel.Settings.forced0[t][i]))
			# 	end
			# end
		#end
	end
	#JuMP.@variable(m, v >= 0)
	for i=1:ATAmodel.Settings.nFS
		JuMP.@constraint(m,  sum(x[i,t] for t=1:(ATAmodel.Settings.T)) <=1)
	end
	for t=1:ATAmodel.Settings.T
		#JuMP.@constraint(m,  sum(y[i,t] for i=1:(ATAmodel.Settings.nFS))  >= 1)
		JuMP.@constraint(m,  sum(x[i,t] for i=1:(ATAmodel.Settings.nFS))  >= 1)
	end

	for pair in Pairs_t
		for i=1:ATAmodel.Settings.nFS
				#JuMP.@constraint(m,  x[i,pair[1]] + y[i,pair[2]] <= 1)
				#JuMP.@constraint(m,  y[i,pair[1]] + x[i,pair[2]] <= 1)
		end
	end

	#expected score
	for t=1:ATAmodel.Settings.T
		if size(ICF_new[t],1) > 0
			for k=1:size(ICF_new[t],1)
				# if t==ATAmodel.Settings.T
				# 	if ATAmodel.Constraints[t].ExS.Min[k] > 0
   				# 	 JuMP.@constraint(m, sum(x[i,t]*round(ATAmodel.Constraints[t].ExS.Val[k,i];digits = 3) for i=1:ATAmodel.Settings.nFS) >= ATAmodel.Constraints[t].ExS.Min[k]*ATAmodel.Constraints[t].length_min)
   				# 	 global c+=1
   				#  end
   				#  if ATAmodel.Constraints[t].ExS.Max[k] < 1
   				# 	 JuMP.@constraint(m, sum(x[i,t]*round(ATAmodel.Constraints[t].ExS.Val[k,i];digits = 3) for i=1:ATAmodel.Settings.nFS)  <= ATAmodel.Constraints[t].ExS.Max[k]*ATAmodel.Constraints[t].length_max)
   				# 	 global c+=1
   				#  end
				# else
				 if ATAmodel.Constraints[t].ExS.Min[k] > 0
					 JuMP.@constraint(m, sum(x[i,t]*round(ATAmodel.Constraints[t].ExS.Val[k,i];digits = 3) for i=1:ATAmodel.Settings.nFS) +sum(y[i,t]*round(ATAmodel.Constraints[t].ExS.Val[k,i];digits = 3) for i=1:ATAmodel.Settings.nFS) >= ATAmodel.Constraints[t].ExS.Min[k]*ATAmodel.Constraints[t].length_min)
					 global c+=1
				 end
				 if ATAmodel.Constraints[t].ExS.Max[k] < 1
					 JuMP.@constraint(m, sum(x[i,t]*round(ATAmodel.Constraints[t].ExS.Val[k,i];digits = 3) for i=1:ATAmodel.Settings.nFS) +sum(y[i,t]*round(ATAmodel.Constraints[t].ExS.Val[k,i];digits = 3) for i=1:ATAmodel.Settings.nFS) <= ATAmodel.Constraints[t].ExS.Max[k]*ATAmodel.Constraints[t].length_max)
					 global c+=1
				 end
			 #end
			end
		end
	end

	# #Constraints
	for t=1:ATAmodel.Settings.T
		if size(ATAmodel.Constraints[t].catConstrA,1)>0
			# if t==ATAmodel.Settings.T
			# 	for constr=1:size(ATAmodel.Constraints[t].catConstrA,1)
			# 		JuMP.@constraint(m, sum(x[i,t]*ATAmodel.Constraints[t].catConstrA[constr,i] for i=1:ATAmodel.Settings.nFS) <= ATAmodel.Constraints[t].catConstrb[constr])
			# 		 global c+=1
			# 	 end
			# else
				for constr=1:size(ATAmodel.Constraints[t].catConstrA,1)
					JuMP.@constraint(m, sum(x[i,t]*ATAmodel.Constraints[t].catConstrA[constr,i] for i=1:ATAmodel.Settings.nFS) + sum(y[i,t]*ATAmodel.Constraints[t].catConstrA[constr,i] for i=1:ATAmodel.Settings.nFS) <= ATAmodel.Constraints[t].catConstrb[constr])
					 global c+=1
				end
			# end
		end
	end


	ncons=copy(c)-1

	if ATAmodel.Settings.OptType == "MAXIMIN"
		#Objective bound
		JuMP.@variable(m, 0 <= w <= 40)
		for t = 1:ATAmodel.Settings.T
			# if t==ATAmodel.Settings.T
			# 	for k = 1:size(ATAmodel.Obj.OptPts[t],1)
			# 		JuMP.@constraint(m, sum(round(IIF[t][k,i];digits = 4)*x[i,t] for i = 1:ATAmodel.Settings.nFS)  >= w)
			# 	end
			# else
				for k = 1:size(ATAmodel.Obj.OptPts[t],1)
					JuMP.@constraint(m, sum(round(IIF[t][k,i];digits = 4)*x[i,t] for i = 1:ATAmodel.Settings.nFS) + sum(round(IIF[t][k,i];digits = 4)*y[i,t] for i = 1:ATAmodel.Settings.nFS) >= w)
				end
			# end
		end
		JuMP.@objective(m, Max, w-z)# + (0.9 * sum(0.0 for global c = 1:ncons)))
	elseif ATAmodel.Settings.OptType == ""
		#JuMP.@objective(m, Min, (sum(0.0 for c=1:ncons)))
	end

	JuMP.optimize!(m)
	message *= string("The model has termination status:", JuMP.termination_status(m))
	println(string("The model has termination status:", JuMP.termination_status(m)))
	design = abs.(JuMP.value.(x)) +  abs.(JuMP.value.(y))
	writedlm(string(results_folder,"/design.csv"),design)
	ATAmodel.Output.Design = design
	JLD2.@save string(results_folder,"/ATAmodel.jld2") ATAmodel
	return message
