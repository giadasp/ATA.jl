function jumpATA!(ATAmodel::model; starting_design=Matrix{Float64}(undef,0,0), optimizer_constructor = GLPK.Optimizer, optimizer_attributes = [("tm_lim",500000),("msg_lev",3)], results_folder= "RESULTS")
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
		message *= "Assembling tests with NO objective function..."
	end

	#OVERLAP
	opMatrix = ATAmodel.Settings.olMax
	Pairs = combinations(ATAmodel.Settings.T)
	nPairs = size(Pairs,1)
	ol_max = Array{Int64,1}(undef,nPairs)
	fInd = [Pairs[pair][1] for pair in 1:nPairs]
	fIndFirst = [Pairs[pair][2] for pair in 1:nPairs]
	for pair in 1:nPairs
		ol_max[pair] = opMatrix[fInd[pair],fIndFirst[pair]]
	end

################################################################################
#                                count constraints
################################################################################
	c = 1
	#length
	for t=1:ATAmodel.Settings.T
		c += ATAmodel.Constraints[t].length_max .> 0 ? 1 : 0
		c += ATAmodel.Constraints[t].length_min .> 0 ? 1 : 0
	end

	# item use
	c += size(ATAmodel.IU.Max,1)
	c += sum(ATAmodel.IU.Min .> 0)
	#overlap classic
	c += nPairs

	#expected score
	for t=1:ATAmodel.Settings.T
		if size(ATAmodel.Constraints[t].ExS.Val,1) > 0
			 c += 2
		end
	end

	#Friend Sets
	if size(ATAmodel.Settings.FS.Items,1)!=nItems
		for t=1:ATAmodel.Settings.T
			for fs=1:size(ATAmodel.Settings.FS.Items,1)
				c += 1
			end
		end
	end

	#enemy sets
	if size(ATAmodel.Settings.ES.Sets,1)>0
		c += ATAmodel.Settings.T * size(ATAmodel.Settings.ES.Sets,1)
	end

	#Constraints
	for t=1:ATAmodel.Settings.T
		if size(ATAmodel.Constraints[t].catConstrA,1)>0
			 c += size(ATAmodel.Constraints[t].catConstrA,1)
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
		if (size(starting_design,1)!=nItems || size(starting_design,2)!=ATAmodel.Settings.T)
			message *= "- Starting design must be of size: (nItems x T).\n"
			return message
		end
		if (any(starting_design!=0 && starting_design!=1))
			message *= "- Starting design must contain only 1 or 0.\n"
			return message
		end
	else
		starting_design = zeros(nItems,ATAmodel.Settings.T)
	end

	#decision variables
	JuMP.@variables(m, begin
           x[i = 1:nItems , t = 1:ATAmodel.Settings.T], Bin, (start = starting_design[i, t])
       end)

	#JuMP.@variable(m, x[i = 1:nItems,t = 1:ATAmodel.Settings.T] , Bin)

	#Overlap Vars
	JuMP.@variable(m, y[i = 1:nItems, p = 1:nPairs], Bin)

	#Constraints vars
	JuMP.@variable(m, z[c = 1:ncons] >= 0)
	c = 1
	#length
	for t=1:ATAmodel.Settings.T
		if ATAmodel.Constraints[t].length_max>0
			JuMP.@constraint(m,  sum(x[i,t] for i=1:nItems) - ATAmodel.Constraints[t].length_max <= z[c])
			c+=1
		end
		if ATAmodel.Constraints[t].length_min>0
			JuMP.@constraint(m,  -sum(x[i,t] for i=1:nItems) + ATAmodel.Constraints[t].length_min<=z[c])
			 c+=1
		end
	end

	# Item Use
	for i in 1:nItems
		if ATAmodel.IU.Min[i] .> 0
			JuMP.@constraint(m,  ATAmodel.IU.Min[i]-sum(x[i,t] for t=1:ATAmodel.Settings.T) <= z[c] )
			 c+=1
		end
		JuMP.@constraint(m,  sum(x[i,t] for t=1:ATAmodel.Settings.T) - ATAmodel.IU.Max[i] <= z[c])
		 c+=1
	end

	#overlap classic
	for p=1:nPairs
		JuMP.@constraint(m, sum(y[i,p] for i=1:nItems)-ol_max[p]<=z[c])
		 c+=1
		JuMP.@constraint(m,  [i = 1:nItems], 2*y[i,p]<=x[i,fInd[p]]+x[i,fIndFirst[p]])
		JuMP.@constraint(m,  [i = 1:nItems], y[i,p]>=x[i,fInd[p]]+x[i,fIndFirst[p]]-1)
	end

	#expected score
	for t=1:ATAmodel.Settings.T
		if size(ATAmodel.Constraints[t].ExS.Val,1)>0
			JuMP.@constraint(m, sum(x[i,t]*ATAmodel.Constraints[t].ExS.Val[i] for i=1:nItems)<=ATAmodel.Constraints[t].ExS.Max*ATAmodel.Constraints[t].length_min+z[c])
			 c+=1
			JuMP.@constraint(m, sum(x[i,t]*ATAmodel.Constraints[t].ExS.Val[i] for i=1:nItems)>=ATAmodel.Constraints[t].ExS.Min*ATAmodel.Constraints[t].length_max-z[c])
			 c+=1
		end
	end

	#Friend Sets
	if size(ATAmodel.Settings.FS.Items,1)!=nItems
		for t=1:ATAmodel.Settings.T
			for fs=1:size(ATAmodel.Settings.FS.Items,1)
				friend_set = ATAmodel.Settings.FS.Items[fs]
				length_fs = size(friend_set,1)
				JuMP.@constraint(m, sum(x[i,t] for i in friend_set[1:(length_fs-1)]) - (length_fs-1)*x[friend_set[length_fs],t] == z[c])
				c+=1
			end
		end
	end

	#Enemy Sets
	if size(ATAmodel.Settings.ES.Sets,1)>0
		for t=1:ATAmodel.Settings.T
			for es=1:size(ATAmodel.Settings.ES.Sets,1)
				JuMP.@constraint(m, sum(x[i,t] for i in ATAmodel.Settings.FS.Items[es])<=1+z[c])
				 c+=1
			end
		end
	end

	#Constraints
	for t=1:ATAmodel.Settings.T
		if size(ATAmodel.Constraints[t].catConstrA,1)>0
			for constr=1:size(ATAmodel.Constraints[t].catConstrA,1)
				JuMP.@constraint(m, sum(x[i,t]*ATAmodel.Constraints[t].catConstrA[constr,i] for i=1:nItems)<=ATAmodel.Constraints[t].catConstrb[constr]+z[c])
				 c+=1
			end
		end
	end

	ncons=copy(c)-1

	if ATAmodel.Settings.OptType == "MAXIMIN"
		#Objective bound
		JuMP.@variable(m, w >= 0)
		for t = 1:ATAmodel.Settings.T
			for k = 1:size(ATAmodel.Obj.OptPts[t],1)
				JuMP.@constraint(m, sum(round(IIF[t][k,i];digits = 4)*x[i,t] for i = 1:nItems) >= w)
			end
		end
		JuMP.@objective(m, Min, (- 0.1 * w) + (0.9 * sum(z[c] for c = 1:ncons)))
	elseif ATAmodel.Settings.OptType == ""
		JuMP.@objective(m, Min, (sum(z[c] for c=1:ncons)))
	end

	JuMP.optimize!(m)
	message *= string("The model has termination status:", JuMP.termination_status(m))
	println(string("The model has termination status:", JuMP.termination_status(m)))
	design = abs.(JuMP.value.(x))
	writedlm(string(results_folder,"/design.csv"),design)
	ATAmodel.Output.Design = design
	JLD2.@save string(results_folder,"/ATAmodel.jld2") ATAmodel
	return message
end
