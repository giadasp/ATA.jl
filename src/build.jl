#cc metho

function StartATA()
	if "OPT" in readdir()
		files = readdir("OPT")
		if size(files, 1)>0
			for f in 1:size(files, 1)
				rm(string("OPT/", files[f]), force = true)
			end
		end
	else
		mkdir("OPT")
	end
	#update model
	ATAmodel = model()
	JLD2.@save "OPT/ATAmodel.jld2" ATAmodel
	return ATAmodel
end

function LoadSettings!(ATAmodel::model; settings_file = "settingsATA.jl", bank_file = "bank.csv", bank_delim = ";")
	message = ["", ""]
	if isfile(bank_file)
		ATAmodel.Settings.Bank = CSV.read(bank_file, delim = bank_delim)
		message[2] = message[2] * "- Item bank file read.\n"
	else
		return ["danger", "Not a valid file for Bank"]
	end
	if !isfile(settings_file)
		return ["danger", "Settings file not valid"]
	else
		include(settings_file)
		ATAmodel.Settings.nItems = Inputs.nItems
		ATAmodel.Settings.T = Int(sum(Inputs.T))
		ATAmodel.Settings.Tg = Inputs.T
		#initialize cosntraints
		ATAmodel.Constraints = [Constraint() for t = 1:ATAmodel.Settings.T]
		ATAmodel.Settings.nGroups = size(Inputs.Groups, 1)
		x_forced0 = Vector{Vector{Bool}}(undef, ATAmodel.Settings.T)
		for t = 1:ATAmodel.Settings.T
			x_forced0[t] = fill(true, ATAmodel.Settings.nItems)
			ATAmodel.Constraints[t].catConstrA = zeros(Float64, 0, ATAmodel.Settings.nItems)
		end

		open("OPT/Settings.jl", "w") do f
			write(f, "#Settings \n\n")

				ATAmodel.Settings.IRT.model = Inputs.IRTmodel
				ATAmodel.Settings.IRT.parameters = DataFrame(ATAmodel.Settings.Bank[!, Symbol.(Inputs.IRTparameters)])
				ATAmodel.Settings.IRT.parametrization = Inputs.IRTparametrization
				ATAmodel.Settings.IRT.D = Inputs.IRTD

				if ATAmodel.Settings.IRT.model == "1PL"
					DataFrames.rename!(ATAmodel.Settings.IRT.parameters, [:b])#nqp values in interval\r\n",
				elseif ATAmodel.Settings.IRT.model == "2PL"
					DataFrames.rename!(ATAmodel.Settings.IRT.parameters, [:a, :b]) #nqp values in interval\r\n",
				elseif ATAmodel.Settings.IRT.model == "3PL"
					DataFrames.rename!(ATAmodel.Settings.IRT.parameters, [:a, :b, :c]) #nqp values in interval\r\n",
				else
					return ["danger","Only 1PL, 2PL and 3PL IRT models are allowed."]
				end
				CSV.write("OPT/IRTparameters.csv", ATAmodel.Settings.IRT.parameters)
				message[2] = message[2] * "- IRT item parameters loaded.\n"
			if Inputs.EnemySetsVar != String[]
				ATAmodel.Settings.ES.Var = Symbol.(Inputs.EnemySetsVar)
				val = Symbol.(Inputs.EnemySetsVar)
				write(f, "EnemySetsVar = $val\n\n")
				message[2] = message[2] * "- Variable for Enemy Sets loaded.\n"
			end
			if Inputs.FriendSetsVar != String[]
				ATAmodel.Settings.FS.Var = Symbol.(Inputs.FriendSetsVar)
				val = Symbol.(Inputs.FriendSetsVar)
				write(f, "ATAmodel.Settings.FS.Var = $val\n\n")
				message[2] = message[2] * "- Variable for Friend Sets loaded.\n"
			end

			if Inputs.length_min != Int64[]
				lengthmin = zeros(Int64, ATAmodel.Settings.T)
				lengthweight = ones(Int64, ATAmodel.Settings.T)
				t1 = 1
				for g = 1:ATAmodel.Settings.nGroups
					for t = 1:ATAmodel.Settings.Tg[g]
						lengthmin[t1] = Int(Inputs.length_min[g])
						lengthweight[t1] = Inputs.length_weight[g]
						t1+= 1
					end
				end
				for t = 1:ATAmodel.Settings.T
					ATAmodel.Constraints[t].catConstrA = vcat(ATAmodel.Constraints[t].catConstrA, (-lengthweight[t]).*ones(Float64, ATAmodel.Settings.nItems)')
					ATAmodel.Constraints[t].catConstrb = vcat(ATAmodel.Constraints[t].catConstrb, -lengthmin[t]*lengthweight[t])
					ATAmodel.Constraints[t].length_min = lengthmin[t]
				end
				write(f, "length_min = $lengthmin\n")
				message[2] = message[2] * "- Minimum length of tests constrained.\n"
			end

			if Inputs.length_max != Int64[]
				lengthmax = zeros(Int64, ATAmodel.Settings.T)
				lengthweight = ones(Int64, ATAmodel.Settings.T)
				t1 = 1
				for g = 1:ATAmodel.Settings.nGroups
					for t = 1:ATAmodel.Settings.Tg[g]
						lengthmax[t1] = Int(Inputs.length_max[g])
						lengthweight[t1] = Inputs.length_weight[g]
						t1+= 1
					end
				end
				for t = 1:ATAmodel.Settings.T
					ATAmodel.Constraints[t].catConstrA = vcat(ATAmodel.Constraints[t].catConstrA, (lengthweight[t]).*ones(ATAmodel.Settings.nItems)')
					ATAmodel.Constraints[t].catConstrb = vcat(ATAmodel.Constraints[t].catConstrb, lengthmax[t]*lengthweight[t])
					ATAmodel.Constraints[t].length_max = lengthmax[t]
				end
				write(f, "length_max = $lengthmax\n")
				message[2] = message[2] * "- Maximum length of tests constrained.\n"
			end

			if Inputs.ExSVar != String[]
				t1 = 1
				for g = 1:ATAmodel.Settings.nGroups
					for t = 1:ATAmodel.Settings.Tg[g]
						ATAmodel.Constraints[t1].ExS.Var = Symbol.(Inputs.ExSVar[g])
						t1+= 1
					end
				end
			end
			t1 = 1
			if Inputs.ExSPts != Float64[]
				for g = 1:ATAmodel.Settings.nGroups
					for t = 1:ATAmodel.Settings.Tg[g]
						ATAmodel.Constraints[t1].ExS.Pts = Inputs.ExSPts[g]
						t1+= 1
					end
				end
			end
			message[2] = message[2] * "- Expected score variable and points loaded.\n"

			if Inputs.ExS_min != Float64[]
				t1 = 1
				for g = 1:ATAmodel.Settings.nGroups
					for t = 1:ATAmodel.Settings.Tg[g]
						ATAmodel.Constraints[t1].ExS.Min = Inputs.ExS_min[g]
						t1+= 1
					end
				end
				message[2] = message[2] * "- Minimum expected score constrained.\n"
			end
			if Inputs.ExS_max != Float64[]
				t1 = 1
				for g = 1:ATAmodel.Settings.nGroups
					for t = 1:ATAmodel.Settings.Tg[g]
						ATAmodel.Constraints[t1].ExS.Max = Inputs.ExS_max[g]
						t1+= 1
					end
				end
				message[2] = message[2] * "- Maximum expected score constrained.\n"
			end
			if Inputs.sumVars != Vector{Vector{String}}(undef,0)
				t1 = 1
				for g = 1:ATAmodel.Settings.nGroups
					if !ismissing(sumVars_min[g])
						for v = 1:length(Inputs.sumVars[g])
							var = ATAmodel.Settings.Bank[Symbol(Inputs.sumVars[g][v])]
							var[ismissing.(var)].= zero(Float64)
							#min
							if !ismissing(sumVars_min[g][v])
								for t = 1:ATAmodel.Settings.Tg[g]
									if g>1
										t1 = (g-1)*Tg[g-1]+t
									else
										t1 = copy(t)
									end
									ATAmodel.Constraints[t1].catConstrA = vcat(ATAmodel.Constraints[t1].catConstrA, .-var')
									ATAmodel.Constraints[t1].catConstrb = vcat(ATAmodel.Constraints[t1].catConstrb, -sumVars_min[g])
								end
							end
							#min
							if !ismissing(sumVars_max[g][v])
								for t = 1:ATAmodel.Settings.Tg[g]
									if g>1
										t1 = (g-1)*Tg[g-1]+t
									else
										t1 = copy(t)
									end
									ATAmodel.Constraints[t1].catConstrA = vcat(ATAmodel.Constraints[t1].catConstrA, var')
									ATAmodel.Constraints[t1].catConstrb = vcat(ATAmodel.Constraints[t1].catConstrb, sumVars_max[g])
								end
							end
						end
					end
				end
				message[2] = message[2] * string("- Sum of variables", Inputs.sumVars , " will be constrained.\n")
			end
			if size(Inputs.ItemUse_min, 1) > 0
				ATAmodel.IU.Min = Inputs.ItemUse_min
				message[2] = message[2] * "- Minimum item use constrained.\n"
			end
			if size(Inputs.ItemUse_max, 1) > 0
				ATAmodel.IU.Max = Inputs.ItemUse_max
				for v = 1:ATAmodel.Settings.T
					x_forced0[v][findall(ATAmodel.IU.Max.<1)] .= false
				end
				message[2] = message[2] * "- Maximum item use constrained.\n"
			end
			if Inputs.OptType != ""
				ATAmodel.Settings.OptType = Inputs.OptType
				message[2] = message[2] * "- Optimization type loaded.\n"
			end
			if size(Inputs.OptPts, 1) > 0
				ATAmodel.Obj.OptPts = Vector{Vector{Float64}}(undef, ATAmodel.Settings.T)
				t1 = 1
				for g = 1:ATAmodel.Settings.nGroups
					for t = 1:ATAmodel.Settings.Tg[g]
						ATAmodel.Obj.OptPts[t1] = Inputs.OptPts[g]
						t1+= 1
					end
				end
				message[2] = message[2] *  "- Points to optimize loaded.\n"
			end
			ATAmodel.Obj.AuxInt = Inputs.AuxInt
			ATAmodel.Obj.AuxFloat = Inputs.AuxFloat
			message[2] = message[2] * "- Auxiliars vars for optimization loaded.\n"
			#fictiuos friendSets
			ATAmodel.Settings.FS.Counts = ones(ATAmodel.Settings.nItems)
			ATAmodel.Settings.FS.Sets = string.(collect(1:ATAmodel.Settings.nItems))
			ATAmodel.Settings.FS.Items = [[i] for i = 1:ATAmodel.Settings.nItems]
			ATAmodel.Settings.forced0 = x_forced0
			if Inputs.CATEGORIES != String[]
				val = Symbol.(Inputs.CATEGORIES)
				ATAmodel.Output.Categories = copy(val)
				write(f, "CATEGORIES = $val\n\n")
				message[2] = message[2] * "- CATEGORIES for output loaded.\n"
			end
		end
		message[2] = message[2] * string("Assemble ", ATAmodel.Settings.T, " forms divided in ", ATAmodel.Settings.nGroups, " groups.\n")
		#update model
		JLD2.@save "OPT/ATAmodel.jld2" ATAmodel

	end
	message[1] = "success"
	return message
end


function AddFriendSets!(ATAmodel::model)
	message = ["", ""]
	if !isfile("OPT/Settings.jl")
		return  ["danger","Run LoadSettings!(model) before!"]
	else
		nItems = ATAmodel.Settings.nItems
		FriendSets = string.(unique(vcat([unique(skipmissing(ATAmodel.Settings.Bank[!, (ATAmodel.Settings.FS.Var[isv])])) for isv = 1:(size(ATAmodel.Settings.FS.Var, 1))]...)))
		nFS = size(FriendSets, 1)
		FScounts = zeros(Int, nFS)
		FSItems = [zeros(Int, 0) for i = 1 : nFS]
		single_items = Vector{Union{Missing,String}}([missing for i = 1 : nItems])
		for i = 1 : nItems
			units = [ATAmodel.Settings.Bank[i, ATAmodel.Settings.FS.Var[isv]] for isv = 1:(size(ATAmodel.Settings.FS.Var, 1))]
			if all(ismissing.(units))
				single_items[i] = string(i)
			else
				nonmissing = units[.!ismissing.(units)]
				fs = Int64[]
				f_i = 1
				for f in FriendSets
					if f in nonmissing
						push!(fs, f_i)
					end
					f_i += 1
				end
				for f in fs
					FSItems[f] = vcat(FSItems[f], i)
					FScounts[f] += 1
				end
			end
		end
		push!(ATAmodel.Settings.FS.Var, :SINGLE_FS)
		DataFrames.insertcols!(ATAmodel.Settings.Bank, size(ATAmodel.Settings.Bank,2), :SINGLE_FS => single_items)
		items_single = findall(.!ismissing.(ATAmodel.Settings.Bank[!,:SINGLE_FS]))
		FSItems = vcat(FSItems, [[i] for i in items_single])
		FriendSets = vcat(FriendSets, ATAmodel.Settings.Bank[items_single,:SINGLE_FS])
		FScounts = vcat(FScounts, ones(Int64,size(items_single,1)))
		nFS = size(FScounts,1)
		writedlm("OPT/FriendSets.csv", FriendSets)

		#update model
		ATAmodel.Settings.nFS = nFS
		ATAmodel.Settings.FS.Sets = FriendSets
		ATAmodel.Settings.FS.Items = FSItems
		ATAmodel.Settings.FS.Counts = FScounts
		JLD2.@save "OPT/ATAmodel.jld2" model
		return ["success", string("- ",nFS, " friend sets added.")]
	end
end

function AddEnemySets!(ATAmodel::model)
	Bank = ATAmodel.Settings.Bank
	nItems = ATAmodel.Settings.nItems
	if size(ATAmodel.Constraints[1].catConstrb, 1)>0
		A = Vector{Matrix{Float64}}(undef, ATAmodel.Settings.T)
		b = Vector{Vector{Float64}}(undef, ATAmodel.Settings.T)
		for t = 1:ATAmodel.Settings.T
			A[t] = ATAmodel.Constraints[t].catConstrA
			b[t] = ATAmodel.Constraints[t].catConstrb
		end

	else
		A = Vector{Matrix{Float64}}(undef, ATAmodel.Settings.T)
		b = Vector{Vector{Float64}}(undef, ATAmodel.Settings.T)
		for t = 1:ATAmodel.Settings.T
			A[t] = Matrix{Float64}(undef, 0, ATAmodel.Settings.nItems)
			b[t] = Vector{Float64}(undef, 0)
		end
	end

	if size(ATAmodel.Settings.ES.Var, 1)>0
		EnemySetsVar = ATAmodel.Settings.ES.Var
		EnemySets = unique(vcat([unique(skipmissing(Bank[!, (EnemySetsVar[isv])])) for isv = 1:(size(EnemySetsVar, 1))]...))
		writedlm("EnemySets.csv", EnemySets)
		nES = size(EnemySets, 1)
		Sets = Vector{Vector{Int64}}(undef, nES)
		for es = 1:nES
			Sets[es] = Int64[]
			for t = 1:ATAmodel.Settings.T
				A[t] = vcat(A[t], zeros(ATAmodel.Settings.nItems)')
			end
			for i = 1:nItems
				if any(skipmissing(vcat([((Bank[!, (EnemySetsVar[isv])][i] == EnemySets[es]) == true) for isv = 1:(size(EnemySetsVar, 1))]...)))
					#ESMaATAmodel.Settings.Tg[nItems, es] = 1
					for t = 1:ATAmodel.Settings.T
						A[t][end, i] = 1
					end
					Sets[es] = vcat(Sets[es], i)
				end
			end
			for t = 1:ATAmodel.Settings.T
				push!(b[t], 1)
			end
		end
		for t = 1:ATAmodel.Settings.T
			writedlm("OPT/A_$t.csv", A[t])
			writedlm("OPT/b_$t.csv", b[t])
		end
		#update model
		ATAmodel.Settings.ES.Names = EnemySets
		ATAmodel.Settings.ES.Sets = Sets
		for t = 1:ATAmodel.Settings.T
			ATAmodel.Constraints[t].catConstrb = b[t]
			ATAmodel.Constraints[t].catConstrA = A[t]
		end
		JLD2.@save "OPT/ATAmodel.jld2" model
		return ["success", string("- ",nES, " enemy sets added. ")]
	end
end

function AddConstr!(ATAmodel::model; constraints_file = "CategoricalConstraints.csv", constraints_delim = ";")
	message=["", ""]
	nItems = ATAmodel.Settings.nItems
	Bank = ATAmodel.Settings.Bank
	A = Vector{Matrix{Float64}}(undef, ATAmodel.Settings.T)
	b = Vector{Vector{Float64}}(undef, ATAmodel.Settings.T)
	if size(ATAmodel.Constraints[1].catConstrb, 1)>0
		for t = 1:ATAmodel.Settings.T
			A[t] = copy(ATAmodel.Constraints[t].catConstrA)
			b[t] = copy(ATAmodel.Constraints[t].catConstrb)
		end
	else
		for t = 1:ATAmodel.Settings.T
			A[t] = zeros(Float64, 0, ATAmodel.Settings.nItems)
			b[t] = Vector{Float64}(undef, 0)
		end

	end
	if !isfile(constraints_file)
		return ["danger",string(constraints_file," doesn't exist.")]
	else
		x_forced0 = ATAmodel.Settings.forced0
		Categoricalconsts = CSV.read(constraints_file, delim = constraints_delim)
		if size(Categoricalconsts, 1)==0
			message[2] = message[2] * string("No lines in file ", constraints_file," nothing added.\n")
		else
			Categoricalconsts.var = Symbol.(Categoricalconsts.var)
			CatCons = Vector{Vector{Float64}}(undef, ATAmodel.Settings.T)
			for g = 1:ATAmodel.Settings.nGroups
				if g>1
					tests = collect((sum(ATAmodel.Settings.Tg[1:(g-1)])+1):(sum(ATAmodel.Settings.Tg[1:(g-1)])+ATAmodel.Settings.Tg[g]))
				else
					tests = collect(1:ATAmodel.Settings.Tg[1])
				end
				GroupCatCons = findall((Categoricalconsts.group.== g))
				for con in GroupCatCons
					for t in tests
						if !ismissing(Categoricalconsts[con, :min])
							indices = string.(Bank[!, Symbol.(Categoricalconsts[con, :var])]) .== Categoricalconsts[con, :value]
							indices[ismissing.(indices)] .= false
							indices = findall(indices .== true)
							A[t] = vcat(A[t], zeros(Float64, 1, nItems))
							if ismissing(Categoricalconsts[con, :weight])
								A[t][end, indices] .= -1.0
								push!(b[t], -Categoricalconsts[con, :min])
							else
								A[t][end, indices] .= -(Categoricalconsts[con, :weight])
								push!(b[t], -Categoricalconsts[con, :weight]*Categoricalconsts[con, :min])
							end
						end
						if !ismissing(Categoricalconsts[con, :max])
							indices = string.(ATAmodel.Settings.Bank[!, Symbol.(Categoricalconsts[con, :var])]) .== Categoricalconsts[con, :value]
							indices[ismissing.(indices)] .= false
							indices = findall(indices .== true)
							if Categoricalconsts[con, :max] == 0
								x_forced0[t][indices] .= false
							else
								A[t] = vcat(A[t], zeros(Float64, 1, nItems))
								if ismissing(Categoricalconsts[con, :weight])
									A[t][end, indices] .= 1.0
									push!(b[t], Categoricalconsts[con, :max])
								else
									A[t][end, indices] .= (Categoricalconsts[con, :weight])
									push!(b[t], Categoricalconsts[con, :weight]*Categoricalconsts[con, :max])
								end
							end
						end
					end
				end
			end
			#add quantitative constraints
			for t = 1:ATAmodel.Settings.T
				for var in 1:size(ATAmodel.Constraints[t].sumVars, 1)
					vals = copy(ATAmodel.Settings.Bank[!, ATAmodel.Constraints[t].sumVars[var]])
					vals[findall(ismissing.(vals))].= 0.0
					if !ismissing(ATAmodel.Constraints[t].sumVars_min[var])
						A[t] = vcat(A[t], .-vals')
						push!(b[t], -ATAmodel.Constraints[t].sumVars_min[var])
					end
					if !ismissing(ATAmodel.Constraints[t].sumVars_max[var])
						A[t] = vcat(A[t], vals')
						push!(b[t], ATAmodel.Constraints[t].sumVars_max[var])
					end
				end
			end
			message[2] = message[2] * "- Sum variables constraints added. \n"
			for t = 1:ATAmodel.Settings.T
				writedlm("OPT/A_$t.csv", A[t])
				writedlm("OPT/b_$t.csv", b[t])
			end
			writedlm("OPT/x_forced0.jl", x_forced0)
		end
		#update model
		ATAmodel.Settings.forced0 = x_forced0
		for t = 1:ATAmodel.Settings.T
			ATAmodel.Constraints[t].catConstrb = b[t]
			ATAmodel.Constraints[t].catConstrA = A[t]
		end
		JLD2.@save "OPT/ATAmodel.jld2" model
		message[1] = "success"
		message[2] = message[2] * "- Constraints added. \n"
	end
	return message
end

function AddOverlaps!(ATAmodel::model; overlap_file = "OverlapMatrix.csv", overlap_delim=";")
	T = ATAmodel.Settings.T
	nItems = ATAmodel.Settings.nItems
	if !isfile(overlap_file)
		return ["danger", string(overlap_file, " not found.")]
	else
		opMatrix = Matrix{Int64}(CSV.read(overlap_file, delim = overlap_delim, header = false))
		if size(opMatrix, 1)>0
			#ol_max = Vector{Vector{Int64}}(undef, T)
			# for t = 1:T
			# 	ol_max[t] = opMatrix[t, setdiff(collect(1:T), t)]
			# end
			opMatrix = opMatrix[1:T, 1:T]
			writedlm("OPT/overlap.jl", opMatrix)
			ATAmodel.Settings.olMax = opMatrix
		else
			ATAmodel.Settings.olMax = ones(ATAmodel.Settings.nItems,ATAmodel.Settings.nItems).*ATAmodel.Settings.nItems
			return ["success", string("- No lines in ", overlap_file, ", Maximum overlap set at ", ATAmodel.Settings.nItems,".\n")]
		end
		JLD2.@save "OPT/ATAmodel.jld2" ATAmodel
		return ["success", "- Maximum overlap constrained.\n"]
	end
end

function AddExpScore!(ATAmodel::model)
	T = ATAmodel.Settings.T
	nItems = ATAmodel.Settings.nItems
	nGroups = ATAmodel.Settings.nGroups
	ExSVar = [ATAmodel.Constraints[t].ExS.Var for t = 1:ATAmodel.Settings.T]
	ICF = Vector{Matrix{Float64}}(undef, T)
	for t = 1:T
		if ExSVar[t] == Symbol("")
			ICF[t] = ItemCharFun(ATAmodel.Settings.IRT.parameters, ATAmodel.Constraints[t].ExS.Pts, model=ATAmodel.Settings.IRT.model, parametrization=ATAmodel.Settings.IRT.parametrization, D=ATAmodel.Settings.IRT.D)[1][:,:,1]# K[t] x I
		else
			ICF[t] = zeros(Float64,1,nItems)
			ICF[t][1,:] = ATAmodel.Settings.Bank[!, ExSVar[t]]
		end
	end
	for t = 1:T
		ATAmodel.Constraints[t].ExS.Val = ICF[t]
	end
	# for t = 1:ATAmodel.Settings.T
	# 	writedlm("OPT/A_$t.csv", ATAmodel.Constraints[t].catConstrA)
	# 	writedlm("OPT/b_$t.csv", ATAmodel.Constraints[t].catConstrb)
	# end
	JLD2.@save "OPT/ICF.jld2" ICF
	return ["success", "- Expected Score constrained."]
end

function GroupByFriendSet!(ATAmodel::model) #last
	nItems = ATAmodel.Settings.nItems
	#only works for categorical variables and item use, all the other contraitns need expansion by FSitems
	if ATAmodel.Settings.nFS == 0
		return ["danger", "No friend sets found, run AddFriendSets!(model) before."]
	end
	if size(ATAmodel.Constraints[1].catConstrb, 1) == 0
		return ["danger", "No constraints to group, run AddConstr!(model) before."]
	else
		A = Vector{Matrix{Float64}}(undef, ATAmodel.Settings.T)
		b = Vector{Vector{Float64}}(undef, ATAmodel.Settings.T)
		A_new = Vector{Matrix{Float64}}(undef, ATAmodel.Settings.T)
		for t = 1:ATAmodel.Settings.T
			for t = 1:ATAmodel.Settings.T
				A[t] = ATAmodel.Constraints[t].catConstrA
				b[t] = ATAmodel.Constraints[t].catConstrb
			end
		end
		nFS = ATAmodel.Settings.nFS
		for t = 1:ATAmodel.Settings.T
			A_new[t] = Matrix{Float64}(undef, size(A[t], 1), nFS)
		end
		for fs = 1 : nFS
			for t = 1:ATAmodel.Settings.T
				A_new[t][:, fs] = sum(A[t][:, ATAmodel.Settings.FS.Items[fs]], dims = 2)
			end
		end
		for t = 1:ATAmodel.Settings.T
			writedlm("OPT/A_$t.csv", A_new[t])
			writedlm("OPT/b_$t.csv", b[t])
			ATAmodel.Constraints[t].catConstrA = A_new[t]
		end
		open("OPT/FSItems.jl", "w") do f
			write(f, "FSItems = Vector{Vector{Int64}}(undef, $nFS)\n")
			for fs = 1:nFS
				write(f, string("FSItems[$fs] =", ATAmodel.Settings.FS.Items[fs], "\n"))
			end
		end
		open("OPT/Settings.jl", "a") do f
			write(f, "group_by_fs = true\n\n")
		end
		#transform forced0
		x_forced0 = ATAmodel.Settings.forced0
		x_forced0_new = Vector{Vector{Bool}}(undef, ATAmodel.Settings.T)
		for v = 1 : ATAmodel.Settings.T
			x_forced0_new[v] = fill(true, nFS)
			for i = 1:ATAmodel.Settings.nItems
				if x_forced0[v][i] == false
					for fs = 1:nFS
						if any(i .== ATAmodel.Settings.FS.Items[fs])
							x_forced0_new[v][fs] = false
						end
					end
				end
			end
		end
		#update ATAmodel.forced0
		ATAmodel.Settings.forced0 = x_forced0_new
		writedlm("OPT/x_forced0.jl", x_forced0_new)
		#item use
		ItemUse_min = ATAmodel.IU.Min
		ItemUse_max = ATAmodel.IU.Max
		ItemUse_min_new = zeros(Int, nFS)
		ItemUse_max_new = zeros(Int, nFS)
		for fs = 1:nFS
			ItemUse_min_new[fs] = Int(maximum(ATAmodel.IU.Min[ATAmodel.Settings.FS.Items[fs]]))
			ItemUse_max_new[fs] = Int(minimum(ATAmodel.IU.Max[ATAmodel.Settings.FS.Items[fs]]))
		end
		#enemy sets

		#update model
		ATAmodel.IU.Max = ItemUse_max_new
		ATAmodel.IU.Min = ItemUse_min_new
		open("OPT/Settings.jl", "a") do f
			write(f, "ItemUse_min = $ItemUse_min_new\n\n")
			write(f, "ItemUse_max = $ItemUse_max_new\n\n")
		end
		return ["success", string("- Grouped in ", nFS, " Friend sets")]
	end
end

function AddObjFun!(ATAmodel::model)
	message = ["",""]
	T = ATAmodel.Settings.T
	nItems = ATAmodel.Settings.nItems
	ThetasOpt = ATAmodel.Obj.OptPts
	if ATAmodel.Settings.OptType == "MAXIMIN" || ATAmodel.Settings.OptType == "CC"
		IIF = Vector{Array{Float64, 2}}(undef, T)
		ICF = Vector{Array{Float64, 2}}(undef, T)
		K = zeros(Int, T)
		for t = 1:T
			K[t] = size(ThetasOpt[t], 1)
			IIF[t] = zeros(K[t], nItems)
			ICF[t] = zeros(K[t], nItems)
			for k = 1:K[t]
				IIF[t][k, :] = ItemInfoFun(ATAmodel.Settings.IRT.parameters, ThetasOpt[t][k], model = ATAmodel.Settings.IRT.model, parametrization = ATAmodel.Settings.IRT.parametrization, D = ATAmodel.Settings.IRT.D)# K[t] x I
				ICF[t][k, :] = ItemCharFun(ATAmodel.Settings.IRT.parameters, ThetasOpt[t][k], model = ATAmodel.Settings.IRT.model, parametrization = ATAmodel.Settings.IRT.parametrization, D = ATAmodel.Settings.IRT.D)[1][:,:,1] # K[t] x I
			end
		end
		JLD2.@save "OPT/IIF.jld2" IIF
		if !isfile("OPT/ICF.jld2")
			JLD2.@save "OPT/ICF.jld2" ICF
		end
		if ATAmodel.Settings.OptType == "CC"
			R = ATAmodel.Obj.AuxInt
			K = zeros(Int, T)
			IIF = Vector{Array{Float64, 3}}(undef, T)
			ICF = Vector{Array{Float64, 3}}(undef, T)
			JLD2.@load "BSPar.jld2" BSPar
			BSa = Matrix(BSPar[2])[:, 2:end]
			BSb = Matrix(BSPar[1])[:, 2:end]
			for t = 1:T
				K[t] = size(ThetasOpt[t], 1)
				IIF[t] = zeros(K[t], nItems, R)
				ICF[t] = zeros(K[t], nItems, R)
				for r = 1:R
					if ATAmodel.Settings.IRT.model == "1PL"
						df = DataFrame(b = BSb[:, r]) #nqp values in interval\r\n",
					elseif ATAmodel.Settings.IRT.model == "2PL"
						df = DataFrame(a = BSa[:, r], b = BSb[:, r]) #nqp values in interval\r\n",
					elseif ATAmodel.Settings.IRT.model == "3PL"
						df = DataFrame(a = BSa[:, r], b = BSb[:, r], c = BSc[:, r])
					end
					for k = 1:K[t]
						IIF[t][k, :, r] = ItemInfoFun(df, ThetasOpt[t][k]; model = (ATAmodel.Settings.IRT.model)) # K[t] x I x R
						ICF[t][k, :, r] = ItemCharFun(df, ThetasOpt[t][k]; model = (ATAmodel.Settings.IRT.model)) # K[t] x I x R
					end
				end
			end
			JLD2.@save "OPT/IIF_CC.jld2" IIF
			JLD2.@save "OPT/ICF_CC.jld2" ICF
			message = ["success","- MAXIMIN CC objective function applied.\n"]
		else
			message = ["success","- MAXIMIN objective function applied.\n"]
		end
	else
		return 	["danger","- Only MAXIMIN and MAXIMIN CC are supported. \n"]
	end
	open("OPT/Settings.jl", "a") do f
		write(f, "K = $K\n\n")
	end
	return message
end

function LoadDesign!(design::Matrix{Any}, ATAmodel::model)
	ATAmodel.Output.Design = design
end
