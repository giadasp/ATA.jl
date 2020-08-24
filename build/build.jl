#cc metho

function start_ATA()
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
	ATAmodel = Model()
	JLD2.@save "OPT/ATAmodel.jld2" ATAmodel
	return ATAmodel
end

function load_settings!(ATAmodel::Model; settings_file = "settingsATA.jl", bank_file = "bank.csv", bank_delim = ";")
	message = ["", ""]
	if isfile(bank_file)
		ATAmodel.settings.bank = CSV.read(bank_file, delim = bank_delim)
		message[2] = message[2] * "- Item bank file read.\n"
	else
		return ["danger", "Not a valid file for bank"]
	end
	if !isfile(settings_file)
		return ["danger", "Settings file not valid"]
	else
		include(settings_file)
		ATAmodel.settings.n_items = Inputs.n_items
		ATAmodel.settings.T = Int(sum(Inputs.T))
		ATAmodel.settings.Tg = Inputs.T
		#initialize cosntraints
		ATAmodel.constraints = [Constraint() for t = 1:ATAmodel.settings.T]
		ATAmodel.settings.n_groups = size(Inputs.groups, 1)
		x_forced0 = Vector{Vector{Bool}}(undef, ATAmodel.settings.T)
		for t = 1:ATAmodel.settings.T
			x_forced0[t] = fill(true, ATAmodel.settings.n_items)
			ATAmodel.constraints[t].constr_A = zeros(Float64, 0, ATAmodel.settings.n_items)
		end

		open("OPT/Settings.jl", "w") do f
			write(f, "#Settings \n\n")

				ATAmodel.settings.IRT.model = Inputs.IRT_model
				ATAmodel.settings.IRT.parameters = DataFrames.DataFrame(ATAmodel.settings.bank[!, Symbol.(Inputs.IRT_parameters)])
				ATAmodel.settings.IRT.parametrization = Inputs.IRT_parametrization
				ATAmodel.settings.IRT.D = Inputs.IRT_D

				if ATAmodel.settings.IRT.model == "1PL"
					DataFrames.DataFrames.rename!(ATAmodel.settings.IRT.parameters, [:b])#nqp values in interval\r\n",
				elseif ATAmodel.settings.IRT.model == "2PL"
					DataFrames.DataFrames.rename!(ATAmodel.settings.IRT.parameters, [:a, :b]) #nqp values in interval\r\n",
				elseif ATAmodel.settings.IRT.model == "3PL"
					DataFrames.DataFrames.rename!(ATAmodel.settings.IRT.parameters, [:a, :b, :c]) #nqp values in interval\r\n",
				else
					return ["danger","Only 1PL, 2PL and 3PL IRT models are allowed."]
				end
				CSV.write("OPT/IRT_parameters.csv", ATAmodel.settings.IRT.parameters)
				message[2] = message[2] * "- IRT item parameters loaded.\n"
			if Inputs.enemy_sets_var != String[]
				ATAmodel.settings.ES.var = Symbol.(Inputs.enemy_sets_var)
				val = Symbol.(Inputs.enemy_sets_var)
				write(f, "enemy_sets_var = $val\n\n")
				message[2] = message[2] * "- Variable for Enemy sets loaded.\n"
			end
			if Inputs.friend_sets_var != String[]
				ATAmodel.settings.FS.var = Symbol.(Inputs.friend_sets_var)
				val = Symbol.(Inputs.friend_sets_var)
				write(f, "ATAmodel.settings.FS.var = $val\n\n")
				message[2] = message[2] * "- Variable for Friend sets loaded.\n"
			end

			if Inputs.length_min != Int64[]
				lengthmin = zeros(Int64, ATAmodel.settings.T)
				lengthweight = ones(Int64, ATAmodel.settings.T)
				t1 = 1
				for g = 1:ATAmodel.settings.n_groups
					for t = 1:ATAmodel.settings.Tg[g]
						lengthmin[t1] = Int(Inputs.length_min[g])
						lengthweight[t1] = Inputs.length_weight[g]
						t1+= 1
					end
				end
				for t = 1:ATAmodel.settings.T
					ATAmodel.constraints[t].constr_A = vcat(ATAmodel.constraints[t].constr_A, (-lengthweight[t]).*ones(Float64, ATAmodel.settings.n_items)')
					ATAmodel.constraints[t].constr_b = vcat(ATAmodel.constraints[t].constr_b, -lengthmin[t]*lengthweight[t])
					ATAmodel.constraints[t].length_min = lengthmin[t]
				end
				write(f, "length_min = $lengthmin\n")
				message[2] = message[2] * "- Minimum length of tests constrained.\n"
			end

			if Inputs.length_max != Int64[]
				lengthmax = zeros(Int64, ATAmodel.settings.T)
				lengthweight = ones(Int64, ATAmodel.settings.T)
				t1 = 1
				for g = 1:ATAmodel.settings.n_groups
					for t = 1:ATAmodel.settings.Tg[g]
						lengthmax[t1] = Int(Inputs.length_max[g])
						lengthweight[t1] = Inputs.length_weight[g]
						t1+= 1
					end
				end
				for t = 1:ATAmodel.settings.T
					ATAmodel.constraints[t].constr_A = vcat(ATAmodel.constraints[t].constr_A, (lengthweight[t]).*ones(ATAmodel.settings.n_items)')
					ATAmodel.constraints[t].constr_b = vcat(ATAmodel.constraints[t].constr_b, lengthmax[t]*lengthweight[t])
					ATAmodel.constraints[t].length_max = lengthmax[t]
				end
				write(f, "length_max = $lengthmax\n")
				message[2] = message[2] * "- Maximum length of tests constrained.\n"
			end

			if Inputs.expected_score_var != String[]
				t1 = 1
				for g = 1:ATAmodel.settings.n_groups
					for t = 1:ATAmodel.settings.Tg[g]
						ATAmodel.constraints[t1].expected_score.var = Symbol.(Inputs.expected_score_var[g])
						t1+= 1
					end
				end
			end
			t1 = 1
			if Inputs.expected_score_pts != Float64[]
				for g = 1:ATAmodel.settings.n_groups
					for t = 1:ATAmodel.settings.Tg[g]
						ATAmodel.constraints[t1].expected_score.pts = Inputs.expected_score_pts[g]
						t1+= 1
					end
				end
			end
			message[2] = message[2] * "- Expected score variable and points loaded.\n"

			if Inputs.expected_score_min != Float64[]
				t1 = 1
				for g = 1:ATAmodel.settings.n_groups
					for t = 1:ATAmodel.settings.Tg[g]
						ATAmodel.constraints[t1].expected_score.min = Inputs.expected_score_min[g]
						t1+= 1
					end
				end
				message[2] = message[2] * "- Minimum expected score constrained.\n"
			end
			if Inputs.expected_score_max != Float64[]
				t1 = 1
				for g = 1:ATAmodel.settings.n_groups
					for t = 1:ATAmodel.settings.Tg[g]
						ATAmodel.constraints[t1].expected_score.max = Inputs.expected_score_max[g]
						t1+= 1
					end
				end
				message[2] = message[2] * "- Maximum expected score constrained.\n"
			end
			if Inputs.sum_vars != Vector{Vector{String}}(undef,0)
				t1 = 1
				for g = 1:ATAmodel.settings.n_groups
					if !ismissing(sum_vars_min[g])
						for v = 1:length(Inputs.sum_vars[g])
							var = ATAmodel.settings.bank[Symbol(Inputs.sum_vars[g][v])]
							var[ismissing.(var)].= zero(Float64)
							#min
							if !ismissing(sum_vars_min[g][v])
								for t = 1:ATAmodel.settings.Tg[g]
									if g>1
										t1 = (g-1)*Tg[g-1]+t
									else
										t1 = copy(t)
									end
									ATAmodel.constraints[t1].constr_A = vcat(ATAmodel.constraints[t1].constr_A, .-var')
									ATAmodel.constraints[t1].constr_b = vcat(ATAmodel.constraints[t1].constr_b, -sum_vars_min[g])
								end
							end
							#min
							if !ismissing(sum_vars_max[g][v])
								for t = 1:ATAmodel.settings.Tg[g]
									if g>1
										t1 = (g-1)*Tg[g-1]+t
									else
										t1 = copy(t)
									end
									ATAmodel.constraints[t1].constr_A = vcat(ATAmodel.constraints[t1].constr_A, var')
									ATAmodel.constraints[t1].constr_b = vcat(ATAmodel.constraints[t1].constr_b, sum_vars_max[g])
								end
							end
						end
					end
				end
				message[2] = message[2] * string("- Sum of variables", Inputs.sum_vars , " will be constrained.\n")
			end
			if size(Inputs.item_use_min, 1) > 0
				ATAmodel.IU.min = Inputs.item_use_min
				message[2] = message[2] * "- Minimum item use constrained.\n"
			end
			if size(Inputs.item_use_max, 1) > 0
				ATAmodel.IU.max = Inputs.item_use_max
				for v = 1:ATAmodel.settings.T
					x_forced0[v][findall(ATAmodel.IU.max.<1)] .= false
				end
				message[2] = message[2] * "- Maximum item use constrained.\n"
			end
			if Inputs.obj_type != ""
				ATAmodel.obj.type = Inputs.obj_type
				ATAmodel.obj.fun = Inputs.obj_fun
				ATAmodel.obj.args = Inputs.obj_args
				ATAmodel.obj.aux_int = Inputs.aux_int
				ATAmodel.obj.aux_float = Inputs.aux_float
				message[2] = message[2] * "- Optimization type and objective function loaded.\n"
			end
			if size(Inputs.points, 1) > 0
				ATAmodel.obj.points = Vector{Vector{Float64}}(undef, ATAmodel.settings.T)
				t1 = 1
				for g = 1:ATAmodel.settings.n_groups
					for t = 1:ATAmodel.settings.Tg[g]
						ATAmodel.obj.points[t1] = Inputs.points[g]
						t1+= 1
					end
				end
				message[2] = message[2] *  "- Points to optimize loaded.\n"
			end
			
			message[2] = message[2] * "- Auxiliars vars for optimization loaded.\n"
			#fictiuos friendSets
			ATAmodel.settings.FS.counts = ones(ATAmodel.settings.n_items)
			ATAmodel.settings.FS.sets = string.(collect(1:ATAmodel.settings.n_items))
			ATAmodel.settings.FS.items = [[i] for i = 1:ATAmodel.settings.n_items]
			ATAmodel.settings.forced0 = x_forced0
			if Inputs.categories != String[]
				val = Symbol.(Inputs.categories)
				ATAmodel.output.categories = copy(val)
				write(f, "categories = $val\n\n")
				message[2] = message[2] * "- categories for output loaded.\n"
			end
		end
		message[2] = message[2] * string("assemble ", ATAmodel.settings.T, " forms divided in ", ATAmodel.settings.n_groups, " groups.\n")
		#update model
		JLD2.@save "OPT/ATAmodel.jld2" ATAmodel

	end
	message[1] = "success"
	return message
end


function add_friends!(ATAmodel::Model)
	message = ["", ""]
	if !isfile("OPT/Settings.jl")
		return  ["danger","Run load_settings!(model) before!"]
	else
		n_items = ATAmodel.settings.n_items
		FriendSets = string.(unique(vcat([unique(skipmissing(ATAmodel.settings.bank[!, (ATAmodel.settings.FS.var[isv])])) for isv = 1:(size(ATAmodel.settings.FS.var, 1))]...)))
		n_FS = size(FriendSets, 1)
		FS_counts = zeros(Int, n_FS)
		FSItems = [zeros(Int, 0) for i = 1 : n_FS]
		single_items = Vector{Union{Missing,String}}([missing for i = 1 : n_items])
		for i = 1 : n_items
			units = [ATAmodel.settings.bank[i, ATAmodel.settings.FS.var[isv]] for isv = 1:(size(ATAmodel.settings.FS.var, 1))]
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
					FS_counts[f] += 1
				end
			end
		end
		push!(ATAmodel.settings.FS.var, :SINGLE_FS)
		DataFrames.DataFrames.insertcols!(ATAmodel.settings.bank, size(ATAmodel.settings.bank,2), :SINGLE_FS => single_items)
		items_single = findall(.!ismissing.(ATAmodel.settings.bank[!,:SINGLE_FS]))
		FSItems = vcat(FSItems, [[i] for i in items_single])
		FriendSets = vcat(FriendSets, ATAmodel.settings.bank[items_single,:SINGLE_FS])
		FS_counts = vcat(FS_counts, ones(Int64,size(items_single,1)))
		n_FS = size(FS_counts,1)
		DelimitedFiles.writedlm("OPT/FriendSets.csv", FriendSets)

		#update model
		ATAmodel.settings.n_FS = n_FS
		ATAmodel.settings.FS.sets = FriendSets
		ATAmodel.settings.FS.items = FSItems
		ATAmodel.settings.FS.counts = FS_counts
		JLD2.@save "OPT/ATAmodel.jld2" Model
		return ["success", string("- ",n_FS, " friend sets added.")]
	end
end

function add_enemies!(ATAmodel::Model)
	bank = ATAmodel.settings.bank
	n_items = ATAmodel.settings.n_items
	if size(ATAmodel.constraints[1].constr_b, 1)>0
		A = Vector{Matrix{Float64}}(undef, ATAmodel.settings.T)
		b = Vector{Vector{Float64}}(undef, ATAmodel.settings.T)
		for t = 1:ATAmodel.settings.T
			A[t] = ATAmodel.constraints[t].constr_A
			b[t] = ATAmodel.constraints[t].constr_b
		end

	else
		A = Vector{Matrix{Float64}}(undef, ATAmodel.settings.T)
		b = Vector{Vector{Float64}}(undef, ATAmodel.settings.T)
		for t = 1:ATAmodel.settings.T
			A[t] = Matrix{Float64}(undef, 0, ATAmodel.settings.n_items)
			b[t] = Vector{Float64}(undef, 0)
		end
	end

	if size(ATAmodel.settings.ES.var, 1)>0
		enemy_sets_var = ATAmodel.settings.ES.var
		EnemySets = unique(vcat([unique(skipmissing(bank[!, (enemy_sets_var[isv])])) for isv = 1:(size(enemy_sets_var, 1))]...))
		DelimitedFiles.writedlm("EnemySets.csv", EnemySets)
		nES = size(EnemySets, 1)
		sets = Vector{Vector{Int64}}(undef, nES)
		for es = 1:nES
			sets[es] = Int64[]
			for t = 1:ATAmodel.settings.T
				A[t] = vcat(A[t], zeros(ATAmodel.settings.n_items)')
			end
			for i = 1:n_items
				if any(skipmissing(vcat([((bank[!, (enemy_sets_var[isv])][i] == EnemySets[es]) == true) for isv = 1:(size(enemy_sets_var, 1))]...)))
					#ESMaATAmodel.settings.Tg[n_items, es] = 1
					for t = 1:ATAmodel.settings.T
						A[t][end, i] = 1
					end
					sets[es] = vcat(sets[es], i)
				end
			end
			for t = 1:ATAmodel.settings.T
				push!(b[t], 1)
			end
		end
		for t = 1:ATAmodel.settings.T
			DelimitedFiles.writedlm("OPT/A_$t.csv", A[t])
			DelimitedFiles.writedlm("OPT/b_$t.csv", b[t])
		end
		#update model
		ATAmodel.settings.ES.names = EnemySets
		ATAmodel.settings.ES.sets = sets
		for t = 1:ATAmodel.settings.T
			ATAmodel.constraints[t].constr_b = b[t]
			ATAmodel.constraints[t].constr_A = A[t]
		end
		JLD2.@save "OPT/ATAmodel.jld2" Model
		return ["success", string("- ",nES, " enemy sets added. ")]
	end
end

function add_constraints!(ATAmodel::Model; constraints_file = "CategoricalConstraints.csv", constraints_delim = ";")
	message=["", ""]
	n_items = ATAmodel.settings.n_items
	bank = ATAmodel.settings.bank
	A = Vector{Matrix{Float64}}(undef, ATAmodel.settings.T)
	b = Vector{Vector{Float64}}(undef, ATAmodel.settings.T)
	if size(ATAmodel.constraints[1].constr_b, 1)>0
		for t = 1:ATAmodel.settings.T
			A[t] = copy(ATAmodel.constraints[t].constr_A)
			b[t] = copy(ATAmodel.constraints[t].constr_b)
		end
	else
		for t = 1:ATAmodel.settings.T
			A[t] = zeros(Float64, 0, ATAmodel.settings.n_items)
			b[t] = Vector{Float64}(undef, 0)
		end

	end
	if !isfile(constraints_file)
		return ["danger",string(constraints_file," doesn't exist.")]
	else
		x_forced0 = ATAmodel.settings.forced0
		Categoricalconsts = CSV.read(constraints_file, delim = constraints_delim)
		if size(Categoricalconsts, 1)==0
			message[2] = message[2] * string("No lines in file ", constraints_file," nothing added.\n")
		else
			Categoricalconsts.var = Symbol.(Categoricalconsts.var)
			CatCons = Vector{Vector{Float64}}(undef, ATAmodel.settings.T)
			for g = 1:ATAmodel.settings.n_groups
				if g>1
					tests = collect((sum(ATAmodel.settings.Tg[1:(g-1)])+1):(sum(ATAmodel.settings.Tg[1:(g-1)])+ATAmodel.settings.Tg[g]))
				else
					tests = collect(1:ATAmodel.settings.Tg[1])
				end
				GroupCatCons = findall((Categoricalconsts.group.== g))
				for con in GroupCatCons
					for t in tests
						if !ismissing(Categoricalconsts[con, :min])
							indices = string.(bank[!, Symbol.(Categoricalconsts[con, :var])]) .== Categoricalconsts[con, :value]
							indices[ismissing.(indices)] .= false
							indices = findall(indices .== true)
							A[t] = vcat(A[t], zeros(Float64, 1, n_items))
							if ismissing(Categoricalconsts[con, :weight])
								A[t][end, indices] .= -1.0
								push!(b[t], -Categoricalconsts[con, :min])
							else
								A[t][end, indices] .= -(Categoricalconsts[con, :weight])
								push!(b[t], -Categoricalconsts[con, :weight]*Categoricalconsts[con, :min])
							end
						end
						if !ismissing(Categoricalconsts[con, :max])
							indices = string.(ATAmodel.settings.bank[!, Symbol.(Categoricalconsts[con, :var])]) .== Categoricalconsts[con, :value]
							indices[ismissing.(indices)] .= false
							indices = findall(indices .== true)
							if Categoricalconsts[con, :max] == 0
								x_forced0[t][indices] .= false
							else
								A[t] = vcat(A[t], zeros(Float64, 1, n_items))
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
			for t = 1:ATAmodel.settings.T
				for var in 1:size(ATAmodel.constraints[t].sum_vars, 1)
					vals = copy(ATAmodel.settings.bank[!, ATAmodel.constraints[t].sum_vars[var]])
					vals[findall(ismissing.(vals))].= 0.0
					if !ismissing(ATAmodel.constraints[t].sum_vars_min[var])
						A[t] = vcat(A[t], .-vals')
						push!(b[t], -ATAmodel.constraints[t].sum_vars_min[var])
					end
					if !ismissing(ATAmodel.constraints[t].sum_vars_max[var])
						A[t] = vcat(A[t], vals')
						push!(b[t], ATAmodel.constraints[t].sum_vars_max[var])
					end
				end
			end
			message[2] = message[2] * "- Sum variables constraints added. \n"
			for t = 1:ATAmodel.settings.T
				DelimitedFiles.writedlm("OPT/A_$t.csv", A[t])
				DelimitedFiles.writedlm("OPT/b_$t.csv", b[t])
			end
			DelimitedFiles.writedlm("OPT/x_forced0.txt", x_forced0)
		end
		#update model
		ATAmodel.settings.forced0 = x_forced0
		for t = 1:ATAmodel.settings.T
			ATAmodel.constraints[t].constr_b = b[t]
			ATAmodel.constraints[t].constr_A = A[t]
		end
		JLD2.@save "OPT/ATAmodel.jld2" Model
		message[1] = "success"
		message[2] = message[2] * "- constraints added. \n"
	end
	return message
end

function add_overlap!(ATAmodel::Model; overlap_file = "OverlapMatrix.csv", overlap_delim=";")
	T = ATAmodel.settings.T
	n_items = ATAmodel.settings.n_items
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
			DelimitedFiles.writedlm("OPT/overlap.txt", opMatrix)
			ATAmodel.settings.ol_max = opMatrix
		else
			ATAmodel.settings.ol_max = ones(ATAmodel.settings.n_items,ATAmodel.settings.n_items).*ATAmodel.settings.n_items
			return ["success", string("- No lines in ", overlap_file, ", Maximum overlap set at ", ATAmodel.settings.n_items,".\n")]
		end
		JLD2.@save "OPT/ATAmodel.jld2" ATAmodel
		return ["success", "- Maximum overlap constrained.\n"]
	end
end

function add_exp_score!(ATAmodel::Model)
	T = ATAmodel.settings.T
	n_items = ATAmodel.settings.n_items
	n_groups = ATAmodel.settings.n_groups
	expected_score_var = [ATAmodel.constraints[t].expected_score.var for t = 1:ATAmodel.settings.T]
	ICF = Vector{Matrix{Float64}}(undef, T)
	for t = 1:T
		if expected_score_var[t] == Symbol("")
			ICF[t] = item_char(ATAmodel.settings.IRT.parameters, ATAmodel.constraints[t].expected_score.pts, model=ATAmodel.settings.IRT.model, parametrization=ATAmodel.settings.IRT.parametrization, D=ATAmodel.settings.IRT.D)[1][:,:,1]# K[t] x I
		else
			ICF[t] = zeros(Float64,1,n_items)
			ICF[t][1,:] = ATAmodel.settings.bank[!, expected_score_var[t]]
		end
	end
	for t = 1:T
		ATAmodel.constraints[t].expected_score.val = ICF[t]
	end
	# for t = 1:ATAmodel.settings.T
	# 	DelimitedFiles.writedlm("OPT/A_$t.csv", ATAmodel.constraints[t].constr_A)
	# 	DelimitedFiles.writedlm("OPT/b_$t.csv", ATAmodel.constraints[t].constr_b)
	# end
	JLD2.@save "OPT/ICF.jld2" ICF
	return ["success", "- Expected Score constrained."]
end

function group_by_friends!(ATAmodel::Model) #last
	n_items = ATAmodel.settings.n_items
	#only works for categorical variables and item use, all the other contraitns need expansion by FS_items
	if ATAmodel.settings.n_FS == 0
		return ["danger", "No friend sets found, run add_friends!(model) before."]
	end
	if size(ATAmodel.constraints[1].constr_b, 1) == 0
		return ["danger", "No constraints to group, run add_constraints!(model) before."]
	else
		A = Vector{Matrix{Float64}}(undef, ATAmodel.settings.T)
		b = Vector{Vector{Float64}}(undef, ATAmodel.settings.T)
		A_new = Vector{Matrix{Float64}}(undef, ATAmodel.settings.T)
		for t = 1:ATAmodel.settings.T
			for t = 1:ATAmodel.settings.T
				A[t] = ATAmodel.constraints[t].constr_A
				b[t] = ATAmodel.constraints[t].constr_b
			end
		end
		n_FS = ATAmodel.settings.n_FS
		for t = 1:ATAmodel.settings.T
			A_new[t] = Matrix{Float64}(undef, size(A[t], 1), n_FS)
		end
		for fs = 1 : n_FS
			for t = 1:ATAmodel.settings.T
				A_new[t][:, fs] = sum(A[t][:, ATAmodel.settings.FS.items[fs]], dims = 2)
			end
		end
		for t = 1:ATAmodel.settings.T
			DelimitedFiles.writedlm("OPT/A_$t.csv", A_new[t])
			DelimitedFiles.writedlm("OPT/b_$t.csv", b[t])
			ATAmodel.constraints[t].constr_A = A_new[t]
		end
		open("OPT/FSItems.jl", "w") do f
			write(f, "FSItems = Vector{Vector{Int64}}(undef, $n_FS)\n")
			for fs = 1:n_FS
				write(f, string("FSItems[$fs] =", ATAmodel.settings.FS.items[fs], "\n"))
			end
		end
		open("OPT/Settings.jl", "a") do f
			write(f, "group_by_fs = true\n\n")
		end
		#transform forced0
		x_forced0 = ATAmodel.settings.forced0
		x_forced0_new = Vector{Vector{Bool}}(undef, ATAmodel.settings.T)
		for v = 1 : ATAmodel.settings.T
			x_forced0_new[v] = fill(true, n_FS)
			for i = 1:ATAmodel.settings.n_items
				if x_forced0[v][i] == false
					for fs = 1:n_FS
						if any(i .== ATAmodel.settings.FS.items[fs])
							x_forced0_new[v][fs] = false
						end
					end
				end
			end
		end
		#update ATAmodel.forced0
		ATAmodel.settings.forced0 = x_forced0_new
		DelimitedFiles.writedlm("OPT/x_forced0.txt", x_forced0_new)
		#item use
		item_use_min = ATAmodel.IU.min
		item_use_max = ATAmodel.IU.max
		item_use_min_new = zeros(Int, n_FS)
		item_use_max_new = zeros(Int, n_FS)
		for fs = 1:n_FS
			item_use_min_new[fs] = Int(maximum(ATAmodel.IU.min[ATAmodel.settings.FS.items[fs]]))
			item_use_max_new[fs] = Int(minimum(ATAmodel.IU.max[ATAmodel.settings.FS.items[fs]]))
		end
		#enemy sets

		#update model
		ATAmodel.IU.max = item_use_max_new
		ATAmodel.IU.min = item_use_min_new
		open("OPT/Settings.jl", "a") do f
			write(f, "item_use_min = $item_use_min_new\n\n")
			write(f, "item_use_max = $item_use_max_new\n\n")
		end
		return ["success", string("- Grouped in ", n_FS, " Friend sets")]
	end
end

function add_obj_fun!(ATAmodel::Model)
	message = ["",""]
	T = ATAmodel.settings.T
	n_items = ATAmodel.settings.n_items
	obj_points = ATAmodel.obj.points
	if ATAmodel.obj.type == "MAXIMIN" || ATAmodel.obj.type == "CC"
		IRT_parameters = ATAmodel.settings.IRT.parameters
		IRT_model = ATAmodel.settings.IRT.model
		IRT_D = ATAmodel.settings.IRT.D
		IRT_parametrization = ATAmodel.settings.IRT.parametrization
		IIF = Vector{Array{Float64, 2}}(undef, T)
		ICF = Vector{Array{Float64, 2}}(undef, T)
		K = zeros(Int, T)
		for t = 1:T
			K[t] = size(obj_points[t], 1)
			IIF[t] = zeros(K[t], n_items)
			ICF[t] = zeros(K[t], n_items)
			for k = 1:K[t]
				IIF[t][k, :] = item_info(IRT_parameters, obj_points[t][k], model = IRT_model, parametrization = IRT_parametrization, D = IRT_D)# K[t] x I
				ICF[t][k, :] = item_char(IRT_parameters, obj_points[t][k], model = IRT_model, parametrization = IRT_parametrization, D = IRT_D)[1][:,:,1] # K[t] x I
			end
		end
		JLD2.@save "OPT/IIF.jld2" IIF
		if !isfile("OPT/ICF.jld2")
			JLD2.@save "OPT/ICF.jld2" ICF
		end
		if ATAmodel.obj.type == "CC"
			R = ATAmodel.obj.aux_int
			K = zeros(Int, T)
			IIF = Vector{Array{Float64, 3}}(undef, T)
			ICF = Vector{Array{Float64, 3}}(undef, T)
			JLD2.@load "BSPar.jld2" BSPar
			BSa = Matrix(BSPar[2])[:, 2:end]
			BSb = Matrix(BSPar[1])[:, 2:end]
			for t = 1:T
				K[t] = size(obj_points[t], 1)
				IIF[t] = zeros(K[t], n_items, R)
				ICF[t] = zeros(K[t], n_items, R)
				for r = 1:R
					if IRT_model == "1PL"
						df = DataFrames.DataFrame(b = BSb[:, r]) #nqp values in interval\r\n",
					elseif IRT_model == "2PL"
						df = DataFrames.DataFrame(a = BSa[:, r], b = BSb[:, r]) #nqp values in interval\r\n",
					elseif IRT_model == "3PL"
						df = DataFrames.DataFrame(a = BSa[:, r], b = BSb[:, r], c = BSc[:, r])
					end
					for k = 1:K[t]
						IIF[t][k, :, r] = item_info(df, obj_points[t][k]; model = IRT_model, parametrization = IRT_parametrization, D = IRT_D) # K[t] x I x R
						ICF[t][k, :, r] = item_char(df, obj_points[t][k]; model = IRT_model, parametrization = IRT_parametrization, D = IRT_D) # K[t] x I x R
					end
				end
			end
			JLD2.@save "OPT/IIF_CC.jld2" IIF
			JLD2.@save "OPT/ICF_CC.jld2" ICF
			message = ["success","- MAXIMIN CC objective function applied.\n"]
		else
			message = ["success","- MAXIMIN objective function applied.\n"]
		end
		open("OPT/Settings.jl", "a") do f
		write(f, "K = $K\n\n")
	end
	elseif !(ATAmodel.obj.type == "custom" || ATAmodel.obj.type == "")
		return 	["danger","- Only \"MAXIMIN\", \"CC\", \"\" (no objective) and \"custom\" are supported. \n"]
	end
	
	return message
end

function load_design!(design::Matrix{Any}, ATAmodel::Model)
	ATAmodel.output.design = design
end
