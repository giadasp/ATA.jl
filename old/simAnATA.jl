include("evalATA.jl")

#optimize with simulated annealing
function siman!(ATAmodel::model; starting_design = Matrix{Float64}(undef,0,0), start_temp= 0.1, geom_temp=0.1, max_time = 1000.00, results_folder= "RESULTS", n_item_sample = 1000, n_test_sample = 2, conv_max = 2, verbosity = 2, opt_feas = 0.0, n_fill = 1 , feas_nh = 0, opt_nh = 5)
	message=""
	if !(results_folder in readdir())
		mkdir(results_folder)
	else
		message *= string("You have already a folder with this name, files in ", results_folder," will be overwritten.\n")
	end

	if ATAmodel.Settings.OptType == "MAXIMIN"
		JLD2.@load "OPT/IIF.jld2" IIF
	elseif ATAmodel.Settings.OptType == "CC"
		JLD2.@load "OPT/IIF_CC.jld2" IIF
	elseif ATAmodel.Settings.OptType == ""
		opt_nh = 0
	end
	##
	nNH = feas_nh + opt_nh
	opt_nh = feas_nh + 1
	if feas_nh == 0
		fF = false
	else
		fF = true
	end
	ATAmodel.Settings.nFS = size(ATAmodel.Settings.FS.Items, 1)
	T = ATAmodel.Settings.T
	nFS = ATAmodel.Settings.nFS
	if nFS == 0
		nFS = ATAmodel.Settings.nItems
	end
	nItems = ATAmodel.Settings.nItems
	FScounts = ATAmodel.Settings.FS.Counts * ones(Float64, T)'
	iu⁺ = 0
	start_time = copy(time())
	nacc = 0 # total accepted trials
	t = copy(start_temp) # temperature - will initially rise or fall to cover parameter space. Then it will fall
	finish = 0 # convergence indicator 0 (failure), 1 (normal success), or 2 (convergence but near bounds)
	f_evals = 0
	hline = " = "^80
	f_star = typemax(Float64) * ones(2)
	NH = 1

	#starting design check
	if size(starting_design,1)>0
		if (size(starting_design,1) != nFS || size(starting_design,2)!=ATAmodel.Settings.T)
			message *= "- Starting design must be of size: (nItems x T).\n"
			return message
		end
		if (any(starting_design!=0 && starting_design!=1))
			message *= "- Starting design must contain only 1 or 0.\n"
			return message
		end
		x₀ = Float64.(starting_design)
	else
		x₀ = zeros(Float64, nFS, ATAmodel.Settings.T)
	end

	NH⁺ = Neighbourhood(x₀, Inf, zeros(T), 1e6 * ones(T), 1e6 * ones(T), zero(Float64))
	bestNH = 1
	if size(x₀, 1) == 0
		NH⁺.x = zeros(Float64, nFS, T)
	end
	#compute f
	iu = sum(NH⁺.x, dims = 2) - ATAmodel.IU.Max
	iu = iu[iu .> 0]
	if size(iu, 1) == 0
		NH⁺.iu = 0
	else
		NH⁺.iu = sum(iu)
	end
	t = copy(start_temp)
	for v2 = 1:T
		NH⁺.infeas[v2], x_Iᵥ = checkFeas(ATAmodel.Settings.FS, ATAmodel.Constraints[v2], NH⁺.x[:, v2], nFS, nItems, v2)
		NH⁺.ol = evalOverlap(NH⁺.x, FScounts, ATAmodel.Settings.olMax, T, NH⁺.ol)
		if fF == false
			if ATAmodel.Settings.OptType == "MAXIMIN"
				NH⁺.obj[v2] = evalTIFMMv(x_Iᵥ, IIF[v2])
			elseif ATAmodel.Settings.OptType == "CC"
				NH⁺.obj[v2] = evalTIFCCv(x_Iᵥ, IIF[v2];α = ATAmodel.Obj.AuxFloat)
			end
		end
	end
	NH⁺.f = comp_f(NH⁺, opt_feas)
	f⁺ = copy(NH⁺.f)
	if (Distributed.nprocs() > 1)
		processors = collect(1:(Distributed.nprocs() - 1))
	else
		processors = [one(Int64)]
	end
	NHs = copy(processors)
	ATAmodel.Output.Neighbourhoods = [Neighbourhood() for n = 1:NHs[end]]
	for nh = 1:NHs[end]
		ATAmodel.Output.Neighbourhoods[nh] = mycopy(NH⁺, ATAmodel.Output.Neighbourhoods[nh])
	end
	round = 1
	while finish == 0
		Distributed.@sync Distributed.@distributed for p in processors
			nh_tot = (NHs[end] * (round - 1)) + NHs[p]
			#println("f was ", ATAmodel.Output.Neighbourhoods[NHs[p]].f)
			NH_proc = Neighbourhood()
			NH_proc = mycopy(ATAmodel.Output.Neighbourhoods[NHs[p]], NH_proc)
			NH_proc = analyseNH(NH_proc, ATAmodel, IIF;fF = fF, n_fill = n_fill, opt_feas = opt_feas, conv_max = conv_max, start_temp = t, geom_temp = geom_temp, n_item_sample = n_item_sample, n_test_sample = n_test_sample, verbosity = verbosity, start_time = start_time, max_time = max_time)
			open(string(results_folder,"/neigh_", nh_tot, "_x.csv"), "w") do io
				DelimitedFiles.writedlm(io, NH_proc.x)
			end
			open(string(results_folder,"/neigh_", nh_tot, "_infeas.csv"), "w") do io
				DelimitedFiles.writedlm(io, NH_proc.infeas)
			end
			open(string(results_folder,"/neigh_", nh_tot, "_obj.csv"), "w") do io
				DelimitedFiles.writedlm(io, NH_proc.obj)
			end
			open(string(results_folder,"/neigh_", nh_tot, "_ol.csv"), "w") do io
				DelimitedFiles.writedlm(io, NH_proc.ol)
			end
			open(string(results_folder,"/neigh_", nh_tot, "_iu.csv"), "w") do io
				DelimitedFiles.writedlm(io, NH_proc.iu)
			end
			open(string(results_folder,"/neigh_", nh_tot, "_f.csv"), "w") do io
				DelimitedFiles.writedlm(io, NH_proc.f)
			end
			println("neighbourhood ", nh_tot, " fully explored, increase temperature")
		end
		for p in processors
			nh_tot = (NHs[end] * (round - 1)) + NHs[p]
			ATAmodel.Output.Neighbourhoods[NHs[p]].x = DelimitedFiles.readdlm(string(results_folder,"/neigh_", nh_tot, "_x.csv"), '\t')
			ATAmodel.Output.Neighbourhoods[NHs[p]].obj = DelimitedFiles.readdlm(string(results_folder,"/neigh_", nh_tot, "_obj.csv"), '\t')[:, 1]
			ATAmodel.Output.Neighbourhoods[NHs[p]].infeas = DelimitedFiles.readdlm(string(results_folder,"/neigh_", nh_tot, "_infeas.csv"), '\t')[:, 1]
			ATAmodel.Output.Neighbourhoods[NHs[p]].ol = DelimitedFiles.readdlm(string(results_folder,"/neigh_", nh_tot, "_ol.csv"), '\t')[:, 1]
			ATAmodel.Output.Neighbourhoods[NHs[p]].iu = DelimitedFiles.readdlm(string(results_folder,"/neigh_", nh_tot, "_iu.csv"), '\t')[1, 1]
			ATAmodel.Output.Neighbourhoods[NHs[p]].f = DelimitedFiles.readdlm(string(results_folder,"/neigh_", nh_tot, "_f.csv"), '\t')[1, 1]
		end
		#NH+= 1
		#find best nh
		f_loc, nh_loc = findmin([ATAmodel.Output.Neighbourhoods[n].f for n = 1:NHs[end]])
		if f_loc<= f⁺
			bestNH = (NHs[end] * (round - 1)) + nh_loc
			f⁺ = copy(f_loc)
			NH⁺ = mycopy(ATAmodel.Output.Neighbourhoods[nh_loc], NH⁺)
		end

		#NHs = NHs.+Distributed.nprocs()
		if NHs[end] * round + 1 > feas_nh
			fF = false
		end
		if NHs[end] * round + 1 > nNH
			finish = 1
		else
			#perturbate the solution and go to warmup
		end
		round+= 1
		# f_evals+= f_evals_NH
		# if f_evals >= max_evals
		# 	finish = 2
		# end
		if verbosity >= 1
			println(hline)
			println("Results")
			if (finish == 1)
				println(" == > Maximum number of Neighbourhoods explored <== ")
				# elseif (finish == 2) #time max reached
				# 	println(" == > Maximum number of evaluations reached <== ")
			elseif (finish == 3) #evals max reached
				println(" == > Maximum time reached <== ")
			end
			Printf.@printf("\n Obj. value:	%16.10f", NH⁺.f)
			Printf.@printf("\n Optimality:	%16.10f", minimum(NH⁺.obj))
			Printf.@printf("\n Infeasibility:	%16.10f", sum(NH⁺.infeas + NH⁺.ol) + NH⁺.iu)
			Printf.@printf("\n Elapsed Time:	%16.10f", time() - start_time)
			Printf.@printf("\n Best Neighbourhood:	%16.10f", bestNH)
			println("\n")
			println(hline)
		end
		if time() - start_time >= max_time
			ATAmodel.Output.ElapsedTime = time() - start_time
			finish = 3
		end
		if finish > 0
			#f⁺, NH = findmin([ATAmodel.Output.Neighbourhoods[n].f for n = 1:nNH])
			# x⁺ = ATAmodel.Output.Neighbourhoods[NH].x
			# infeas⁺ = ATAmodel.Output.Neighbourhoods[NH].infeas+ATAmodel.Output.Neighbourhoods[NH].ol
			# iu⁺ = ATAmodel.Output.Neighbourhoods[NH].iu
			# TIF⁺ = ATAmodel.Output.Neighbourhoods[NH].obj
			ATAmodel.Output.Design = NH⁺.x
			ATAmodel.Output.Feas = NH⁺.infeas + NH⁺.ol
			ATAmodel.Output.f = NH⁺.f
			ATAmodel.Output.TIF = NH⁺.obj
		end
		# j = 1
		# for v = 1:ATAmodel.Settings.T
		# 	for i = 1:ATAmodel.Settings.nFS
		# 		if NH⁺.x[i, v] == one(Float64)
		# 			if j == 1
		# 				NH⁺.x[i, v] = zero(Float64)
		# 				j+= 1
		# 			else
		# 				j - = 1
		# 			end
		# 		end
		# 	end
		# end
	end #end of finish
	JLD2.@save string(results_folder,"/ATAmodel.jld2") ATAmodel
	writedlm(string(results_folder,"/design.csv"), reshape(ATAmodel.Output.Design, nFS, T))
	open(string(results_folder,"/ResultsATA.txt"), "w") do io
		write(io, "tests")
		write(io, "\r\n")
		writedlm(io, collect(1:T)', ", ")
		write(io, "\r\n")
		write(io, "f⁺")
		write(io, "\r\n")
		writedlm(io, ATAmodel.Output.f)
		write(io, "\r\n")
		write(io, "infeasibility")
		write(io, "\r\n")
		writedlm(io, ATAmodel.Output.Feas', ", ")
		write(io, "\r\n")
		write(io, "Item use")
		write(io, "\r\n")
		write(io, string(NH⁺.iu))
		write(io, "\r\n")
		write(io, "TIF")
		write(io, "\r\n")
		writedlm(io, ATAmodel.Output.TIF', ", ")
		write(io, "\r\n")
		write(io, "elapsed Time for optimization")
		write(io, "\r\n")
		write(io, string(ATAmodel.Output.ElapsedTime))
		write(io, "\r\n")
	end
	return message
end
