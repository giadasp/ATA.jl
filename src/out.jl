
function PrintResults!(ATAmodel; group_by_fs=false, results_folder="RESULTS", plots_out=false)
	if plots_out==true
		pgfplot()
	end
	if size(ATAmodel.Output.Design,1)>0
		IRTpars=ATAmodel.Settings.IRT.parameters
		T=ATAmodel.Settings.T
		nItems=ATAmodel.Settings.nItems
		length_max=[ATAmodel.Constraints[t].length_max  for t=1:T]
		CATEGORIES=ATAmodel.Output.Categories
		#summarize=Main.summarize
		model=ATAmodel.Settings.IRT.model
		nFS=size(ATAmodel.Settings.FS.Items,1)
		if nFS<ATAmodel.Settings.nItems
			n=nFS*T
		else
			n=nItems*T
		end

		if ATAmodel.Settings.OptType=="CC"
			JLD2.@load "OPT/IIF_CC.jld2" IIF
			JLD2.@load "OPT/ICF_CC.jld2" ICF
			IIF_CC=copy(IIF)
			ICF_CC=copy(ICF)
		end
		if ATAmodel.Settings.OptType=="MAXIMIN" || ATAmodel.Settings.OptType=="CC"
			JLD2.@load "OPT/IIF.jld2" IIF
			JLD2.@load "OPT/ICF.jld2" ICF
		end
		if group_by_fs == true
			design = reshape(ATAmodel.Output.Design, nFS, T)
			new_design = zeros(Float64, nItems, T)
			for t=1:T
				new_design[vcat(ATAmodel.Settings.FS.Items[findall(design[:, t] .== one(Float64))]...), t] .= one(Float64)
			end
			design = new_design
			writedlm(string(results_folder, "/designItemLevel.csv"), design)
		else
			design=reshape(ATAmodel.Output.Design,nItems,T)
		end
		#overlap
		olMatrixOUT=zeros(Int64,T,T)
		for t=1:T
			for v=1:T
				olMatrixOUT[t,v]=design[:,t]'*design[:,v]
			end
		end
		#writedlm(string(results_folder,"/olMatrixOUT.csv"),olMatrixOUT)
		#TIF e ICF
		if ATAmodel.Settings.OptType == "MAXIMIN" ||  ATAmodel.Settings.OptType == "CC"
			if isfile("simPool.csv")
				simPool = CSV.read("simPool.csv")
			else
				simPool = Float64[]
			end
			if plots_out == true
                ThetasPlot=collect(range(-4, stop = 4, length=101)) #nqp values in interval/r/n",
                IIFplot=ItemInfoFun(ATAmodel.Settings.IRT.parameters, ThetasPlot,
				model = ATAmodel.Settings.IRT.model)
                ICFplot=ItemCharFun(ATAmodel.Settings.IRT.parameters, ThetasPlot,
				model = ATAmodel.Settings.IRT.model)[1][:, :, 1]
                IIFf=Array{Float64,2}(undef,T,101)
                ICFf=Array{Float64,2}(undef,T,101)
                for t in 1:T
                    IIFf[t,:] = [sum(IIFplot[findall(design[:,t] .> 0 ), i])   for i in 1:101]
                    ICFf[t,:] = [sum(ICFplot[findall(design[:,t] .> 0 ), i])   for i in 1:101]
                end
				PlotATA!(ATAmodel,
				IIFf,
				ICFf,
				design,
				simPool = simPool,
				results_folder = results_folder)
				if ATAmodel.Settings.OptType == "CC"
					PlotATA_CC!(ATAmodel, IIFf, ICFf, design,
					simPool = simPool,
					results_folder = results_folder)
				end
			end
		end
		#save values
		if ATAmodel.Settings.OptType=="CC"
			Min=Float64[]
			First=Float64[]
			Median=Float64[]
			Third=Float64[]
			Max=Float64[]
			Alpha=Float64[]
			Estimated=Float64[]
			if isfile("simPool.csv")
				True=Float64[]
			end
			Distributed.@sync for t=1:T
				K=size(ATAmodel.Obj.OptPts[t],1)
				for k=1:K
					IIF_t = (IIF_CC[t][k,:,:]'*design[:, t])
					Min = vcat(Min, minimum(IIF_t))
					First = vcat(First, StatsBase.quantile(IIF_t, [0.25]))
					Median = vcat(Median, StatsBase.quantile(IIF_t, [0.5]))
					Third = vcat(Third, StatsBase.quantile(IIF_t, [0.75]))
					Max = vcat(Max, maximum(IIF_t))
					Alpha = vcat(Alpha, StatsBase.quantile(IIF_t, [ATAmodel.Obj.AuxFloat]))
					Estimated = vcat(Estimated, IIF[t][k,:]'*design[:,t] )
					if isfile("simPool.csv")
						True=vcat(True, ItemInfoFun(simPool,
						ATAmodel.Obj.OptPts[t][k],
						model = ATAmodel.Settings.IRT.model)'*design[:,t])
					end
				end
			end
			TIFatTheta_k = DataFrame(hcat(Min, First, Median, Third, Max, Alpha, Estimated))

			if isfile("simPool.csv")
				TIFatTheta_k[:True]=True
				DataFrames.names!(TIFatTheta_k, Symbol.(["Min", "0.25-Qle", "Median", "0.75-Qle", "Max", L"{\alpha}-Qle", "Estimated", "True"]))
			else
				DataFrames.names!(TIFatTheta_k, Symbol.(["Min", "0.25-Qle", "Median", "0.75-Qle", "Max", L"{\alpha}-Qle", "Estimated"]))
			end
			CSV.write(string(results_folder, "/TIFatTheta_k.csv"), TIFatTheta_k)

		else
			Estimated=Float64[]
			if isfile("simPool.csv")
				True = Float64[]
			end
			for t=1:T
				K=size(ATAmodel.Obj.OptPts[t],1)
				for k=1:K
					Estimated=vcat(Estimated,IIF[t][k,:]'*design[:,t])
					if isfile("simPool.csv")
						True = vcat(True, ItemInfoFun(simPool, ATAmodel.Obj.OptPts[t][k], model = ATAmodel.Settings.IRT.model)'*design[:,t])
					end
				end
			end
			TIFatTheta_k=DataFrame(Estimated=Estimated)
			if isfile("simPool.csv")
				TIFatTheta_k[:True] = True
				DataFrames.names!(TIFatTheta_k, Symbol.(["Estimated", "True"]))
			end
			CSV.write(string(results_folder, "/TIFatTheta_k.csv"), TIFatTheta_k)

		end

		#expected score
		ESprintIRT=Vector{Vector{Float64}}(undef, T)
		for t in 1:T
			ESprintIRT[t]=zeros(size(ICF[t], 1))
			for k = 1 : size(ICF[t], 1)
				ESprintIRT[t][k]=((ICF[t][k, :]'*design[:, t]))[1] / sum(design[:, t])
			end
		end
		#number of items
		n = sum(design, dims = 1)
		designItems=zeros(Int64, Int(maximum(n)), T)
		#designItems=Array{String,2}(nothing,n_max,T)
		for t in 1:T
			designItems[collect(1:Int(n[t])), t] .= sort!(findall(design[:,t] .== 1))#sort!(CODEVar[findall(design[:,t].==1)])
		end
		open(string(results_folder,"/Results.txt"), "w") do io
			write(io, "Length")
			write(io, "\r\n")
			writedlm(io, [1:T, Int.(n)'], "\t")
			write(io, "\r\n")
			write(io, "Design")
			write(io,"\r\n")
			writedlm(io, designItems, "\t")
			write(io, "\r\n")
			write(io, "expected score")
			write(io, "\r\n")
			for t = 1:T
				writedlm(io, vcat(t, round.(ESprintIRT[t], digits = 3))', "\t")
				write(io, "\r\n")
			end
			if size(CATEGORIES,1)>0
				write(io, "Categorical variables")
				write(io, "\r\n")
				for cats in CATEGORIES
					write(io, cats)
					write(io, "\r\n")
					vals1 = setdiff(unique(ATAmodel.Settings.Bank[!, Symbol(cats)]), [missing])
					vals = copy(vals1)
					for t = 1:T
						valscat = Int64[]
						for val in vals1
							valscat_t = ATAmodel.Settings.Bank[!, Symbol(cats)][design[:, t] .== 1] .== val
							valscat_t = valscat_t[.!ismissing.(valscat_t)]
							if size(valscat_t, 1) > 0
								valscat = vcat(valscat, sum(valscat_t))
							else
								valscat = vcat(valscat, 0)
							end
						end
						vals = hcat(vals, valscat)
					end
					writedlm(io, vals, "\t")
					#writedlm(io,hcat(vals,[sum(skipmissing(ATAmodel.Settings.Bank[Symbol(cats)][design[:,t].==1].== i))   for i in vals,t=1:T]),", ")
					write(io, "\r\n")
				end
				write(io, "Sum variables")
				write(io, "\r\n")
				for var in unique(vcat([ATAmodel.Constraints[t].sumVars for t = 1:T]...))
					write(io, var)
					write(io, "\r\n")
					writedlm(io, round.([sum(skipmissing(ATAmodel.Settings.Bank[var] .* design[:,t])) for t = 1:T], digits = 3)', "\t")
					write(io, "\r\n")
				end
				write(io, "\r\n")
				write(io, "Overlap matrix")
				write(io, "\r\n")
				olMatrixOUT = hcat(vcat("", string.(collect(1:T))), vcat(collect(1:T)', olMatrixOUT))
				writedlm(io, olMatrixOUT, "\t")
				write(io, "\r\n")
			end
			# if size(summarize,1)>0
			# write(io,"Quantitative variables")
			# write(io,"\r\n")
			#   for summ in 1:size(summarize[:,1],1)
			# vals=[summarize[summ][2](skipmissing(ATAmodel.Settings.Bank[Symbol(summarize[summ][1])][design[:,t].==1]))   for t=1:T]
			# write(io,string(summarize[summ][1]))
			# write(io,",")
			# write(io,string(summarize[summ][2]))
			# write(io,",")
			# writedlm(io,vals',",")
			# write(io,"\n")
			# end
			# end
			write(io,"Item use")
			write(io,"\r\n")
			iu=[sum(design[i,:])   for i=1:ATAmodel.Settings.nItems]
			ids=sortperm(iu, rev=true)
			sort!(iu, rev=true)
			for i=1:ATAmodel.Settings.nItems
				write(io,string(ids[i]))
				write(io,"\t")
				write(io,string(iu[i]))
				write(io,"\r\n")
			end
		end
	else
		println("No solution found")
	end
end
