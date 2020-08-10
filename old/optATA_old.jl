function tranformFStoItems(nFS::Int64,nItems::Int64,xv::Vector{Float64},xv_taken::Vector{Int64},FSitems::Vector{Vector{Int64}})
	if size(xv_taken,1)==0
		xv_taken=findall(xv.==one(Float64))
	end
	x_I=zeros(Float64,nItems)
	for i in xv_taken
		x_I[FSitems[i]].=one(Float64)
	end
	return x_I
end
#filling OptFeas==0
function WarmUp(x::Vector{Float64},cons::Vector{Float64},f::Float64,v::Int64,idxv::Vector{Int64},idxnotv::Vector{Vector{Int64}},IU::IU,FS::FS,Constraints::Constraint,x_forced0::Vector{Int64})
	best_i=0
	x_old, f_old, infeas_old = copy(x), Inf, copy(cons)
	nFS=size(FS.Counts,1)
	nItems=size(IIFv,2)
	T=size(Constraints.ol_max,1)+1
	xv=x[idxv]
	idxv2=setdiff(setdiff(idxv,x_forced0),idxv[xv.==one(Float64) ])
	worstv=one(Float64)
	idxv2=Random.shuffle!(idxv2)
	rows=size(Constraints.catConstrA,1)
	for i2 in idxv2
		x_new, infeas_new = copy(x), copy(infeas_old)
		x_new[i2]=1.0
		xv=x_new[idxv]
		es=0
		cons=0
		if size(Constraints.catConstrA,1)>0
			cons=LinearAlgebra.BLAS.gemv('N',Constraints.catConstrA,xv)-Constraints.catConstrb
		end
		ItemUse=(FS.Counts .*(evalItemUse(x_new,nFS,T)))-IU.Max
		ol=evalOverlapv(xv,x_new,FS.Counts,Constraints.ol_max,idxv,idxnotv)
		if nItems!=nFS
			x_I=tranformFStoItems(nFS,nItems,xv,Int64[],FS.Items)
		else
			x_I=xv
		end
		if size(Constraints.ExS.Val,1)>0
			es=dot(Constraints.ExS.Val,x_I)/sum(x_I)
			es=max(es-Constraints.ExS.Max,-es+Constraints.ExS.Min)
		end
		cons=vcat(cons,ol,es)
		cons=cons[cons.>0]
		if size(cons,1)>0
			infeas_new[v]=sum(cons)
		else
			infeas_new[v]=0
		end
		f_new=sum(infeas_new)+max(sum(ItemUse[ItemUse.>0]),0)
		if (f_new<f_old)
			x_old, f_old, infeas_old = copy(x_new), copy(f_new), copy(infeas_new)
		end
	end
	Printf.@printf "."
	return x_old, f_old, zeros(T), infeas_old
end

#filling opt>0 maximin
function WarmUp(x::Vector{Float64},cons::Vector{Float64},f::Float64,IIFv::Matrix{Float64},TIF::Vector{Float64},findFeasible::Bool,OptFeas::Float64,v::Int64,idxv::Vector{Int64},idxnotv::Vector{Vector{Int64}},IU::IU,FS::FS,Constraints::Constraint,x_forced0::Vector{Int64})
	x_old, f_old, TIF_old, infeas_old=copy(x),Inf,copy(TIF),copy(cons)
	nFS=size(FS.Counts,1)
	nItems=size(IIFv,2)
	T=size(Constraints.ol_max,1)+1
	xv=x[idxv]
	idxv2=setdiff(setdiff(idxv,x_forced0),idxv[xv.==one(Float64) ])
	worstv=one(Float64)
	idxv2=Random.shuffle!(idxv2)
	rows=size(Constraints.catConstrA,1)
	for i2 in idxv2
		x_new, TIF_new, infeas_new=copy(x), copy(TIF_old), copy(infeas_old)
		x_new[i2]=one(Float64)
		xv=x_new[idxv]
		es=0
		cons=0
		if rows>0
			cons=Constraints.catConstrA*xv-Constraints.catConstrb
		end
		ItemUse=(FS.Counts .*evalItemUse(x_new,nFS,T))-IU.Max
		overlap=evalOverlapv(xv,x_new,FS.Counts,Constraints.ol_max,idxv,idxnotv)
		if nItems!=nFS
			x_I=tranformFStoItems(nFS,nItems,xv,Int64[],FS.Items)
		else
			x_I=xv
		end
		if size(Constraints.ExS.Val,1)>0
			es=(Constraints.ExS.Val'*x_I)/sum(x_I)
			es=max(es-Constraints.ExS.Max,-es+Constraints.ExS.Min)
		end
		cons=vcat(cons,overlap,es)
		cons=cons[cons.>0]
		if size(cons,1)>0
			infeas_new[v]=sum(cons)
		else
			infeas_new[v]=0
		end
		TIF_new[v]=evalTIFMMv(x_I,IIFv)
		f_new=(1-OptFeas)*(sum(infeas_new)+max(sum(ItemUse[ItemUse.>0]),0))-OptFeas*minimum(TIF_new)
		if (f_new<f_old)
			x_old, f_old, TIF_old,infeas_old = copy(x_new), copy(f_new), copy(TIF_new), copy(infeas_new)
		end
	end
	Printf.@printf "."
	return x_old, f_old, TIF_old, infeas_old
end
#filling opt>0 CC
function WarmUp(x::Vector{Float64},cons::Vector{Float64},f::Float64,IIFv::Array{Float64,3},TIF::Vector{Float64},α::Float64,findFeasible::Bool,OptFeas::Float64,v::Int64,idxv::Vector{Int64},idxnotv::Vector{Vector{Int64}},IU::IU,FS::FS,Constraints::Constraint,x_forced0::Vector{Int64})
	best_i=0
	x_old, f_old, TIF_old, infeas_old=copy(x),Inf,copy(TIF),copy(cons)
	nFS=size(FS.Counts,1)
	nItems=size(IIFv,2)
	T=size(Constraints.ol_max,1)+1
	xv=x[idxv]
	idxv2=setdiff(setdiff(idxv,x_forced0),idxv[xv.==one(Float64)])
	worstv=1
	idxv2=Random.shuffle!(idxv2)
	for i2 in idxv2
		x_new, TIF_new, infeas_new=copy(x), copy(TIF_old), copy(infeas_old)
		x_new[i2]=one(Float64)
		xv=x_new[idxv]
		rows=size(Constraints.catConstrA,1)
		cons=0
		es=0
		if rows>0
			cons=Constraints.catConstrA*xv-Constraints.catConstrb
		end
		ItemUse=(FS.Counts .*evalItemUse(x_new,nFS,T))-IU.Max
		overlap=evalOverlapv(xv,x_new,FS.Counts,Constraints.ol_max,idxv,idxnotv)
		if nItems!=nFS
			x_I=tranformFStoItems(nFS,nItems,xv,Int64[],FS.Items)
		else
			x_I=xv
		end
		if size(Constraints.ExS.Val,1)>0
			es=(Constraints.ExS.Val'*x_I)/sum(x_I)
			es=max(es-Constraints.ExS.Max,-es+Constraints.ExS.Min)
		end
		cons=vcat(cons,overlap,es)
		cons=cons[cons.>0]
		if size(cons,1)>0
			infeas_new[v]=sum(cons)
		else
			infeas_new[v]=0
		end
		TIF_new[v]=evalTIFCCv(x_I,IIFv;α=α)
		f_new=(1-OptFeas)*(sum(infeas_new)+max(sum(ItemUse[ItemUse.>0]),0))-OptFeas*minimum(TIF_new)
		if (f_new<f_old)
			x_old, f_old, TIF_old,infeas_old = copy(x_new), copy(f_new), copy(TIF_new), copy(infeas_new)
		end
	end
	return x_old, f_old, TIF_old, infeas_old
end

function Optim(ATAmodel::model,x_old,temp_0::Float64,geom_temp::Float64;max_evals=1e6,max_time=1000.00,nItemSample=1000,nTestSample=2,verbosity=2,GroupByFS=false,OptFeas=0.0,nRand=1 ,feasNH=0,optNH=5)
	T=ATAmodel.Settings.T
	nItems=ATAmodel.Settings.nItems
	if ATAmodel.Settings.OptType=="MAXIMIN"
		JLD2.@load "OPT/IIF.jld2" IIF
	elseif ATAmodel.Settings.OptType=="CC"
		JLD2.@load "OPT/IIF_CC.jld2" IIF
		α=ATAmodel.Obj.AuxFloat
	end
	##
	nNH=feasNH+optNH
	ATAmodel.Output.Neighborhoods=[Neighborhood() for n=1:nNH]
	optNH=feasNH+1
	if optNH==1
		findFeasible=false
	else
		findFeasible=true
	end
	nFS=size(ATAmodel.Settings.FS.Items,1)
	startTime=time()
	nacc = 0 # total accepted trials
	t = copy(temp_0) # temperature - will initially rise or fall to cover parameter space. Then it will fall
	converge = 0 # convergence indicator 0 (failure), 1 (normal success), or 2 (convergence but near bounds)
	f_evals=0
	hline = "="^80
	fstar = typemax(Float64)*ones(2)
	NH=1
	ATAmodel.Output.ElapsedTime=time()-startTime
	n=T*nFS
	if size(x_old,1)==0
		x_old=zeros(Float64,n)
	end
	f_old,infeas_old,TIF_old=Inf,1000*ones(T),zeros(T)
	x_opt, f_opt, infeas_opt, TIF_opt=copy(x_old),copy(f_old),copy(infeas_old),copy(TIF_old)
	startTime=time()
	idx=Vector{Vector{Int64}}(undef,T)
	for v in 1:T
		idx[v]=collect(1+(v-1)*nFS:nFS+(v-1)*nFS)
	end
	v=1
	v2=2
	x_new,f_new,TIF_new,infeas_new=copy(x_opt),copy(f_opt),copy(TIF_opt),copy(infeas_opt)
	notv=[setdiff(collect(1:T),t)   for t=1:T]

	while converge==0
		#warm up
		println("warm up starting...")

		round=1
		for round=1:nRand
			w_iter=zeros(T)
			x_new,f_new,TIF_new,infeas_new=copy(x_opt),copy(f_opt),copy(TIF_opt),copy(infeas_opt)
			warmup=fill(true,T)
			while any(warmup) #filling forms
				v=findfirst(warmup.==true)
				Consts=ATAmodel.Constraints[v]
				f_v=((1-OptFeas).*(infeas_new) ).-(OptFeas.*TIF_new)
				mm=f_v[v]
				for v2 in findall(warmup.==true)
					if f_v[v2]>mm
						mm=copy(f_v[v2])
						v=copy(v2)
					end
				end
				idxv=idx[v]
				idxnotv=idx[notv[v]]
				#filling
				if OptFeas==0
					x_new, f_new, TIF_new, infeas_new = WarmUp(x_new,infeas_new,f_new,v,idxv,idxnotv,ATAmodel.IU,ATAmodel.Settings.FS,Consts,ATAmodel.Settings.forced0)
				else
					if ATAmodel.Settings.OptType=="MAXIMIN"
						x_new, f_new, TIF_new, infeas_new = WarmUp(x_new,infeas_new,f_new,IIF[v],TIF_new,true,OptFeas,v,idxv,idxnotv,ATAmodel.IU,ATAmodel.Settings.FS,Consts,ATAmodel.Settings.forced0)
					elseif ATAmodel.Settings.OptType=="CC"
						x_new, f_new, TIF_new, infeas_new = WarmUp(x_new,infeas_new,f_new,IIF[v],TIF_new,α,true,OptFeas,v,idxv,idxnotv,ATAmodel.IU,ATAmodel.Settings.FS,Consts,ATAmodel.Settings.forced0)
					end
				end
				n_t=x_new[idxv]'*ATAmodel.Settings.FS.Counts
				if n_t>=Consts.length_min
					if n_t<=Consts.length_max
						#try to add other items, the first time it goes over n_max it stops
						while n_t<Consts.length_max
							if OptFeas==0
								x_add, f_add, TIF_add, infeas_add = WarmUp(x_new,infeas_new,f_new,v,idxv,idxnotv,ATAmodel.IU,ATAmodel.Settings.FS,Consts,ATAmodel.Settings.forced0)
							else
								if ATAmodel.Settings.OptType=="MAXIMIN"
									x_add, f_add, TIF_add, infeas_add = WarmUp(x_new,infeas_new,f_new,IIF[v],TIF_new,true,OptFeas,v,idxv,idxnotv,ATAmodel.IU,ATAmodel.Settings.FS,Consts,ATAmodel.Settings.forced0)
								elseif ATAmodel.Settings.OptType=="CC"
									x_add, f_add, TIF_add, infeas_add = WarmUp(x_new,infeas_new,f_new,IIF[v],TIF_new,α,true,OptFeas,v,idxv,idxnotv,ATAmodel.IU,ATAmodel.Settings.FS,Consts,ATAmodel.Settings.forced0)
								end
							end
							n_t=x_add[idxv]'*ATAmodel.Settings.FS.Counts
							if n_t<=Consts.length_max
								x_new, f_new, TIF_new, infeas_new=copy(x_add), copy(f_add), copy(TIF_add), copy(infeas_add)
							end
						end
						warmup[v]=false
					else
						w_iter[v]+=1
						x_new[idxv]=zeros(nFS)
						TIF_old[v]=0
						if w_iter[v]>T
							vs=sample(setdiff(collect(1:T),v),1)[1]
							x_new[idx[vs]]=zeros(nFS)
							warmup[vs]=true
							w_iter[vs]=0
							TIF_old[vs]=0
						end
					end
				end
			end
			if f_new<=f_opt
				f_opt, x_opt, TIF_opt, infeas_opt = copy(f_new), copy(x_new), copy(TIF_new), copy(infeas_new)
			end
		end
		if NH>=optNH
			for v=1:T
				xv=x_opt[idx[v]]
				x_I=tranformFStoItems(nFS,nItems,xv,Int64[],ATAmodel.Settings.FS.Items)
				if ATAmodel.Settings.OptType=="MAXIMIN"
					TIF_opt[v]=evalTIFMMv(x_I,IIF[v])
					minTIF=minimum(TIF_new)#,worstv=findmin(TIF_new)
				elseif ATAmodel.Settings.OptType=="CC"
					TIF_opt[v]=evalTIFCCv(x_I,IIF[v];α=α)
					minTIF=minimum(TIF_new)#,worstv=findmin(TIF_new)
				end
			end
		end
		f_opt=sum(infeas_opt)-minimum(TIF_opt)
		f_old, x_old, TIF_old, infeas_old = copy(f_opt), copy(x_opt), copy(TIF_opt), copy(infeas_opt)
		c_iter=0
		println("end of warmup")
		if sum(infeas_old)==0
			println("Feasible solution found in warm up")
		else
			println("Feasible solution not found in warm up")
		end
		c_iter=0
		v2=1
		# statistics to repogeom_temp at each temp change, set back to zero
		coverage_ok=0
		nT=copy(nTestSample)
		nI=copy(nItemSample)
		convergence=0
		notv=[setdiff(collect(1:T),t)   for t=1:T]
		while coverage_ok==0
			nrej=0
			testOrder=sortperm(infeas_old-TIF_old, rev=true)
			if nT>T
				nT=T
			end
			#testOrder=testOrder[1:nT]
			t2=0
			exit=0
			while exit==0 && t2<nT#t2<=(size(testOrder,1)-1)
				t2+=1
				v=testOrder[t2]
				#  for v in testOrder
				TIF_new, infeas_new, x_new =copy(TIF_old), copy(infeas_old),copy(x_old)
				idxv=idx[v]
				idxnotv=idx[notv[v]]
				Consts=ATAmodel.Constraints[v]
				IIFv=IIF[v]
				xv_bin=x_new[idxv]
				xv_int=findall(xv_bin.==one(Float64))
				n_t=sum(xv_bin)
				takenItems=idxv[xv_int]
				if nI>n_t
					nI=Int(n_t)
				end
				#takenItems=sample(takenItems,nI,replace=false)
				idxvav=Random.shuffle!(setdiff(setdiff(idxv,ATAmodel.Settings.forced0),takenItems))
				takenItems=Random.shuffle!(takenItems)
				exit=0
				h2=0
				while exit==0 && h2<nI
					h2+=1
					#TIF_new, infeas_new, x_new =copy(TIF_old), copy(infeas_old),copy(x_old)
					#try to remove h
					#h=sample(takenItems)
					h_bin=takenItems[h2]
					x_new[h_bin]=zero(Float64)
					xv_bin=x_new[idxv]
					xv_int=findall(xv_bin.==one(Float64))
					cons=0
					es=0
					if size(Consts.catConstrA,1)>0
						cons=LinearAlgebra.BLAS.gemv('N',Consts.catConstrA,xv_bin)-Consts.catConstrb#Av*xv-bv
					end
					ItemUse=(ATAmodel.Settings.FS.Counts .*(evalItemUse(x_new,nFS,T)))-ATAmodel.IU.Max
					overlap=evalOverlapv(xv_bin,x_new,ATAmodel.Settings.FS.Counts,Consts.ol_max,idxv,idxnotv)
					if nItems!=nFS
						x_I=tranformFStoItems(nFS,nItems,Float64[],xv_int,ATAmodel.Settings.FS.Items)
					else
						x_I=xv_bin
					end
					if size(Consts.ExS.Val,1)>0
						es=dot(Consts.ExS.Val,x_I)/sum(x_I)
						es=max(es-Consts.ExS.Max,-es+Consts.ExS.Min)
					end
					cons=vcat(cons,overlap,es)
					cons=cons[cons.>0]
					if size(cons,1)>0
						infeas_new[v]=sum(cons)
					else
						infeas_new[v]=0
					end
					f_new=(1-OptFeas)*(sum(infeas_new)+max(sum(ItemUse[ItemUse.>0]),0))
					if findFeasible==false
						if ATAmodel.Settings.OptType=="MAXIMIN"
							TIF_new[v]=evalTIFMMv(x_I,IIFv)
							#TIF_new[v]=evalTIFMMv(x_I,IIFv)
							minTIF=minimum(TIF_new)
						elseif ATAmodel.Settings.OptType=="CC"
							TIF_new[v]=evalTIFCCv(x_I,IIFv;α=α)
							minTIF=minimum(TIF_new)
						end
						f_new=f_new-OptFeas*minTIF
					end
					f_evals+=1
					if (f_new <= f_old)
						#remove item
						x_old, f_old,  TIF_old, infeas_old = copy(x_new), copy(f_new), copy(TIF_new) , copy(infeas_new)
						if (f_new < f_opt)
							println("better with remove, new f_opt= ",f_new)
							exit=1
							convergence=0
							nrej=0
							nT=copy(nTestSample)
							nI=copy(nItemSample)
							x_opt, f_opt,  TIF_opt, infeas_opt = copy(x_new),copy(f_new),copy(TIF_new), copy(infeas_new)
							# nnew +=one(Float64)
						end
					else
						p = exp_c(-(f_new - f_old) / t)
						if (rand() < p)
							#remove
							exit=1
							x_old, f_old, TIF_old, infeas_old= copy(x_new),copy(f_new),copy(TIF_new), copy(infeas_new)
						else
							nrej += 1
						end
					end
					#try to switch h with i2
					x_start=copy(x_new)
					i3=1
					while exit==0 && i3<=size(idxvav,1)
						i2=idxvav[i3]
						if i2!=h_bin
							TIF_new, infeas_new, x_new =copy(TIF_old), copy(infeas_old),copy(x_start)
							x_new[i2]=one(Float64)
							xv_bin=x_new[idxv]
							xv_int=findall(xv_bin.==one(Float64))
							es=0
							cons=0
							if size(Consts.catConstrA,1)>0
								cons=LinearAlgebra.BLAS.gemv('N', Consts.catConstrA, xv_bin)-Consts.catConstrb# phi=New_pars*X1', if A'*B then 'T', 'N'
							end
							iu=(ATAmodel.Settings.FS.Counts .* (evalItemUse(x_new,nFS,T)))-ATAmodel.IU.Max
							ol=evalOverlapv(xv_bin,x_new,ATAmodel.Settings.FS.Counts,Consts.ol_max,idxv,idxnotv)
							if nItems>nFS
								x_I=tranformFStoItems(nFS,nItems,Float64[],xv_int,ATAmodel.Settings.FS.Items)
							else
								x_I=xv_bin
							end
							if size(Consts.ExS.Val,1)>0
								es=LinearAlgebra.dot(Consts.ExS.Val, x_I)/sum(x_I)#(ESv'*x_I)/sum(x_I)
								es=max(es-Consts.ExS.Max,-es+Consts.ExS.Min)
							end
							cons=vcat(cons,ol,es)
							cons=cons[cons.>0]
							if size(cons,1)>0
								infeas_new[v]=sum(cons)
							else
								infeas_new[v]=0
							end
							f_new=(1-OptFeas)*(sum(infeas_new)+max(sum(iu[iu.>0]),0))
							if findFeasible==false
								if ATAmodel.Settings.OptType=="MAXIMIN"
									TIF_new[v]=evalTIFMMv(x_I,IIFv)
									minTIF=minimum(TIF_new)
								elseif ATAmodel.Settings.OptType=="CC"
									TIF_new[v]=evalTIFCCv(x_I,IIFv;α=α)
									minTIF=minimum(TIF_new)
								end
								f_new=f_new-OptFeas*minTIF
							end
							f_evals+=1
							if (f_new <= f_old)
								#println("switch with ",i2)
								x_old, f_old, TIF_old, infeas_old = copy(x_new), copy(f_new), copy(TIF_new) , copy(infeas_new)
								if (f_new < f_opt)
									println("better with switch, new f_opt= ",f_new)
									exit=1
									convergence=0
									nrej=0
									nT=copy(nTestSample)
									nI=copy(nItemSample)
									x_opt, f_opt,  TIF_opt, infeas_opt = copy(x_new),copy(f_new),copy(TIF_new), copy(infeas_new)
									#nnew +=one(Float64)
								end
							else
								#remove if p>rand()
								p = exp_c(-(f_new - f_old) / t)
								if (rand() < p)
									exit=1
									x_old, f_old, TIF_old, infeas_old = copy(x_new), copy(f_new), copy(TIF_new) , copy(infeas_new)
								else
									nrej += 1
								end
							end
						end #end of i2!=h
						i3+=1
					end #end of i2
					h2+=1
				end #end of exit==0 (h)
			end #end of testorder
			if verbosity > 1
				println(hline)
				println("Intermediate results before next temperature change")
				println("temperature: ", t)
				println("current local Infeasibility: ", sum(infeas_opt))
				println("current local Optimality: ", minimum(TIF_opt))
				println("total evaluations so far: ", f_evals)
				println("total time elapsed: ", ATAmodel.Output.ElapsedTime)
				# println("total moves since last temperature reduction: ", nacc+nrej)#nup + ndown + nrej)
				# println("downhill: ", ndown)
				# println("accepted uphill: ", nupacc)
				# println("rejected uphill: ", nrej)
				# println("new minima this temperature: ", nnew)
				println(hline*"\n") #print results
			end
			fstar[1] = f_old #it was f_old
			if fstar[2]==fstar[1]
				#nI+=one(Float64)
				#nT+=one(Float64)
				convergence+=1
			end
			println("convergence is ", convergence)
			#how many equal f_old in the last iterations?
			if fstar[2]==fstar[1] && convergence==100
				println("neighbourhood ",NH," fully explored, increase temperature")
				nT=copy(nTestSample)
				nI=copy(nItemSample)
				coverage_ok=1
				t = copy(temp_0)
				ATAmodel.Output.Neighborhoods[NH].f=copy(f_opt)
				ATAmodel.Output.Neighborhoods[NH].Design=copy(x_opt)
				ATAmodel.Output.Neighborhoods[NH].Feas=copy(infeas_opt)
				ATAmodel.Output.Neighborhoods[NH].TIF=copy(TIF_opt)
				fstar = typemax(Float64)*ones(2)
				NH+=1
				if NH>nNH
					converge=1
				else
					if NH>=optNH
						findFeasible=false
					end
					#perturbate the solution
					j=1
					for i=1:n
						if x_opt[i]==one(Float64)
							if j==1
								x_opt[i]=0
								j+=1
							else
								j-=1
							end
						end
					end
				end
				for v=1:T
					Consts=ATAmodel.Constraints[v]
					xv_bin=x_opt[idx[v]]
					cons=0
					es=0
					if size(Consts.catConstrA,1)>0
						cons=Consts.catConstrA*xv_bin-Consts.catConstrb
					end
					ItemUse=(ATAmodel.Settings.FS.Counts .*(evalItemUse(x_opt,nFS,T)))-ATAmodel.IU.Max
					overlap=evalOverlapv(xv_bin,x_opt,ATAmodel.Settings.FS.Counts,Consts.ol_max,idx[v],idx[setdiff(collect(1:T),v)])
					if nItems!=nFS
						x_I=tranformFStoItems(nFS,nItems,xv_bin,Int64[],ATAmodel.Settings.FS.Items)
					else
						x_I=xv_bin
					end
					if size(Consts.ExS.Val,1)>0
						es=(Consts.ExS.Val'*x_I)/sum(x_I)
						es=max(es-Consts.ExS.Max,-es+Consts.ExS.Min)
					end
					cons=vcat(cons,overlap,es)
					cons=cons[cons.>0]
					if size(cons,1)>0
						infeas_opt[v]=sum(cons)
					else
						infeas_opt[v]=0
					end
					f_opt=sum(infeas_opt)+max(sum(ItemUse[ItemUse.>0]),0)
					if findFeasible==false
						if ATAmodel.Settings.OptType=="MAXIMIN"
							TIF_opt[v]=evalTIFMMv(x_I,IIF[v])
							minTIF,worstv=findmin(TIF_opt)
						elseif ATAmodel.Settings.OptType=="CC"
							TIF_opt[v]=evalTIFCCv(x_I,IIF[v];α=α)
							minTIF,worstv=findmin(TIF_opt)
						end
						f_opt=(1-OptFeas)*f_opt-OptFeas*minTIF
					end
				end
			else
				println("neighbourhood ",NH," not entirely explored, decrease temperature")
				coverage_ok=0
				if fstar[1]<=fstar[2]
					fstar[2]=fstar[1]
				end
			end
			if f_evals>=max_evals
				converge=2
			end
			ATAmodel.Output.ElapsedTime=time()-startTime
			if ATAmodel.Output.ElapsedTime>=max_time
				converge=3
			end
			if coverage_ok==0
				#change neighbourhood
				# Reduce temperature, record current function value in the
				t *= geom_temp
				x_old, f_old,  TIF_old, infeas_old = copy(x_opt), copy(f_opt), copy(TIF_opt) , copy(infeas_opt)
			else  # coverage not ok - increase temperature quickly to exp_cand search area
				t += 1
				x_old, f_old,  TIF_old, infeas_old = copy(x_opt), copy(f_opt), copy(TIF_opt) , copy(infeas_opt)
			end
		end #end of NH coverage
		if converge>0
				f_opt,NH=findmin([ATAmodel.Output.Neighborhoods[n].f   for n=1:nNH])
				x_opt=ATAmodel.Output.Neighborhoods[NH].Design
				infeas_opt=ATAmodel.Output.Neighborhoods[NH].Feas
				TIF_opt=ATAmodel.Output.Neighborhoods[NH].TIF
				ATAmodel.Output.Design=x_opt
				ATAmodel.Output.Feas=infeas_opt
				ATAmodel.Output.f=f_opt
				ATAmodel.Output.TIF=TIF_opt
				if verbosity >1
					println(hline)
					println("Results")
					if (converge == 1)
						println("==> Maximum number of neighborhoods explored <==")
					elseif (converge==2) #time max reached
						println("==> Maximum number of evaluations reached <==")
					elseif (converge==3) #evals max reached
						println("==> Maximum time reached <==")
					end
					Printf.@printf("\n     Obj. value:	%16.10f", f_opt)
					Printf.@printf("\n     Optimality:	%16.10f", minimum(TIF_opt))
					Printf.@printf("\n     Infeasibility:	%16.10f", sum(infeas_opt))
					Printf.@printf("\n     Obj. function evaluations:	%16.10f", f_evals)
					Printf.@printf("\n     Elapsed Time:	%16.10f", ATAmodel.Output.ElapsedTime)
					println("\n")
					println(hline)
				end
		end
	end #end of converge


	if !("RESULTS" in readdir())
		mkdir("RESULTS")
	end
	JLD2.@save "RESULTS/ATAmodel.jld2" ATAmodel
	writedlm("RESULTS/design.csv",reshape(x_opt,nFS,T))
	open("RESULTS/ResultsATA.txt", "w") do io
		write(io,"tests")
		write(io,"\r\n")
		writedlm(io, collect(1:T)',",")
		write(io,"\r\n")
		write(io,"f_opt")
		write(io,"\r\n")
		writedlm(io, f_opt)
		write(io,"\r\n")
		write(io,"infeasibility")
		write(io,"\r\n")
		writedlm(io, infeas_opt',",")
		write(io,"\r\n")
		write(io,"TIF")
		write(io,"\r\n")
		writedlm(io, TIF_opt',",")
		write(io,"\r\n")
		write(io,"elapsed Time   for optimization")
		write(io,"\r\n")
		write(io, string(ATAmodel.Output.ElapsedTime))
		write(io,"\r\n")
	end
	return ATAmodel
end


function PrintResults(ATAmodel; GroupByFS=false, BSfolder="BS")
	Plots.pgfplots()
	if size(ATAmodel.Output.Design,1)>0
		IRTpars=ATAmodel.Settings.IRTparameters
		T=ATAmodel.Settings.T
		nItems=ATAmodel.Settings.nItems
		length_max=[ATAmodel.Constraints[t].length_max   for t=1:T]
		CATEGORIES=ATAmodel.Output.Categories
		#summarize=Main.summarize
		model=ATAmodel.Settings.IRTmodel.model
		nFS=size(ATAmodel.Settings.FS.Items,1)
		if nFS<ATAmodel.Settings.nItems
			n=nFS*T
		else
			n=nItems*T
		end

		if ATAmodel.Settings.OptType=="CC"
			JLD2.@load "OPT/IIF_CC.jld2" IIF
			JLD2.@load "OPT/ICF_CC.jld2" ICF
			α=ATAmodel.Obj.AuxFloat
			IIF_CC=copy(IIF)
			ICF_CC=copy(ICF)
		end
		if ATAmodel.Settings.OptType=="MAXIMIN" || ATAmodel.Settings.OptType=="CC"
			JLD2.@load "OPT/IIF.jld2" IIF
			JLD2.@load "OPT/ICF.jld2" ICF
		end
		design=reshape(ATAmodel.Output.Design,nFS,T)
		if GroupByFS==true
			new_design=zeros(Float64,nItems,T)
			for t=1:T
				new_design[vcat(ATAmodel.Settings.FS.Items[findall(design[:,t].==1.0)]...),t].=1.0
			end
			design=new_design
		end
		#overlap
		olMatrixOUT=zeros(Int64,T,T)
		for t=1:T
			for v=1:T
				olMatrixOUT[t,v]=design[:,t]'*design[:,v]
			end
		end
		writedlm("RESULTS/olMatrixOUT.csv",olMatrixOUT)
		#TIF e ICF
		if ATAmodel.Settings.OptType=="MAXIMIN" ||  ATAmodel.Settings.OptType=="CC"
			ThetasPlot=collect(range(-4,stop=4,length=101)) #nqp values in interval/r/n",
			IIFplot=ItemInfoFun(IRTpars,ThetasPlot,model=ATAmodel.Settings.IRTmodel.model)
			ICFplot=ItemCharFun(IRTpars,ThetasPlot,model=ATAmodel.Settings.IRTmodel.model)
			IIFf=Array{Float64,2}(undef,T,101)
			ICFf=Array{Float64,2}(undef,T,101)
			for t in 1:T
				IIFf[t,:]=[sum(IIFplot[findall(design[:,t].==1 ),i])   for i in 1:101]
				ICFf[t,:]=[sum(ICFplot[findall(design[:,t].==1 ),i])   for i in 1:101]
				#IIFf[t,:]=vec(colwise(sum,DataFrame(IIFplot[ATADesign[:,t].==1 ,:])))
			end
			Plots.plot(ThetasPlot,IIFf',xlims=(-4,4),xticks = -4:1:4,
			titlefontsize=16, size=(500,400),yticks=0:2:maximum(IIFf)+2,ylims=(0,maximum(IIFf)+2),tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.7,foreground_color_border=:black, linewidth=0.5, label=[string("t",t)   for t=1:T],
			);
			Plots.savefig("RESULTS/TIFPlot.pdf")
			Plots.plot(ThetasPlot,ICFf',xlims=(-4,4), xticks = -4:1:4,
			titlefontsize=16, size=(500,400),tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.7,foreground_color_border=:black, linewidth=0.5,label=[string("t",t)   for t=1:T],
			);
			Plots.savefig("RESULTS/ICFPlot.pdf")
			if isfile("simPool.csv")
				simPool=CSV.read("simPool.csv")
				IIFtrue=Vector{Vector{Float64}}(undef,T)
				for t=1:T
					IIFtrue[t]=ItemInfoFun(simPool,ThetasPlot,model=ATAmodel.Settings.IRTmodel.model)'*design[:,t]
				end
			end
			if ATAmodel.Settings.OptType=="CC"
				IIF_plot=Vector{Array{Float64,3}}(undef,T)
				ICF_CC_plot=Vector{Array{Float64,3}}(undef,T)
				ThetasPlot=collect(range(-4,stop=4,length=101)) #nqp values in interval/r/n",
				JLD2.@load "BSPar.jld2" BSPar
				BSa=Matrix(BSPar[2])[:,2:end]
				BSb=Matrix(BSPar[1])[:,2:end]
				R=ATAmodel.Obj.AuxInt
				for t=1:T
					println(t)
					IIF_plot[t]=zeros(101,ATAmodel.Settings.nItems,R)
					ICF_CC_plot[t]=zeros(101,ATAmodel.Settings.nItems,R)
					for r=1:R
						if model=="1PL"
							df=DataFrame(b=BSb[:,r]) #nqp values in interval\r\n",
						elseif model=="2PL"
							df=DataFrame(a=BSa[:,r],b=BSb[:,r]) #nqp values in interval\r\n",
						elseif model=="3PL"
							df=DataFrame(a=BSa[:,r],b=BSb[:,r],c=BSc[:,r])
						end
						for k=1:101
							IIF_plot[t][k,:,r]=ItemInfoFun(df,ThetasPlot[k],model=ATAmodel.Settings.IRTmodel.model) # IxK[t]
							ICF_CC_plot[t][k,:,r]=ItemCharFun(df,ThetasPlot[k],model=ATAmodel.Settings.IRTmodel.model)# IxK[t]
						end
					end
				end
				#TIF=Array{Array{Float64,2},1}(undef,T)
				IIFdesigntoplot=Array{Array{Float64,2},1}(undef,T)
				for t in 1:T
					TIF=Array{Float64,2}(undef,101,R)
					IIFdesigntoplot[t]=Array{Float64,2}(undef,101,6)
					for k=1:101
						for r in 1:R
							TIF[k,r]=(IIF_plot[t][k,:,r]'*design[:,t])#,[0,0.25,0.5,0.75,1,α])[1:6]
						end
						IIFdesigntoplot[t][k,:]=quantile(TIF[k,:],[0,0.25,0.5,0.75,1,α])
					end
					CSV.write(string("RESULTS/IIFdesigntoplot_",t,".csv"),DataFrame(IIFdesigntoplot[t]))
				end

				for t in 1:T
					#plot(IIFdesigntoplot[t][:,1],IIFdesigntoplot[t][:,2],seriestype=:scatter)
					Plots.plot(size=(500,400),yticks=0.0:2.0:maximum(IIFf)+2,ylims=(0,maximum(IIFf)+2),xlims=(-4.0,4.0));
					if isfile("simPool.csv")
						Plots.plot!(ThetasPlot,IIFtrue[t],size=(500,400),yticks=0:2:maximum(IIFf)+2,ylims=(0,maximum(IIFf)+2),tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.7,foreground_color_border=:black, linewidth=0.5,linestyle = :solid,linecolor=:violetred4,label="True");
					end
					Plots.plot!(ThetasPlot,IIFdesigntoplot[t][:,5],
					titlefontsize=16, size=(500,400),tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.7,foreground_color_border=:black, linewidth=0.5,
					linestyle = :dot,linecolor=:darkcyan,label="Max");
					Plots.plot!(ThetasPlot,IIFdesigntoplot[t][:,4],tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.7,foreground_color_border=:black, linewidth=0.5,yticks=0:2:maximum(IIFf)+2,ylims=(0,maximum(IIFf)+2),linestyle = :dash,linecolor=:darkcyan,label="75-Qle");
					Plots.plot!(ThetasPlot,IIFdesigntoplot[t][:,3],tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.7,foreground_color_border=:black, linewidth=0.5,linestyle = :solid,linecolor=:darkcyan,label="Median");
					Plots.plot!(ThetasPlot,IIFdesigntoplot[t][:,2],tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.7,foreground_color_border=:black, linewidth=0.5,linestyle = :dash,linecolor=:darkcyan,label="25-Qle");
					Plots.plot!(ThetasPlot,IIFdesigntoplot[t][:,6],tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.7,foreground_color_border=:black, linewidth=0.5,linestyle = :dashdotdot,linecolor=:indigo,label=L"{\alpha}-Qle");
					Plots.plot!(ThetasPlot,IIFdesigntoplot[t][:,1],tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.7,foreground_color_border=:black, linewidth=0.5,linestyle = :dot,linecolor=:darkcyan,label="Min");

					Plots.plot!(ThetasPlot,IIFf[t,:],tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.7,foreground_color_border=:black, linewidth=0.5,colour=[:black],label="estimated")
					Plots.savefig(string("RESULTS/",t,"_infoplot.pdf"))
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
			for t=1:T
				K=size(ATAmodel.Obj.OptPts[t],1)
				for k=1:K
					IIF_t=(IIF_CC[t][k,:,:]'*design[:,t])
					Min=vcat(Min,minimum(IIF_t))
					First=vcat(First,StatsBase.quantile(IIF_t,0.25))
					Median=vcat(Median,StatsBase.median(IIF_t))
					Third=vcat(Third,StatsBase.quantile(IIF_t,0.75))
					Max=vcat(Max,maximum(IIF_t))
					Alpha=vcat(Alpha,StatsBase.quantile(IIF_t,α))
					Estimated=vcat(Estimated,IIF[t][k,:]'*design[:,t])
					if isfile("simPool.csv")
						True=vcat(True,ItemInfoFun(simPool,ATAmodel.Obj.OptPts[t][k],model=ATAmodel.Settings.IRTmodel.model)'*design[:,t])
					end
				end
			end
			TIFatTheta_k=DataFrame(hcat(Min,First,Median,Third,Max,Alpha,Estimated))

			if isfile("simPool.csv")
				TIFatTheta_k[:True]=True
				DataFrames.names!(TIFatTheta_k,Symbol.(["Min", "25-Qle", "Median", "75-Qle", "Max" ,L"{\alpha}-Qle","Estimated","True"]))
			else
				DataFrames.names!(TIFatTheta_k,Symbol.(["Min", "25-Qle", "Median", "75-Qle", "Max" ,L"{\alpha}-Qle","Estimated"]))
			end
			CSV.write("RESULTS/TIFatTheta_k.csv", TIFatTheta_k)

		else
			Estimated=Float64[]
			if isfile("simPool.csv")
				True=Float64[]
			end
			for t=1:T
				K=size(ATAmodel.Obj.OptPts[t],1)
				for k=1:K
					Estimated=vcat(Estimated,IIF[t][k,:]'*design[:,t])
					if isfile("simPool.csv")
						True=vcat(True,ItemInfoFun(simPool,ATAmodel.Obj.OptPts[t][k],model=ATAmodel.Settings.IRTmodel.model)'*design[:,t])
					end
				end
			end
			TIFatTheta_k=DataFrame(Estimated=Estimated)
			if isfile("simPool.csv")
				TIFatTheta_k[:True]=True
				DataFrames.names!(TIFatTheta_k,Symbol.(["Estimated","True"]))
			end
			CSV.write("RESULTS/TIFatTheta_k.csv", TIFatTheta_k)

		end

		#expected score
		ESprintIRT=Array{Float64,1}(undef,T)
		for t in 1:T
			ESprintIRT[t]=((ICF[t]*design[:,t]))[1]/sum(design[:,t])
		end
		#number of items
		n=Array{Float64,1}(undef,T)
		for t in 1:T
			n[t]=sum(design[:,t])
		end
		designItems=Matrix{Int64}(undef,Int(maximum(n)),T)
		#designItems=Array{String,2}(nothing,n_max,T)
		for t in 1:T
			designItems[collect(1:Int(n[t])),t].=sort!(findall(design[:,t].==1))#sort!(CODEVar[findall(design[:,t].==1)])
		end
		open("RESULTS/Results.txt", "w") do io
			write(io,"Length")
			write(io,"\r\n")
			writedlm(io, [1:T,n'],",")
			write(io, "\r\n")
			write(io, "Design")
			write(io,"\r\n")
			writedlm(io, designItems,",")
			write(io,"\r\n")
			write(io, "expected score")
			write(io,"\r\n")
			writedlm(io, [1:T,ESprintIRT'],",")
			if size(CATEGORIES,1)>0
				write(io,"Categorical variables")
				write(io,"\r\n")
				for cats in CATEGORIES
					write(io,cats)
					write(io,"\r\n")
					vals=setdiff(unique(ATAmodel.Settings.Bank[Symbol(cats)]),[missing])
					println(vals)
					writedlm(io,hcat(vals,[sum(skipmissing(ATAmodel.Settings.Bank[Symbol(cats)][design[:,t].==1].== i))   for i in vals,t=1:T]),", ")
					write(io,"\r\n")
				end
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
			ItemUse=[sum(design[i,:])   for i=1:ATAmodel.Settings.nItems]
			ids=sortperm(ItemUse, rev=true)
			sort!(ItemUse, rev=true)
			for i=1:ATAmodel.Settings.nItems
				write(io,string(ids[i]))
				write(io,",")
				write(io,string(ItemUse[i]))
				write(io,"\r\n")
			end
		end
	end
end
