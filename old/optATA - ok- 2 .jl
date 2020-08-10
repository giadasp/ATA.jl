function checkFeas(infeas_new::Float64,FS::FS,Cons_v::Constraint,x::Matrix{Float64},x_v::Vector{Float64},nFS::Int64,nItems::Int64,v::Int64)
	es=0
	cons=0
	if size(Cons_v.catConstrA,1)>0
		cons=copy(Cons_v.catConstrb)
		#LinearAlgebra.BLAS.gemv('N',Cons_v.catConstrA,x_v)-Cons_v.catConstrb
		cons=gemvblas(Cons_v.catConstrA,x_v,cons,size(x_v,1))
	end
	#ol=evalOverlapv(x_v,x,FS.Counts,Cons_v.ol_max,v)#idx_v,idx_notv)
	if nItems>nFS
		x_I_v=tranformFStoItems(nFS,nItems,x_v,FS.Items)
	else
		x_I_v=x_v
	end
	if size(Cons_v.ExS.Val,1)>0
		es=dot(Cons_v.ExS.Val,x_I_v)/sum(x_I_v)
		es=max(es-Cons_v.ExS.Max,-es+Cons_v.ExS.Min)
	end
	cons=vcat(cons,es)#,0.5 .*ol)
	cons=cons[cons.>0]
	#println("cons: ",cons)
	if size(cons,1)>0
		infeas_new=sum(cons)
	else
		infeas_new=zero(Float64)
	end
	return infeas_new::Float64, x_I_v::Vector{Float64}
end
function tranformFStoItems(nFS::Int64,nItems::Int64,x_v::Vector{Float64},FSitems::Vector{Vector{Int64}})
	x_v_taken=findall(x_v.==one(Float64))
	x_I_v=zeros(Float64,nItems)
	for i in x_v_taken
		x_I_v[FSitems[i]].=one(Float64)
	end
	return x_I_v
end
#filling OptFeas==0
function WarmUp(x::Matrix{Float64},cons::Vector{Float64},v::Int64,IU::IU,FS::FS,Constraints::Constraint,x_forced0_v::Vector{Bool},nItems::Int64,nFS::Int64,FScounts::Matrix{Float64},olMax::Matrix{Float64},ol_old::Vector{Float64})
	f_old,  infeas_old =  Inf, copy(cons)
	T=size(Constraints.ol_max,1)+1
	x_v=x[:,v]
	idx_v2=setdiff(collect(1:nFS),findall(x_v.==one(Float64))) #i2=1,...,I
	idx_v2=Random.shuffle!(idx_v2)
	iu_old=sum(x,dims=2)-IU.Max
	iu_opt=copy(iu_old)
	i_best=1
	ol_old_v=0
	for i2 in idx_v2
		if x_forced0_v[i2]
			x_new, infeas_new, iu=copy(x), copy(infeas_old), copy(iu_old)
			x_new[i2,v]=one(Float64)
			iu[i2]=iu[i2]+1
			iu=iu[iu.>0]
			if size(iu,1)==0
				iu=0
			else
				iu=sum(iu)
			end
			x_v=copy(x_new[:,v])
			infeas_new[v], x_I_v=checkFeas(infeas_new[v],FS,Constraints,x_new,x_v,nFS,nItems,v)
			#ol=evalOverlap(x_new,FScounts,olMax,T)
			# for v=1:T
			# 	ol_v=ol[setdiff(collect(1:T),v),v]
			# 	ol_new[v]=sum(ol_v[ol_v.>0])
			# end
			ol=evalOverlapv(x_v,x_new,FS.Counts,olMax[:,v],v)
			f_new=sum(infeas_new)+iu+ol
			if (f_new<f_old)
				f_old, infeas_old ,  iu_opt =  copy(f_new),  copy(infeas_new),  copy(iu)
				ol_old_v=copy(ol)
				i_best=copy(i2)
			end
		end
	end
	ol_old[v]=copy(ol_old_v)
	x[i_best,v]=one(Float64)
	return x, f_old, infeas_old, ol_old, iu_opt
end

#filling opt>0 maximin
function WarmUp(x::Matrix{Float64},cons::Vector{Float64},IIFv::Matrix{Float64},TIF::Vector{Float64},OptFeas::Float64,v::Int64,IU::IU,FS::FS,Constraints::Constraint,x_forced0_v::Vector{Bool},nItems::Int64,nFS::Int64,FScounts::Matrix{Float64},olMax::Matrix{Float64},ol_old::Vector{Float64})
	f_old, TIF_old, infeas_old=Inf,copy(TIF),copy(cons)
	T=size(Constraints.ol_max,1)+1
	x_v=x[:,v]
	idx_v2=setdiff(collect(1:nFS),findall(x_v.==one(Float64))) #i2=1,...,I
	idx_v2=Random.shuffle!(idx_v2)
	iu_old=sum(x,dims=2)-IU.Max
	iu_opt=copy(iu_old)
	i_best=1
	ol_old_v=0
	for i2 in idx_v2
		if x_forced0_v[i2]
			x_new, TIF_new, infeas_new, iu=copy(x), copy(TIF_old), copy(infeas_old), copy(iu_old)
			x_new[i2,v]=one(Float64)
			x_v=copy(x_new[:,v])
			iu[i2]= iu[i2]+1
			iu=iu[iu.>0]
			if size(iu,1)==0
				iu=0
			else
				iu=sum(iu)
			end
			infeas_new[v], x_I_v=checkFeas(infeas_new[v],FS,Constraints,x_new,x_v,nFS,nItems,v)
			TIF_new[v]=evalTIFMMv(x_I_v,IIFv)
			# ol=evalOverlap(x_new,FScounts,olMax,T)
			# for v=1:T
			# 	ol_v=ol[Int.(setdiff(collect(1:T),v)),v]
			# 	ol_new[v]=sum(ol_v[ol_v.>0])
			# end
			ol=evalOverlapv(x_v,x_new,FS.Counts,olMax[:,v],v)
			f_new=(1-OptFeas)*(sum(infeas_new)+iu+ol)-OptFeas*minimum(TIF_new)
			if (f_new<f_old)
				i_best=copy(i2)
				f_old, TIF_old, infeas_old, iu_opt =  copy(f_new), copy(TIF_new), copy(infeas_new),  copy(iu)
				ol_old_v = copy(ol)
			end
		end
	end
	ol_old[v]=copy(ol_old_v)
	x[i_best,v]=one(Float64)
	return x, f_old, TIF_old, infeas_old, ol_old, iu_opt
end
#filling opt>0 CC
function WarmUp(x::Matrix{Float64},cons::Vector{Float64},IIFv::Array{Float64,3},TIF::Vector{Float64},α::Float64,OptFeas::Float64,v::Int64,IU::IU,FS::FS,Constraints::Constraint,x_forced0_v::Vector{Bool},nItems::Int64,nFS::Int64,FScounts::Matrix{Float64},olMax::Matrix{Float64},ol_old::Vector{Float64})
	f_old, TIF_old, infeas_old=Inf,copy(TIF),copy(cons)
	T=size(Constraints.ol_max,1)+1
	x_v=x[:,v]
	idx_v2=setdiff(collect(1:nFS),findall(x_v.==one(Float64))) #i2=1,...,I
	idx_v2=Random.shuffle!(idx_v2)
	iu_old=sum(x,dims=2)-IU.Max
	iu_opt=copy(iu_old)
	i_best=1
	ol_old_v=0
	for i2 in idx_v2
		if x_forced0_v[i2]
			x_new, TIF_new, infeas_new, iu=copy(x), copy(TIF_old), copy(infeas_old), copy(iu_old)
			x_new[i2,v]=one(Float64)
			x_v=copy(x_new[:,v])
			iu[i2]= iu[i2]+1
			iu=iu[iu.>0]
			if size(iu,1)==0
				iu=0
			else
				iu=sum(iu)
			end
			infeas_new[v], x_I_v=checkFeas(infeas_new[v],FS,Constraints,x_new,x_v,nFS,nItems,v)
			TIF_new[v]=evalTIFCCv(x_I_v,IIFv;α=α)
			# ol=evalOverlap(x_new,FScounts,olMax,T)
			# for v=1:T
			# 	ol_v=ol[Int.(setdiff(collect(1:T),v)),v]
			# 	ol_new[v]=sum(ol_v[ol_v.>0])
			# end
			ol=evalOverlapv(x_v,x_new,FS.Counts,olMax[:,v],v)
			f_new=(1-OptFeas)*(sum(infeas_new)+iu+ol)-OptFeas*minimum(TIF_new)
			if (f_new<f_old)
				i_best=copy(i2)
				f_old, TIF_old, infeas_old,  iu_opt =  copy(f_new), copy(TIF_new), copy(infeas_new),  copy(iu)
				ol_old_v=copy(ol)
			end
		end
	end
	ol_old[v]=copy(ol_old_v)
	x[i_best,v]=one(Float64)
	return x, f_old, TIF_old, infeas_old, ol_old, iu_opt

end

function Optim(ATAmodel::model,x_old::Matrix{Float64},temp_0::Float64,geom_temp::Float64;max_evals=1e6,max_time=1000.00,nItemSample=1000,nTestSample=2,convMax=1,verbosity=2,GroupByFS=false,OptFeas=0.0,nRand=1 ,feasNH=0,optNH=5)
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
	if feasNH==0
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
		x_old=zeros(Float64,nFS,T)
	end
	f_old,infeas_old, ol_old,TIF_old=Inf,1e6*ones(T),zeros(T),zeros(T)
	x_opt, f_opt, infeas_opt, ol_opt, TIF_opt=copy(x_old),copy(f_old),copy(infeas_old), copy(ol_old), copy(TIF_old)
	startTime=time()
	FScounts=ATAmodel.Settings.FS.Counts*ones(Float64,T)'
	x_new,f_new,TIF_new, infeas_new, ol_new=copy(x_opt),copy(f_opt),copy(TIF_opt),copy(infeas_opt), copy(ol_opt)
	olMax=Float64.(readdlm("overlapMatrix.csv",';'))
	iu=copy(ATAmodel.IU.Max)
	while converge==0
		#warm up
		println("warm up starting")
		round=1
		for round=1:nRand
			x_new, f_new, TIF_new, infeas_new, ol_new=copy(x_opt), copy(f_opt), copy(TIF_opt), copy(infeas_opt), copy(ol_opt)
			warmup=fill(true,T)
			while any(warmup) #filling forms
				ol=evalOverlap(x_new,FScounts,olMax,ATAmodel.Settings.T)
				for v=1:ATAmodel.Settings.T
					ol_v=ol[setdiff(collect(1:T),v),v]
					ol_v=ol_v[ol_v.>0]
					if size(ol_v,1)==0
						ol_new[v]=0
					else
						ol_new[v]=sum(ol_v)
					end
				end
				v=findfirst(warmup.==true)
				Constraints=ATAmodel.Constraints[v]
				f_v=((1-OptFeas).*(infeas_new+ol_new))-(OptFeas.*TIF_new)
				mm=f_v[v]
				for v2 in findall(warmup.==true)
					if f_v[v2]>mm
						mm=copy(f_v[v2])
						v=copy(v2)
					end
				end
				#filling
				n_t=dot(x_new[:,v],ATAmodel.Settings.FS.Counts)
				#try to add other items, the first time it goes over n_max it stops
				TIF_add=zeros(ATAmodel.Settings.T)
				if n_t<Constraints.length_max
					if OptFeas==0 || findFeasible==true
						x_add, f_add, infeas_add, ol_add, iu = WarmUp(x_new,infeas_new,v,ATAmodel.IU,ATAmodel.Settings.FS,Constraints,ATAmodel.Settings.forced0[v],nItems,nFS,FScounts,olMax,ol_new)
					else
						if ATAmodel.Settings.OptType=="MAXIMIN"
							x_add, f_add, TIF_add, infeas_add, ol_add, iu = WarmUp(x_new,infeas_new,IIF[v],TIF_new,OptFeas,v,ATAmodel.IU,ATAmodel.Settings.FS,Constraints,ATAmodel.Settings.forced0[v],nItems,nFS,FScounts,olMax,ol_new)
						elseif ATAmodel.Settings.OptType=="CC"
							x_add, f_add, TIF_add, infeas_add, ol_add, iu = WarmUp(x_new,infeas_new,IIF[v],TIF_new,α,OptFeas,v,ATAmodel.IU,ATAmodel.Settings.FS,Constraints,ATAmodel.Settings.forced0[v],nItems,nFS,FScounts,olMax,ol_new)
						end
					end
					Printf.@printf "."
					n_v=dot(x_add[:,v],ATAmodel.Settings.FS.Counts)
					#println("length for test ",v,": ",n_v)
					if n_t<=Constraints.length_max
						x_new, f_new, TIF_new, infeas_new, ol_new=copy(x_add), copy(f_add), copy(TIF_add), copy(infeas_add), copy(ol_add)
					else
						#x_new, f_new, TIF_new, infeas_new=copy(x_add), copy(f_add), copy(TIF_add), copy(infeas_add)
						warmup[v]=false
						println("f_new:", f_new)
						println("-Test ",v," filled up with ",n_t," items,")
					end
				else
					warmup[v]=false
					println("-Test ",v," filled up with ",n_t," items,")
				end
			end
		end
		x_opt, TIF_opt, infeas_opt, ol_opt = copy(x_new), copy(TIF_new), copy(infeas_new), copy(ol_new)
		f_opt=(1-OptFeas)*(sum(infeas_opt+ol_opt)+iu)-OptFeas*minimum(TIF_opt)
		println("ol: ",ol_opt)
		println("iu: ", iu)
		f_old, x_old, TIF_old, infeas_old, ol_old = copy(f_opt), copy(x_opt), copy(TIF_opt), copy(infeas_opt), copy(ol_opt)
		println("end of warmup")
		if sum(infeas_old)==0
			println("Feasible solution found in warm up")
		else
			println("Feasible solution not found in warm up")
		end
		# statistics to repogeom_temp at each temp change, set back to zero
		coverage_ok=0
		if nTestSample>T
			nTestSample=T
		end
		nI=copy(nItemSample)
		convergence=0
		while coverage_ok==0
			testOrder=sortperm(((1-OptFeas).*(infeas_old+ol_old))-((OptFeas).*TIF_old), rev=true)
			t2=0
			exit=0
			while exit==0 && t2<nTestSample#t2<=(size(testOrder,1)-1)
				t2+=1
				v=testOrder[t2]
				x_forced0_v=ATAmodel.Settings.forced0[v]
				x_new, TIF_new, infeas_new, ol_new = copy(x_old), copy(TIF_old), copy(infeas_old), copy(ol_old)
				Constraints=ATAmodel.Constraints[v]
				IIFv=IIF[v]
				x_v=x_new[:,v]
				takenItems=findall(x_v.==1)
				if nI>size(takenItems,1)
					nI=Int(size(takenItems,1))
				end
				takenItems=Random.shuffle!(takenItems)
				exit=0
				h2=0
				while exit==0 && h2<nI
					h2+=1
					x_new, TIF_new, infeas_new, ol_new = copy(x_old), copy(TIF_old), copy(infeas_old),copy(ol_old)
					#try to remove h
					h=takenItems[h2]
					x_new[h,v]=zero(Float64)
					x_v=x_new[:,v]
					infeas_new[v], x_I_v=checkFeas(infeas_new[v],ATAmodel.Settings.FS,Constraints,x_new,x_v,nFS,nItems,v)
					iu=sum(x_new,dims=2)-ATAmodel.IU.Max
					iu=iu[iu.>0]
					if size(iu,1)==0
						iu=0
					else
						iu=sum(iu)
					end
					ol=evalOverlap(x_new,FScounts,olMax,T)
					for v=1:ATAmodel.Settings.T
						ol_v=ol[setdiff(collect(1:T),v),v]
						ol_v=ol_v[ol_v.>0]
						if size(ol_v,1)==0
							ol_new[v]=0
						else
							ol_new[v]=sum(ol_v)
						end
					end
					if findFeasible==false
						if ATAmodel.Settings.OptType=="MAXIMIN"
							TIF_new[v]=evalTIFMMv(x_I_v,IIFv)
						elseif ATAmodel.Settings.OptType=="CC"
							TIF_new[v]=evalTIFCCv(x_I_v,IIFv;α=α)
						end
					end
					f_new=(1-OptFeas)*(sum(infeas_new+ol_new)+iu)-OptFeas*minimum(TIF_new)
					f_evals+=1

					if (f_new <= f_old)
						#remove item
						x_old, f_old,  TIF_old, infeas_old, ol_old = copy(x_new), copy(f_new), copy(TIF_new) , copy(infeas_new), copy(ol_new)
						if (f_new < f_opt)
							println("better to remove, new f_opt= ",f_new)
							exit=1
							convergence=0
							nI=copy(nItemSample)
							x_opt, f_opt,  TIF_opt, infeas_opt, ol_opt = copy(x_new),copy(f_new),copy(TIF_new), copy(infeas_new), copy(ol_new)
						end
					else
						p = exp_c(-(f_new - f_old) / t)
						if (rand() < p)
							#remove
							exit=1
							x_old, f_old,  TIF_old, infeas_old, ol_old = copy(x_new), copy(f_new), copy(TIF_new) , copy(infeas_new), copy(ol_new)
						end
					end
					#try to switch h with i2, warmpu strategy, take the best to swwitch
					# if OptFeas==0 || findFeasible==true
					# 	x_new, f_new, infeas_new, ol_new, iu= WarmUp(x_new,infeas_new,v,ATAmodel.IU,ATAmodel.Settings.FS,Constraints,x_forced0_v,nItems,nFS,FScounts,olMax,ol_new)
					# else
					# 	if ATAmodel.Settings.OptType=="MAXIMIN"
					# 		x_new, f_new, TIF_new, infeas_new, ol_new, iu = WarmUp(x_new,infeas_new,IIF[v],TIF_new,OptFeas,v,ATAmodel.IU,ATAmodel.Settings.FS,Constraints,x_forced0_v,nItems,nFS,FScounts,olMax,ol_new)
					# 	elseif ATAmodel.Settings.OptType=="CC"
					# 		x_new, f_new, TIF_new, infeas_new, ol_new, iu = WarmUp(x_new,infeas_new,IIF[v],TIF_new,α,OptFeas,v,ATAmodel.IU,ATAmodel.Settings.FS,Constraints,x_forced0_v,nItems,nFS,FScounts,olMax,ol_new)
					# 	end
					# end
					idx_v2=setdiff(collect(1:nFS),findall(x_v.==one(Float64))) #i2=1,...,I
					idx_v2=Random.shuffle!(idx_v2)
					x_start=copy(x_new)
					iu_old=sum(x_start,dims=2)-ATAmodel.IU.Max
					betterFound=0
					i3=1
					# infeas_old_v=copy(infeas_old_v)
					# ol_old_v=copy(ol_old[v])
					# f_old_v=(1-OptFeas)*(sum(infeas_old_v+ol_old_v)-OptFeas*TIF_old[v]

					while betterFound==0 && i3<=size(idx_v2,1)
						i2=idx_v2[i3]
						if x_forced0_v[i2]
							x_new, TIF_new, infeas_new, iu=copy(x_start), copy(TIF_old), copy(infeas_old), copy(iu_old)
							x_new[i2,v]=one(Float64)
							x_v=copy(x_new[:,v])
							iu[i2]= iu[i2]+1
							iu=iu[iu.>0]
							if size(iu,1)==0
								iu=0
							else
								iu=sum(iu)
							end
							infeas_new[v], x_I_v=checkFeas(infeas_new[v],ATAmodel.Settings.FS,Constraints,x_new,x_v,nFS,nItems,v)
							if findFeasible==false
								if ATAmodel.Settings.OptType=="MAXIMIN"
									TIF_new[v]=evalTIFMMv(x_I_v,IIFv)
								elseif ATAmodel.Settings.OptType=="CC"
									TIF_new[v]=evalTIFCCv(x_I_v,IIFv;α=α)
								end
							end
							ol=evalOverlap(x_new,FScounts,olMax,T)
							for v2=1:ATAmodel.Settings.T
								ol_v=ol[setdiff(collect(1:T),v2),v2]
								ol_v=ol_v[ol_v.>0]
								if size(ol_v,1)==0
									ol_new[v2]=0
								else
									ol_new[v2]=sum(ol_v)
								end
							end
							f_new=(1-OptFeas)*(sum(infeas_new+ol_new)+iu)-OptFeas*minimum(TIF_new)
							if (f_new <= f_old)
								#remove item
								x_old, f_old,  TIF_old, infeas_old, ol_old = copy(x_new), copy(f_new), copy(TIF_new) , copy(infeas_new), copy(ol_new)
								if (f_new < f_opt)
									println("better to switch, new f_opt= ",f_new)
									exit=1
									betterFound=1
									convergence=0
									nI=copy(nItemSample)
									x_opt, f_opt,  TIF_opt, infeas_opt, ol_opt = copy(x_new),copy(f_new),copy(TIF_new), copy(infeas_new), copy(ol_new)
								end
							else
								p = exp_c(-(f_new - f_old) / t)
								if (rand() < p)
									#remove
									exit=1
									betterFound=1
									x_old, f_old,  TIF_old, infeas_old, ol_old = copy(x_new), copy(f_new), copy(TIF_new) , copy(infeas_new), copy(ol_new)
								end
							end
						end
						i3+=1
					end
					# if nItems>nFS
					# 	x_I_v=tranformFStoItems(nFS,nItems,x_new[:,v],ATAmodel.Settings.FS.Items)
					# else
					# 	x_I_v=x_new[:,v]
					# end
					#warmpu strategy from here
					# ol=evalOverlap(x_new,FScounts,olMax,T)
					# for v=1:ATAmodel.Settings.T
					# 	ol_v=ol[setdiff(collect(1:T),v),v]
					# 	ol_v=ol_v[ol_v.>0]
					# 	if size(ol_v,1)==0
					# 		ol_new[v]=0
					# 	else
					# 		ol_new[v]=sum(ol_v)
					# 	end
					# end
					# iu=sum(x_new,dims=2)-ATAmodel.IU.Max
					# iu=iu[iu.>0]
					# if size(iu,1)==0
					# 	iu=0
					# else
					# 	iu=sum(iu)
					# end
					# if findFeasible==false
					# 	if ATAmodel.Settings.OptType=="MAXIMIN"
					# 		TIF_new[v]=evalTIFMMv(x_I_v,IIFv)
					# 	elseif ATAmodel.Settings.OptType=="CC"
					# 		TIF_new[v]=evalTIFCCv(x_I_v,IIFv;α=α)
					# 	end
					# end
					# f_new=(1-OptFeas)*(sum(infeas_new+ol_new)+iu)-OptFeas*minimum(TIF_new)
					# if f_new <= f_old
					# 		x_old, f_old,  TIF_old, infeas_old, ol_old = copy(x_new), copy(f_new), copy(TIF_new) , copy(infeas_new), copy(ol_new)
					# 		if (f_new < f_opt)
					# 			println("better to switch, new f_opt= ",f_new)
					# 			exit=1
					# 			convergence=0
					# 			nI=copy(nItemSample)
					# 			x_opt, f_opt,  TIF_opt, infeas_opt, ol_opt = copy(x_old),copy(f_old),copy(TIF_old), copy(infeas_old), copy(ol_old)
					# 		end
					# 	else
					# 		#remove if p>rand()
					# 		p = exp_c(-(f_new - f_old) / t)
					# 		if (rand() < p)
					# 			exit=1
					# 			x_old, f_old,  TIF_old, infeas_old, ol_old = copy(x_new), copy(f_new), copy(TIF_new) , copy(infeas_new), copy(ol_new)
					# 		end
					# end
				    #warmpu strategy to here

					#alternative
					# x_start=copy(x_new)
					# i3=0
					# f_v_old=(1-OptFeas)*(infeas_old[v]+ol_old[v]+iu)-(OptFeas*TIF_old[v])
					# println("f_v_old: ",f_v_old)
					# iu_old=sum(x_old,dims=2)-ATAmodel.IU.Max
					# while exit==0 && i3<size(idx_vav,1)
					# 	i3+=1
					# 	i2=idx_vav[i3]
					# 	if i2!=h && x_forced0_v[i2]
					# 		x_new=copy(x_start)
					# 		x_new[i2,v]=one(Float64)
					# 		x_v=x_new[:,v]
					# 		infeas_new[v], x_I_v=checkFeas(infeas_new[v],ATAmodel.Settings.FS,Constraints,x_new,x_v,nFS,nItems,v)
					# 		ol=evalOverlapv(x_v,x_new,FS.Counts,olMax[:,v],v)
					# 		iu=copy(iu_old)
					# 		iu[i2]=iu[i2]+1-ATAmodel.IU.Max[i2]
					# 		iu=iu[iu.>0]
					# 		if size(iu,1)==0
					# 			iu=0
					# 		else
					# 			iu=sum(iu)
					# 		end
					# 		if findFeasible==false
					# 			if ATAmodel.Settings.OptType=="MAXIMIN"
					# 				TIF_new[v]=evalTIFMMv(x_I_v,IIFv)
					# 			elseif ATAmodel.Settings.OptType=="CC"
					# 				TIF_new[v]=evalTIFCCv(x_I_v,IIFv;α=α)
					# 			end
					# 		end
					# 		f_v_new=(1-OptFeas)*(infeas_new[v]+ol+iu)-(OptFeas*TIF_new[v])
					# 		if f_v_new <= f_v_old
					# 			ol=evalOverlap(x_new,FScounts,olMax,T)
					# 			for v=1:T
					# 				ol_v=ol[setdiff(collect(1:T),v),v]
					# 				ol_new[v]=0.5*sum(ol_v[ol_v.>0])
					# 			end
					# 			f_old=(1-OptFeas)*(sum(infeas_new+ol_new)+iu)-OptFeas*minimum(TIF_new)
					# 			x_old,  TIF_old, infeas_old, ol_old = copy(x_new), copy(TIF_new) , copy(infeas_new), copy(ol_new)
					# 			if (f_v_new < f_v_old)
					# 				println("f_v_new: ",f_v_new)
					# 				println("better to switch, new f_opt= ",f_new)
					# 				#Printf.@printf "="
					# 				exit=1
					# 				convergence=0
					# 				nT=copy(nTestSample)
					# 				nI=copy(nItemSample)
					# 				x_opt, f_opt,  TIF_opt, infeas_opt, ol_opt = copy(x_old),copy(f_old),copy(TIF_old), copy(infeas_old), copy(ol_old)
					# 				#nnew +=one(Float64)
					# 			end
					# 		else
					# 			#remove if p>rand()
					# 			p = exp_c(-(f_v_new - f_v_old) / t)
					# 			if (rand() < p)
					# 				exit=1
					# 				ol=evalOverlap(x_new,FScounts,olMax,T)
					# 				for v=1:T
					# 					ol_v=ol[setdiff(collect(1:T),v),v]
					# 					ol_new[v]=0.5*sum(ol_v[ol_v.>0])
					# 				end
					# 				f_old=(1-OptFeas)*(sum(infeas_new+ol_new)+iu)-OptFeas*minimum(TIF_new)
					# 				x_old,  TIF_old, infeas_old, ol_old = copy(x_new), copy(TIF_new) , copy(infeas_new), copy(ol_new)
					# 			else
					# 			end
					# 		end#end of if f_v_new<=f_v_old
					# 	end #end of i2!=h
					# end #end of i3
				end #end of itemorder h2
			end #end of testorder t2

			if verbosity > 1
				println(hline)
				println("Intermediate results before next temperature change")
				println("temperature: ", t)
				println("current local Infeasibility: ", sum(infeas_opt))
				println("current local Optimality: ", minimum(TIF_opt))
				println("total evaluations so far: ", f_evals)
				println("total time elapsed: ", ATAmodel.Output.ElapsedTime)
				println(hline*"\n") #print results
			end
			fstar[1] = f_old #it was f_old
			if fstar[2]==fstar[1]
				convergence+=1
			end
			#println("convergence is ", convergence)
			#how many equal f_old in the last iterations?
			if fstar[2]==fstar[1] && convergence==convMax#00
				println("neighbourhood ",NH," fully explored, increase temperature")
				nI=copy(nItemSample)
				coverage_ok=1
				t = copy(temp_0)
				ATAmodel.Output.Neighborhoods[NH].f=copy(f_opt)
				ATAmodel.Output.Neighborhoods[NH].Design=copy(x_opt)
				ATAmodel.Output.Neighborhoods[NH].Feas=copy(infeas_opt+ol_opt)
				ATAmodel.Output.Neighborhoods[NH].TIF=copy(TIF_opt)
				fstar = typemax(Float64)*ones(2)
				NH+=1
				if NH>nNH
					converge=1
				else
					if NH>feasNH
						findFeasible=false
					end
					#perturbate the solution and go to warmup
					j=1
					for v=1:ATAmodel.Settings.T
						for i=1:nFS
							if x_opt[i,v]==one(Float64)
								if j==1
									x_opt[i,v]=zero(Float64)
									j+=1
								else
									j-=1
								end
							end
						end
					end

				end
			else
				if coverage_ok==1
					println("neighbourhood ",NH," not entirely explored, decrease temperature")
				end
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
		length_max=[ATAmodel.Constraints[t].length_max  for t=1:T]
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
		if GroupByFS==true
			design=reshape(ATAmodel.Output.Design,nFS,T)
			new_design=zeros(Float64,nItems,T)
			for t=1:T
				new_design[vcat(ATAmodel.Settings.FS.Items[findall(design[:,t].==one(Float64))]...),t].=one(Float64)
			end
			design=new_design
			writedlm("RESULTS/designItemLevel.csv",design)
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
				DataFrames.names!(TIFatTheta_k,Symbol.(["Min", "0.25-Qle", "Median", "0.75-Qle", "Max" ,L"{\alpha}-Qle","Estimated","True"]))
			else
				DataFrames.names!(TIFatTheta_k,Symbol.(["Min", "0.25-Qle", "Median", "0.75-Qle", "Max" ,L"{\alpha}-Qle","Estimated"]))
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
			ESprintIRT[t]=((ICF[t]'*design[:,t]))[1]/sum(design[:,t])
		end
		#number of items
		n=Array{Float64,1}(undef,T)
		for t in 1:T
			n[t]=sum(design[:,t])
		end
		designItems=zeros(Int64,Int(maximum(n)),T)
		#designItems=Array{String,2}(nothing,n_max,T)
		for t in 1:T
			designItems[collect(1:Int(n[t])),t].=sort!(findall(design[:,t].==1))#sort!(CODEVar[findall(design[:,t].==1)])
		end
		open("RESULTS/Results.txt", "w") do io
			write(io,"Length")
			write(io,"\r\n")
			writedlm(io, [1:T,Int.(n)'],"\t")
			write(io, "\r\n")
			write(io, "Design")
			write(io,"\r\n")
			writedlm(io, designItems,"\t")
			write(io,"\r\n")
			write(io, "expected score")
			write(io,"\r\n")
			writedlm(io, [1:T,round.(ESprintIRT,digits=3)'],"\t")
			if size(CATEGORIES,1)>0
				write(io,"Categorical variables")
				write(io,"\r\n")
				for cats in CATEGORIES
					write(io,cats)
					write(io,"\r\n")
					vals1=setdiff(unique(ATAmodel.Settings.Bank[Symbol(cats)]),[missing])
					vals=copy(vals1)
					for t=1:T
						valscat=Int64[]
						for val in vals1
							valscat_t=ATAmodel.Settings.Bank[Symbol(cats)][design[:,t].==1].== val
							valscat_t=valscat_t[.!ismissing.(valscat_t)]
							if size(valscat_t,1)>0
								valscat=vcat(valscat,sum(valscat_t))
							else
								valscat=vcat(valscat,0)
							end
						end
						vals=hcat(vals,valscat)
					end
					writedlm(io,vals,"\t")
					#writedlm(io,hcat(vals,[sum(skipmissing(ATAmodel.Settings.Bank[Symbol(cats)][design[:,t].==1].== i))   for i in vals,t=1:T]),", ")
					write(io,"\r\n")
				end
				write(io,"Sum variables")
				write(io,"\r\n")
				for var in unique(vcat([ATAmodel.Constraints[t].sumVars for t=1:T]...))
					write(io,var)
					write(io,"\r\n")
					writedlm(io,round.([sum(skipmissing(ATAmodel.Settings.Bank[var].*design[:,t])) for t=1:T],digits=3)',"\t")
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
	end
end
