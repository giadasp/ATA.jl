function eval_overlapv(x_v::Vector{Float64}, x::Matrix{Float64}, FScounts::Vector{Int64}, ol_max_v::Vector{Float64}, v::Int64)#idx_v::Vector{Int64}, idxnotv::Vector{Vector{Int64}}) #x = I*T , f1 = idx[f1] and f2 = idx[f2]
	x_v = x_v.*FScounts
	cons = LinearAlgebra.BLAS.gemv('T', x, x_v)-ol_max_v
	cons[v] = 0
	# for testnotv = 1:Tminus1
	# 	cons[testnotv] = LinearAlgebra.dot(x_v, x[idxnotv[testnotv]])-ol_max_v[testnotv]#x[idx_v]'*(FScounts.*x[idxnotv[testnotv]])-ol_max_v[testnotv] #	cons[testnotv] = ((x[idx_v])'*(FScounts.*x[idxnotv[testnotv]]))-ol_max_v[testnotv]
	# end
	return sum(cons[cons.>0])
end

function eval_overlap(x::Matrix{Float64}, FScounts::Matrix{Float64}, ol_max::Matrix{Float64}, T::Int64, ol2::Vector{Float64})
	FScounts = FScounts.*x
	ol = gemmblas(x, FScounts, ol_max)
	for v2 = 1:T
		ol_v = view(ol, : , v2)
		ol_v[v2] = zero(Float64)
		ol_v = ol_v[ol_v .> 0]
		if size(ol_v, 1) == 0
			ol2[v2] = 0
		else
			ol2[v2] = sum(ol_v)
		end
	end
	return ol2::Vector{Float64}
end

function eval_TIF_CC_v(x_v::Vector{Float64}, IIF::Array{Float64, 3}; α = 0.1) # x = I
	K, I, R = size(IIF)
	alphaR = Int(ceil(R*(α)))
	αQle = Inf
	if K>1
		for k = 1:K
			αQle = min(αQle, sort(LinearAlgebra.BLAS.gemv('T', IIF[k, :, :], x_v))[alphaR])
		end
	else
		αQle = sort(gemvblasT(IIF[1, :, :], x_v, R))[alphaR]
	end
	return αQle
end

function eval_Exp_Score_v(x_v::Vector{Float64}, FSItems::Vector{Vector{Int64}}, ES::Array{Float64, 1}) #x = I
	cons = LinearAlgebra.dot(ES, x_v)
	return cons
end

function eval_TIF_MM_v(x_v::Vector{Float64}, IIF::Matrix{Float64}) # x = I
	K, I = size(IIF)
	if K>1
		TIF = Inf
		for k = 1:K
			TIF = min(TIF, LinearAlgebra.dot(IIF[1, :], x_v))#IIF[k, :]'*x) #minimum for all thetasopt
		end
	else
		TIF = LinearAlgebra.BLAS.gemv('N', IIF, x_v)[1]#IIF[1, :]'*x
	end
	return TIF::Float64
end

function eval_IU(x::Matrix{Float64}, nFS::Int64, T::Int64) #x = I*T
	ItemUse = x * ones(T)
	return ItemUse
end


function check_feas(FS::FS, Consᵥ::Constraint, xᵥ::Vector{Float64}, nFS::Int64, n_items::Int64, v::Int64)
	es = 0
	cons = 0
	if size(Consᵥ.constr_A, 1)>0
		cons = copy(Consᵥ.constr_b)
		cons = gemvblas(Consᵥ.constr_A, xᵥ, cons, size(xᵥ, 1))
	end
	if n_items>nFS
		x_Iᵥ = FS_to_items(nFS, n_items, xᵥ, FS.items)
	else
		x_Iᵥ = copy(xᵥ)
	end
	if size(Consᵥ.expected_score.val, 1) > 0
		es = gemvblas(Consᵥ.expected_score.val, x_Iᵥ, zeros(size(Consᵥ.expected_score.val,1)),size(x_Iᵥ, 1)) / sum(x_Iᵥ)
		es = vcat((es - Consᵥ.expected_score.max)..., (- es + Consᵥ.expected_score.min)...)
	end
	cons = vcat(cons, es)#, 0.5 .*ol)
	cons = cons[cons .> 0]
	if size(cons, 1) > 0
		cons = sum(cons)
	else
		cons = zero(Float64)
	end
	return cons::Float64, x_Iᵥ::Vector{Float64}
end

function FS_to_items(nFS::Int64, n_items::Int64, xᵥ::Vector{Float64}, FSitems::Vector{Vector{Int64}})
	xᵥ_taken = findall(xᵥ .== one(Float64))
	x_Iᵥ = zeros(Float64, n_items)
	for i in xᵥ_taken
		x_Iᵥ[FSitems[i]].= one(Float64)
	end
	return x_Iᵥ::Vector{Float64}
end

#warmup feasibility
function fill_up(NH::Neighbourhood, v::Int64, IU::IU, FS::FS, constraints::Constraint, x_forced0ᵥ::Vector{Bool}, n_items::Int64, nFS::Int64, FScounts::Matrix{Float64}, ol_maxᵥ::Vector{Float64})
	f₀, infeas₀ = Inf, Inf
	T = size(NH.x, 2)
	x₀ = copy(NH.x)
	xᵥ = x₀[:, v]
	idxᵥ₂ = setdiff(collect(1:nFS), findall(xᵥ .== one(Float64))) #i₂ = 1, ..., I
	#idxᵥ₂ = Random.shuffle!(idxᵥ₂) #removed
	i⁺ = copy(idxᵥ₂[1])
	iu⁺ = 0
	iu₀ = sum(x₀, dims = 2) - IU.max
	ol₀ᵥ = 0
	for i₂ in idxᵥ₂
		if x_forced0ᵥ[i₂]
			x₁, infeas₁, iu = copy(x₀), copy(infeas₀), copy(iu₀)
			x₁[i₂, v] = one(Float64)
			iu[i₂] = iu[i₂] + 1
			iu = iu[iu .> 0]
			if size(iu, 1) == 0
				iu = 0
			else
				iu = sum(iu)
			end
			xᵥ = copy(x₁[:, v])
			infeas₁, x_Iᵥ = check_feas(FS, constraints, xᵥ, nFS, n_items, v)
			ol = eval_overlapv(xᵥ, x₁, FS.counts, ol_maxᵥ, v)
			f₁ = infeas₁ + iu + ol
			if (f₁<f₀)
				f₀, infeas₀ = copy(f₁), copy(infeas₁)
				iu⁺ = copy(iu)
				ol₀ᵥ = copy(ol)
				i⁺ = copy(i₂)
			end
		end
	end
	NH.ol[v] = copy(ol₀ᵥ)
	NH.infeas[v] = copy(infeas₀)
	NH.iu = copy(iu⁺)
	NH.x[i⁺, v] = one(Float64)
	return NH::Neighbourhood
end

#warmup MAXIMIN
function fill_up(NH::Neighbourhood, IIFv::Matrix{Float64}, opt_feas::Float64, v::Int64, IU::IU, FS::FS, constraints::Constraint, x_forced0ᵥ::Vector{Bool}, n_items::Int64, nFS::Int64, FScounts::Matrix{Float64}, ol_maxᵥ::Vector{Float64})
	f₀, TIF₀, infeas₀ = Inf, copy(NH.obj), Inf
	T = size(NH.x, 2)
	x₀ = copy(NH.x)
	idxᵥ₂ = setdiff(collect(1:nFS), findall(x₀[:, v] .== one(Float64))) #i₂ = 1, ..., I
	#idxᵥ₂ = Random.shuffle!(idxᵥ₂) #removed
	i⁺ = 1
	iu₀ = sum(NH.x, dims = 2) - IU.max
	ol₀ᵥ = 0
	iu⁺ = 0
	for i₂ in idxᵥ₂
		if x_forced0ᵥ[i₂]
			x₁, TIF₁, infeas₁, iu = copy(x₀), copy(TIF₀), copy(infeas₀), copy(iu₀)
			x₁[i₂, v] = one(Float64)
			xᵥ = copy(x₁[:, v])
			iu[i₂] = iu[i₂] + 1
			iu = iu[iu.>0]
			if size(iu, 1) == 0
				iu = 0
			else
				iu = sum(iu)
			end
			infeas₁, x_Iᵥ = check_feas(FS, constraints, xᵥ, nFS, n_items, v)
			TIF₁[v] = eval_TIF_MM_v(x_Iᵥ, IIFv)
			ol = eval_overlapv(xᵥ, x₁, FS.counts, ol_maxᵥ, v)
			f₁ = (1-opt_feas) * (infeas₁+iu+ol) - opt_feas * minimum(TIF₁)
			if (f₁<f₀)
				i⁺ = copy(i₂)
				f₀, TIF₀, infeas₀ = copy(f₁), copy(TIF₁), copy(infeas₁)
				iu⁺ = copy(iu)
				ol₀ᵥ = copy(ol)
			end
		end
	end
	NH.ol[v] = copy(ol₀ᵥ)#0
	NH.infeas[v] = copy(infeas₀)
	NH.iu = copy(iu⁺)
	NH.obj = copy(TIF₀)
	NH.x[i⁺, v] = one(Float64)
	return NH::Neighbourhood
end

#warmup MAXIMIN CC
function fill_up(NH::Neighbourhood, IIFv::Array{Float64, 3}, α::Float64, opt_feas::Float64, v::Int64, IU::IU, FS::FS, constraints::Constraint, x_forced0ᵥ::Vector{Bool}, n_items::Int64, nFS::Int64, FScounts::Matrix{Float64}, ol_maxᵥ::Vector{Float64})
	f₀, TIF₀, infeas₀ = Inf, copy(NH.obj), Inf
	T = size(NH.x, 2)
	x₀ = copy(NH.x)
	idxᵥ₂ = setdiff(collect(1:nFS), findall(x₀[:, v] .== one(Float64))) #i₂ = 1, ..., I
	#idxᵥ₂ = Random.shuffle!(idxᵥ₂) #removed
	i⁺ = 1
	ol₀ᵥ = 0
	iu₀ = sum(NH.x, dims = 2) - IU.max
	iu⁺ = 0
	for i₂ in idxᵥ₂
		if x_forced0ᵥ[i₂]
			x₁, TIF₁, infeas₁, iu = copy(x₀), copy(TIF₀), copy(infeas₀), copy(iu₀)
			x₁[i₂, v] = one(Float64)
			xᵥ = copy(x₁[:, v])
			iu[i₂] = iu[i₂] + 1
			iu = iu[iu .> 0]
			if size(iu, 1) == 0
				iu = 0
			else
				iu = sum(iu)
			end
			infeas₁, x_Iᵥ = check_feas(FS, constraints, xᵥ, nFS, n_items, v)
			TIF₁[v] = eval_TIF_CC_v(x_Iᵥ, IIFv; α = α)
			ol = eval_overlapv(xᵥ, x₁, FS.counts, ol_maxᵥ, v)
			f₁ = (1-opt_feas) * (infeas₁ + iu + ol) - opt_feas*minimum(TIF₁)
			if (f₁ < f₀)
				i⁺ = copy(i₂)
				f₀, TIF₀, infeas₀ = copy(f₁), copy(TIF₁), copy(infeas₁)
				iu⁺ = copy(iu)
				ol₀ᵥ = copy(ol)
			end
		end
	end
	NH.ol[v] = copy(ol₀ᵥ)
	NH.infeas[v] = copy(infeas₀)
	NH.iu = copy(iu⁺)
	NH.obj = copy(TIF₀)
	NH.x[i⁺, v] = one(Float64)
	return NH::Neighbourhood
end

#MAXIMIN CC neighbourhood
function analyse_NH(NH_start::Neighbourhood, ATAmodel::Model, IIF::Vector{Array{Float64, 3}}; fF = true, n_fill = 1, opt_feas = 0.9, conv_max = 1, start_temp = 1000.0, geom_temp = 0.1, n_item_sample = 1, n_test_sample = 1, verbosity = 1, start_time = 0, max_time = 1000)
	if fF == true
		NH_start.obj = zeros(Float64, ATAmodel.settings.T)
	end
	NH₁ = Neighbourhood()
	NH₁ = mycopy(NH_start, NH₁)
	NH₀ = Neighbourhood()
	NH⁺ = Neighbourhood()
	f_star = ones(2).*Inf
	f_evals = 0
	t = copy(start_temp)
	T = ATAmodel.settings.T
	n_items = ATAmodel.settings.n_items
	nFS = ATAmodel.settings.nFS
	FScounts = ATAmodel.settings.FS.counts * ones(Float64, T)'
	switches = 0
	removes = 0
	adds = 0
	#warm up
	println("fill-up starting")
	round = 1

	if n_fill>0
		for round = 1:n_fill
			NH₁ = mycopy(NH_start, NH₁)
			warmup = fill(true, T)
			while any(warmup) #filling forms
				v = findfirst(warmup .== true)
				constraints = ATAmodel.constraints[v]
				#println("ol = ", NH₁.ol)
				fᵥ = (1-opt_feas) * (NH₁.infeas+NH₁.ol) - (opt_feas * NH₁.obj)
				mm = fᵥ[v]
				for v2 in findall(warmup .== true)
					if fᵥ[v2] > mm
						mm = copy(fᵥ[v2])
						v = copy(v2)
					end
				end
				#filling
				n_t = LinearAlgebra.dot(NH₁.x[:, v], ATAmodel.settings.FS.counts)
				#try to add other items, the first time it goes over n_max it stops
				if n_t<constraints.length_max
					if opt_feas == 0 || fF == true
						NH_add = fill_up(NH₁, v, ATAmodel.IU, ATAmodel.settings.FS, constraints, ATAmodel.settings.forced0[v], n_items, nFS, FScounts, ATAmodel.settings.ol_max[:, v])
						NH_add.f = opt_feas * NH_add.f
					else
						if ATAmodel.settings.opt_type == "MAXIMIN"
							NH_add = fill_up(NH₁, IIF[v], opt_feas, v, ATAmodel.IU, ATAmodel.settings.FS, constraints, ATAmodel.settings.forced0[v], n_items, nFS, FScounts, ATAmodel.settings.ol_max[:, v])
						elseif ATAmodel.settings.opt_type == "CC"
							NH_add = fill_up(NH₁, IIF[v], ATAmodel.obj.aux_float, opt_feas, v, ATAmodel.IU, ATAmodel.settings.FS, constraints, ATAmodel.settings.forced0[v], n_items, nFS, FScounts, ATAmodel.settings.ol_max[:, v])
						end
					end
					NH_add.ol = eval_overlap(NH_add.x, FScounts, ATAmodel.settings.ol_max, T, NH_add.ol)
					Printf.@printf "."
					nᵥ = LinearAlgebra.dot(NH_add.x[:, v], ATAmodel.settings.FS.counts)
					#println("length for test ", v, ": ", nᵥ)
					if n_t<= constraints.length_max
						NH₁ = mycopy(NH_add, NH₁)
					else
						warmup[v] = false
						println("f₁:", NH₁.f)
						println("-Test ", v, " filled up with ", n_t, " items, ")
					end
				else
					warmup[v] = false
					println("-Test ", v, " filled up with ", n_t, " items, ")
				end
			end
		end #end of round, round = nRound
		NH₀ = mycopy(NH₁, NH₀)
		NH₀.f = comp_f(NH₀, opt_feas)
	end
	NH⁺ = mycopy(NH₀, NH⁺)
	println("end of fill-up")
	if sum(NH₀.infeas + NH₀.ol) + NH₀.iu == 0
		println("Feasible solution found in fill-up")
		fF = false
	else
		println("Feasible solution not found in fill-up")
	end
	println("f = ", NH⁺.f)
	println("infeas = ", NH⁺.infeas)
	println("ol = ", NH⁺.ol)
	println("iu = ", NH⁺.iu)
	# statistics to repogeom_temp at each temp change, set back to zero
	coverage_ok = 0
	if n_test_sample>T
		n_test_sample = T
	end
	convergence = 0

	while coverage_ok == 0
		NH₁ = mycopy(NH₀, NH₁)
		weights = (1 - opt_feas) .* (NH₀.infeas + NH₀.ol) - opt_feas .* (NH₀.obj)
		weights = (weights .- minimum(weights)) .+ 1
		weights = weights ./ sum(weights)
		weights = StatsBase.ProbabilityWeights(weights)
		#determine test order
		#testOrder = StatsBase.sample(collect(1:T), weights, n_test_sample, replace = false)
		#testOrder = Random.shuffle!(collect(1:ATAmodel.settings.T))
		testOrder = sortperm(weights, rev = true)
		#println("test order = ", Int.(testOrder))
		v₂ = 0
		exit = 0
		#iteratorTestItem = vec(collect(Iterators.product(collect(1:n_item_sample), testOrder[1:n_test_sample])))
		#println(iteratorTestItem[1:nItemoStatsBase.sample])
		#it = 0
		#xnew = copy(NH₀.x)
		while exit == 0 && v₂ < n_test_sample #it<size(iteratorTestItem, 1) #
			#it+= 1
			v₂ += 1
			exit = 0
			v = testOrder[v₂]
			x_forced0ᵥ = ATAmodel.settings.forced0[v]
			Constraintsᵥ = ATAmodel.constraints[v]
			IIFᵥ = IIF[v]
			ol_maxᵥ = ATAmodel.settings.ol_max[:, v]
			#it<size(iteratorTestItem, 1) #
			#v = iteratorTestItem[it][2]
			#NH₀.x = copy(xnew)
			taken_items = findall(NH₀.x[:, v] .== 1)
			if n_item_sample > size(taken_items, 1)
				nI = Int(size(taken_items, 1))
			else
				nI = n_item_sample
			end
			#taken_items = Random.shuffle!(taken_items) #reset #removed
			# exit2 = 0
			# if iteratorTestItem[it][1]>size(taken_items, 1)
			# 	exit2 = 1
			# end
			#if exit2 == 0
			# if v!= v2
			# 	v = copy(v2)
			#
			# end
			#h = taken_items[iteratorTestItem[it][1]]
			#println("test ", v₂, " of ", size(testOrder, 1))
			#fix test features
			h₂ = 0
			while exit == 0 && h₂ < nI
				add_remove = 0
				while  exit == 0 && add_remove < 2
					add_remove += 1
					NH₁ = mycopy(NH₀, NH₁)
					h₂ += 1
					#println("item ", h₂, " of ", size(taken_items, 1))
					#try to remove h
					h = taken_items[h₂]
					#iu = copy(iu₀)
					#if rand()>0.0
					if add_remove==1 #remove
						NH₁.x[h, v] = zero(Float64)
					end#else add (not remove)
					#iu[h]-= 1
					taken = 0
					#else
					#	taken = 1
					#end
					if (add_remove==1 && sum(NH₁.x[:, v] .* ATAmodel.settings.FS.counts) >= Constraintsᵥ.length_min) || (add_remove==2 && sum(NH₁.x[:, v] .* ATAmodel.settings.FS.counts) < Constraintsᵥ.length_max)
						NH₁.infeas[v], x_Iᵥ = check_feas(ATAmodel.settings.FS, Constraintsᵥ, NH₁.x[:, v], nFS, n_items, v)
						iu = sum(NH₁.x, dims = 2) - ATAmodel.IU.max
						iu = iu[iu .> 0]
						if size(iu, 1) == 0
							NH₁.iu = 0
						else
							NH₁.iu = sum(iu)
						end
						NH₁.ol = eval_overlap(NH₁.x, FScounts, ATAmodel.settings.ol_max, T, NH₁.ol)
						#NH₁.ol[v] = eval_overlapv(NH₁.x[:, v], NH₁.x, ATAmodel.settings.FS.counts, ol_maxᵥ, v)
						if fF == false
							if ATAmodel.settings.opt_type == "MAXIMIN"
								NH₁.obj[v] = eval_TIF_MM_v(x_Iᵥ, IIFᵥ)
							elseif ATAmodel.settings.opt_type == "CC"
								NH₁.obj[v] = eval_TIF_CC_v(x_Iᵥ, IIFᵥ; α = ATAmodel.obj.aux_float)
							end
						end
						NH₁.f = comp_f(NH₁, opt_feas)
						f_evals+= 1
						if (NH₁.f <= NH₀.f)
							#switch item
							#NH₀.ol = eval_overlap(NH₀.x, FScounts, ATAmodel.settings.ol_max, T, NH₀.ol)
							#NH₀.f = comp_f(NH₀, opt_feas)
							NH₀ = mycopy(NH₁, NH₀)
							# println("better to remove, new f₀ = ", NH₀.f)
							if (NH₁.f < NH⁺.f)
								exit = 1
								convergence = 0
								NH⁺ = mycopy(NH₁, NH⁺)
								if verbosity == 2
									println("f⁺ = ", NH⁺.f, ", t = ", v)
								else
									Printf.@printf "+"
								end
							end
						else
							p = sig_c(-(NH₁.f - NH₀.f) / t) # -(NH₁.f - NH₀.f) / t always negative
							if p<1e-10 p = zeros(Float64)[1]; end

							if (rand() < p)
								#remove
								#exit = 1
								NH₀ = mycopy(NH₁, NH₀)
								#NH₀.ol = eval_overlap(NH₀.x, FScounts, ATAmodel.settings.ol_max, T, NH₀.ol)
								#NH₀.f = comp_f(NH₀, opt_feas)
								if verbosity == 2
									println("SA: f₀ = ", NH₀.f, ", t = ", v)
								else
									Printf.@printf "_"
								end
							end
						end
					end
					idxᵥ₂ = findall(NH₁.x[:, v] .== zero(Float64))#Random.shuffle!(findall(NH₁.x[:, v] .== zero(Float64)))#i₂ = 1, ..., I
					betterFound = 0
					i₃ = 0
					x₋₁ = copy(NH₁.x)
					while betterFound == 0 && i₃<size(idxᵥ₂, 1)
						i₃+= 1
						i₂ = idxᵥ₂[i₃]
						if x_forced0ᵥ[i₂] && i₂!= h
							NH₁ = mycopy(NH₀, NH₁)
							NH₁.x = copy(x₋₁)
							#NH₁.x[h, v] = zero(Float64)
							NH₁.x[i₂, v] = one(Float64)
							# iu = copy(iu₀)
							# iu[i₂]+= 1
							# iu[h]-= 1
							iu = sum(NH₁.x, dims = 2) - ATAmodel.IU.max
							iu = iu[iu .> 0]
							if size(iu, 1) == 0
								NH₁.iu = 0
							else
								NH₁.iu = sum(iu)
							end
							NH₁.infeas[v], x_Iᵥ = check_feas(ATAmodel.settings.FS, Constraintsᵥ, NH₁.x[:, v], nFS, n_items, v)
							if fF == false
								if ATAmodel.settings.opt_type == "MAXIMIN"
									NH₁.obj[v] = eval_TIF_MM_v(x_Iᵥ, IIFᵥ)
								elseif ATAmodel.settings.opt_type == "CC"
									NH₁.obj[v] = eval_TIF_CC_v(x_Iᵥ, IIFᵥ; α = ATAmodel.obj.aux_float)
								end
							end
							NH₁.ol = eval_overlap(NH₁.x, FScounts, ATAmodel.settings.ol_max, T, NH₁.ol)
							#NH₁.ol[v] = eval_overlapv(NH₁.x[:, v], NH₁.x, ATAmodel.settings.FS.counts, ol_maxᵥ, v)
							NH₁.f = comp_f(NH₁, opt_feas)
							if (NH₁.f <= NH₀.f)
								#switch item
								#NH₀ = mycopy(NH₁, NH₀)
								# NH₀.ol = eval_overlap(NH₀.x, FScounts, ATAmodel.settings.ol_max, T, NH₀.ol)
								# NH₀.f = comp_f(NH₀, opt_feas)
								NH₀ = mycopy(NH₁, NH₀)
								# println("better to switch, new f₀ = ", NH₀.f)
								if (NH₁.f < NH⁺.f)
									exit = 1
									betterFound = 1
									convergence = 0
									NH⁺ = mycopy(NH₁, NH⁺)
									if verbosity == 2
										println("f⁺ = ", NH⁺.f, ", t = ", v)
									else
										Printf.@printf "+"
									end
								end
							else
								p = sig_c( - (NH₁.f - NH₀.f) / t)
								if p < 1e-10  p = zeros(Float64)[1]; end
								if (rand() < p)
									#remove
									NH₀ = mycopy(NH₁, NH₀)
									# NH₀.ol = eval_overlap(NH₀.x, FScounts, ATAmodel.settings.ol_max, T, NH₀.ol)
									# NH₀.f = comp_f(NH₀, opt_feas)
									if verbosity == 2
										println("SA: f₀ = ", NH₀.f, ", t = ", v)
									else
										Printf.@printf "_"
									end
								end
							end
						end
					end #end of betterFound (betterFound = 1)
				end #end of add_remove
			end #end of itemorder h₂ (exit = 1)

		end #end of testorder v₂ (exit = 1)
		if sum(NH₀.infeas + NH₀.ol) + NH₀.iu <= 0 && fF == true
			fF = false
			for v = 1:T
				x_Iᵥ = FS_to_items(ATAmodel.settings.nFS, ATAmodel.settings.n_items, NH₀.x[:, v], ATAmodel.settings.FS.items)
				if ATAmodel.settings.opt_type == "MAXIMIN"
					NH₀.obj[v] = eval_TIF_MM_v(x_Iᵥ, IIF[v])
				elseif ATAmodel.settings.opt_type == "CC"
					NH₀.obj[v] = eval_TIF_CC_v(x_Iᵥ, IIF[v]; α = ATAmodel.obj.aux_float)
				end
			end
			NH₀.f = comp_f(NH₀, opt_feas)
			NH⁺ = mycopy(NH₀, NH⁺)
		end
		f_star[1] = copy(NH₀.f)
		#if exit == 0
		if f_star[2] == f_star[1]
			convergence += 1
			Printf.@printf("	%16.10f", convergence)
		end
		#println("convergence is ", convergence)
		#how many equal f₀ in the last iterations?
		if (f_star[2] == f_star[1] && convergence == conv_max) || time() - start_time>= max_time
			coverage_ok = 1
			#t += 1
		else
			t *= geom_temp
			coverage_ok = 0
			if f_star[1] <= f_star[2]
				f_star[2] = copy(f_star[1])
			end
		end
		if verbosity == 2
			println("\n")
			Printf.@printf("\n  f⁺:	%16.10f", NH⁺.f)
			Printf.@printf("\n  f₀:	%16.10f", NH₀.f)
			Printf.@printf("\n  Local convergence:	%16.10f", convergence)
			println("\n")
		end
	end #end of NH coverage (coverage_ok = 1)
	return NH⁺::Neighbourhood
end

#MAXIMIN neighbourhood
function analyse_NH(NH_start::Neighbourhood, ATAmodel::Model, IIF::Vector{Array{Float64, 2}}; fF = true, n_fill = 1, opt_feas = 0.9, conv_max = 1, start_temp = 1000.0, geom_temp = 0.1, n_item_sample = 1, n_test_sample = 1, verbosity = 1, start_time = 0, max_time = 1000)
	if fF == true
		NH_start.obj = zeros(Float64, ATAmodel.settings.T)
	end
	NH₁ = Neighbourhood()
	NH₁ = mycopy(NH_start, NH₁)
	NH₀ = Neighbourhood()
	NH⁺ = Neighbourhood()
	f_star = ones(2).*Inf
	f_evals = 0
	t = copy(start_temp)
	T = ATAmodel.settings.T
	n_items = ATAmodel.settings.n_items
	nFS = ATAmodel.settings.nFS
	FScounts = ATAmodel.settings.FS.counts*ones(Float64, T)'
	#warm up
	println("Fill-up starting")
	round = 1

	if n_fill>0
		for round = 1:n_fill
			NH₁ = mycopy(NH_start, NH₁)
			warmup = fill(true, T)
			while any(warmup) #filling forms
				v = findfirst(warmup .== true)
				constraints = ATAmodel.constraints[v]
				#println("ol = ", NH₁.ol)
				fᵥ = (1 - opt_feas) * (NH₁.infeas + NH₁.ol) - (opt_feas * NH₁.obj)
				mm = fᵥ[v]
				for v2 in findall(warmup .== true)
					if fᵥ[v2] > mm
						mm = copy(fᵥ[v2])
						v = copy(v2)
					end
				end
				#filling
				n_t = LinearAlgebra.dot(NH₁.x[:, v], ATAmodel.settings.FS.counts)
				#println("test ", v, " chosen, length was: ", n_t)
				#try to add other items, the first time it goes over n_max it stops
				if n_t<constraints.length_max
					if opt_feas == 0 || fF == true
						NH_add = fill_up(NH₁, v, ATAmodel.IU, ATAmodel.settings.FS, constraints, ATAmodel.settings.forced0[v], n_items, nFS, FScounts, ATAmodel.settings.ol_max[:, v])
						NH_add.f = opt_feas * NH_add.f
					else
						if ATAmodel.settings.opt_type == "MAXIMIN"
							NH_add = fill_up(NH₁, IIF[v], opt_feas, v, ATAmodel.IU, ATAmodel.settings.FS, constraints, ATAmodel.settings.forced0[v], n_items, nFS, FScounts, ATAmodel.settings.ol_max[:, v])
						end
					end
					NH_add.ol = eval_overlap(NH_add.x, FScounts, ATAmodel.settings.ol_max, T, NH_add.ol)
					Printf.@printf "."
					nᵥ = LinearAlgebra.dot(NH_add.x[:, v], ATAmodel.settings.FS.counts)
					#println("length for test ", v, ": ", nᵥ)
					if n_t<= constraints.length_max
						NH₁ = mycopy(NH_add, NH₁)
					else
						warmup[v] = false
						#NH₁ = mycopy(NH_add, NH₁)
						println("f₁:", NH₁.f)
						println("-Test ", v, " filled up with ", n_t, " items, ")
					end
				else
					warmup[v] = false
					println("-Test ", v, " filled up with ", n_t, " items, ")
				end
			end
		end #end of round, round = nRound
		NH₀ = mycopy(NH₁, NH₀)
		NH₀.f = comp_f(NH₀, opt_feas)
	end
	NH⁺ = mycopy(NH₀, NH⁺)
	println("end of fill-up")
	if sum(NH₀.infeas+NH₀.ol)+NH₀.iu == 0
		println("Feasible solution found in fill-up")
		fF = false
	else
		println("Feasible solution not found in fill-up")
	end
	println("f = ", NH⁺.f)
	println("infeas = ", NH⁺.infeas)
	println("ol = ", NH⁺.ol)
	println("iu = ", NH⁺.iu)
	# statistics to repogeom_temp at each temp change, set back to zero
	coverage_ok = 0
	if n_test_sample > T
		n_test_sample = T
	end
	convergence = 0

	while coverage_ok == 0
		NH₁ = mycopy(NH₀, NH₁)
		weights = (1-opt_feas) .* (NH₀.infeas + NH₀.ol) - opt_feas .* (NH₀.obj)
		weights = (weights .- minimum(weights)).+1
		weights = weights ./ sum(weights)
		weights = StatsBase.ProbabilityWeights(weights)
		#determine test order
		#testOrder = StatsBase.sample(collect(1:T), weights, n_test_sample, replace = false)
		#testOrder = Random.shuffle!(collect(1:ATAmodel.settings.T))
		testOrder = sortperm(weights, rev = true)
		#println("test order = ", Int.(testOrder))
		v₂ = 0
		exit = 0
		#iteratorTestItem = vec(collect(Iterators.product(collect(1:n_item_sample), testOrder[1:n_test_sample])))
		#println(iteratorTestItem[1:nItemoStatsBase.sample])
		#it = 0
		#xnew = copy(NH₀.x)
		while exit == 0 && v₂ < n_test_sample #it<size(iteratorTestItem, 1) #
			#it+= 1
			v₂+= 1
			exit = 0
			v = testOrder[v₂]
			x_forced0ᵥ = ATAmodel.settings.forced0[v]
			Constraintsᵥ = ATAmodel.constraints[v]
			IIFᵥ = IIF[v]
			ol_maxᵥ = ATAmodel.settings.ol_max[:, v]
			#it<size(iteratorTestItem, 1) #
			#v = iteratorTestItem[it][2]
			#NH₀.x = copy(xnew)
			taken_items = findall(NH₀.x[:, v] .== 1)
			if n_item_sample > size(taken_items, 1)
				nI = Int(size(taken_items, 1))
			else
				nI = n_item_sample
			end
			taken_items = taken_items #Random.shuffle!(taken_items) #reset
			# exit2 = 0
			# if iteratorTestItem[it][1]>size(taken_items, 1)
			# 	exit2 = 1
			# end
			#if exit2 == 0
			# if v!= v2
			# 	v = copy(v2)
			#
			# end
			#h = taken_items[iteratorTestItem[it][1]]
			#println("test ", v₂, " of ", size(testOrder, 1))
			#fix test features
			h₂ = 0
			while exit == 0 && h₂<nI
				NH₁ = mycopy(NH₀, NH₁)
				h₂+= 1
				#println("item ", h₂, " of ", size(taken_items, 1))
				#try to remove h
				h = taken_items[h₂]
				#iu = copy(iu₀)
				#if rand()>0.0
				NH₁.x[h, v] = zero(Float64)
				#iu[h]-= 1
				taken = 0
				#else
				#	taken = 1
				#end
				if sum(NH₁.x[:, v].*ATAmodel.settings.FS.counts)>= Constraintsᵥ.length_min
					NH₁.infeas[v], x_Iᵥ = check_feas(ATAmodel.settings.FS, Constraintsᵥ, NH₁.x[:, v], nFS, n_items, v)
					iu = sum(NH₁.x, dims = 2) - ATAmodel.IU.max
					iu = iu[iu .> 0]
					if size(iu, 1) == 0
						NH₁.iu = 0
					else
						NH₁.iu = sum(iu)
					end
					NH₁.ol = eval_overlap(NH₁.x, FScounts, ATAmodel.settings.ol_max, T, NH₁.ol)
					#NH₁.ol[v] = eval_overlapv(NH₁.x[:, v], NH₁.x, ATAmodel.settings.FS.counts, ol_maxᵥ, v)
					if fF == false
						if ATAmodel.settings.opt_type == "MAXIMIN"
							NH₁.obj[v] = eval_TIF_MM_v(x_Iᵥ, IIFᵥ)
						end
					end
					NH₁.f = comp_f(NH₁, opt_feas)
					f_evals+= 1
					if (NH₁.f <= NH₀.f)
						#switch item
						#NH₀.ol = eval_overlap(NH₀.x, FScounts, ATAmodel.settings.ol_max, T, NH₀.ol)
						#NH₀.f = comp_f(NH₀, opt_feas)
						NH₀ = mycopy(NH₁, NH₀)
						# println("better to remove, new f₀ = ", NH₀.f)
						if (NH₁.f < NH⁺.f)
							exit = 1
							convergence = 0
							NH⁺ = mycopy(NH₁, NH⁺)
							if verbosity == 2
								println("f⁺ = ", NH⁺.f, ", t = ", v)
							else
								Printf.@printf "+"
							end
						end
					else
						p = sig_c( - (NH₁.f - NH₀.f) / t) # -(NH₁.f - NH₀.f) / t always negative
						if p<1e-10 p = zeros(Float64)[1]; end
						if (rand() < p)
							#remove
							#exit = 1
							NH₀ = mycopy(NH₁, NH₀)
							#NH₀.ol = eval_overlap(NH₀.x, FScounts, ATAmodel.settings.ol_max, T, NH₀.ol)
							#NH₀.f = comp_f(NH₀, opt_feas)
							if verbosity == 2
								println("SA: f₀ = ", NH₀.f, ", t = ", v)
							else
								Printf.@printf "_"
							end
						end
					end
				end
				#try to switch
				idxᵥ₂ = findall(NH₁.x[:, v] .== zero(Float64)) #Random.shuffle!(findall(NH₁.x[:, v] .== zero(Float64)))#i₂ = 1, ..., I
				betterFound = 0
				i₃ = 0
				x₋₁ = copy(NH₁.x)
				while betterFound == 0 && i₃ < size(idxᵥ₂, 1)
					i₃ += 1
					i₂ = idxᵥ₂[i₃]
					if x_forced0ᵥ[i₂] && i₂ != h
						NH₁ = mycopy(NH₀, NH₁)
						NH₁.x = copy(x₋₁)
						#NH₁.x[h, v] = zero(Float64)
						NH₁.x[i₂, v] = one(Float64)
						# iu = copy(iu₀)
						# iu[i₂]+= 1
						# iu[h]-= 1
						iu = sum(NH₁.x, dims = 2) - ATAmodel.IU.max
						iu = iu[iu .> 0]
						if size(iu, 1) == 0
							NH₁.iu = 0
						else
							NH₁.iu = sum(iu)
						end
						NH₁.infeas[v], x_Iᵥ = check_feas(ATAmodel.settings.FS, Constraintsᵥ, NH₁.x[:, v], nFS, n_items, v)
						if fF == false
							if ATAmodel.settings.opt_type == "MAXIMIN"
								NH₁.obj[v] = eval_TIF_MM_v(x_Iᵥ, IIFᵥ)
							end
						end
						NH₁.ol = eval_overlap(NH₁.x, FScounts, ATAmodel.settings.ol_max, T, NH₁.ol)
						#NH₁.ol[v] = eval_overlapv(NH₁.x[:, v], NH₁.x, ATAmodel.settings.FS.counts, ol_maxᵥ, v)
						NH₁.f = comp_f(NH₁, opt_feas)
						if (NH₁.f <= NH₀.f)
							#switch item
							#NH₀ = mycopy(NH₁, NH₀)
							# NH₀.ol = eval_overlap(NH₀.x, FScounts, ATAmodel.settings.ol_max, T, NH₀.ol)
							# NH₀.f = comp_f(NH₀, opt_feas)
							NH₀ = mycopy(NH₁, NH₀)
							# println("better to switch, new f₀ = ", NH₀.f)
							if (NH₁.f < NH⁺.f)
								exit = 1
								betterFound = 1
								convergence = 0
								NH⁺ = mycopy(NH₁, NH⁺)
								if verbosity == 2
									println("f⁺ = ", NH⁺.f, ", t = ", v)
								else
									Printf.@printf "+"
								end
							end
						else
							p = sig_c(-(NH₁.f - NH₀.f) / t)
							if p<1e-10 p = zeros(Float64)[1]; end
							if (rand() < p)
								#remove
								NH₀ = mycopy(NH₁, NH₀)
								# NH₀.ol = eval_overlap(NH₀.x, FScounts, ATAmodel.settings.ol_max, T, NH₀.ol)
								# NH₀.f = comp_f(NH₀, opt_feas)
								if verbosity == 2
									println("SA: f₀ = ", NH₀.f, ", t = ", v)
								else
									Printf.@printf "_"
								end
							end
						end
					end
				end #end of betterFound (betterFound = 1)
			end #end of itemorder h₂ (exit = 1)

		end #end of testorder v₂ (exit = 1)
		if sum(NH₀.infeas+NH₀.ol)+NH₀.iu<= 0 && fF == true
			fF = false
			for v = 1:T
				x_Iᵥ = FS_to_items(ATAmodel.settings.nFS, ATAmodel.settings.n_items, NH₀.x[:, v], ATAmodel.settings.FS.items)
				if ATAmodel.settings.opt_type == "MAXIMIN"
					NH₀.obj[v] = eval_TIF_MM_v(x_Iᵥ, IIF[v])
				end
			end
			NH₀.f = comp_f(NH₀, opt_feas)
			NH⁺ = mycopy(NH₀, NH⁺)
		end
		f_star[1] = copy(NH₀.f)
		#if exit == 0
		if f_star[2] == f_star[1]
			convergence+= 1
			println(convergence)
		end
		#println("convergence is ", convergence)
		#how many equal f₀ in the last iterations?
		if (f_star[2] == f_star[1] && convergence == conv_max) || time()-start_time >= max_time
			coverage_ok = 1
			#t += 1
		else
			t *= geom_temp
			coverage_ok = 0
			if f_star[1] <= f_star[2]
				f_star[2] = copy(f_star[1])
			end
		end
		if verbosity == 2
			println("\n")
			Printf.@printf("\n  f⁺:	%16.10f", NH⁺.f)
			Printf.@printf("\n  f₀:	%16.10f", NH₀.f)
			Printf.@printf("\n  Convergence:	%16.10f", convergence)
			println("\n")
		end
	end #end of NH coverage (coverage_ok = 1)
	return NH⁺::Neighbourhood
end
