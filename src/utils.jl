function exp_c(x::Float64)
	ccall(:exp, Float64, (Float64,), x)
end

function log_c(x::Float64)
	ccall(:log, Float64, (Float64,), x)
end

function log1p_c(x::Float64)
	ccall(:log1p, Float64, (Float64,), x)
end

function sig_c(x::Float64)
	1 / (1 + exp_c(-x))
end

function sig_cplus(x::Float64)
	1 / ( 1 + exp_c(x))
end

function gemvblas(A::Matrix{Float64}, x::Vector{Float64}, cons::Vector{Float64}, sizex::Int64)
	ccall(("dgemv_64_", "libopenblas64_"), Cvoid, (Ref{UInt8}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ptr{Float64}, Ref{Int64}, Ptr{Float64}, Ref{Int64}, Ref{Float64}, Ptr{Float64}, Ref{Int64}), 'N', size(A, 1), sizex, 1.0, A, max(1, stride(A, 2)), x, 1, -1.0, cons, 1)
	return cons
end

function gemvblasT(A::Matrix{Float64}, x::Vector{Float64}, sizeA2::Int64)
	cons = zeros(Float64, sizeA2)
	ccall(("dgemv_64_", "libopenblas64_"), Cvoid, (Ref{UInt8}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ptr{Float64}, Ref{Int64}, Ptr{Float64}, Ref{Int64}, Ref{Float64}, Ptr{Float64}, Ref{Int64}), 'T', size(A, 1), size(A, 2), 1.0, A, max(1, stride(A, 2)), x, 1, -1.0, cons, 1)
	return cons
end

function combinations(T::Int64)
	comb = vcat([[[t, i] for i = (t + 1) : T] for t = 1 : T]...)
end

function myqle(x::Vector{Float64}, lx::Int64, k::Int64)
	y = ones(k) * Inf
	n = 1
	y[1] = copy(x[1])
	for i = 2 : lx
		current = copy(x[i])
		k2 = 1
		while k2 <= k
			if y[k2] > current
				y[(k2 + 1) : n] .= y[k2 : (n - 1)]
				y[k2] = copy(current)
				k2 = k + 1
				if n < k
					n += 1
				end
			end
			k2 += 1
		end
	end
	return y[k]
end

function mycopy(NH::Neighbourhood, NH_new::Neighbourhood)
	NH_new.x = copy(NH.x)
	NH_new.f = copy(NH.f)
	NH_new.obj = copy(NH.obj)
	NH_new.infeas = copy(NH.infeas)
	NH_new.ol = copy(NH.ol)
	NH_new.iu = copy(NH.iu)
	return NH_new::Neighbourhood
end

function comp_f(NH::Neighbourhood, OptFeas::Float64)
	return (1 - OptFeas) * (sum(NH.infeas + NH.ol) + NH.iu) - OptFeas * (minimum(NH.obj))
end

function myqleSimp(x::Vector{Float64}, ind::Vector{Float64})
	R = size(x, 1)
	alphaR = Int.(ceil.(R .* (ind)))
	alphaR[findall(alphaR .< 1)] .= 1
	alphaR[findall(alphaR .> R)] .= R
	αQle = Inf
	xsorted = sort(x)
	return sort(x)[alphaR]
end

function gemmblas(A::Matrix{Float64}, B::Matrix{Float64}, olMax::Matrix{Float64})
	sa = size(A)
	ol = copy(olMax)
	ccall(("dgemm_64_", "libopenblas64_"), Cvoid,
	(Ref{UInt8}, Ref{UInt8}, Ref{Int64}, Ref{Int64},
	Ref{Int64}, Ref{Float64}, Ptr{Float64}, Ref{Int64},
	Ptr{Float64}, Ref{Int64}, Ref{Float64}, Ptr{Float64},
	Ref{Int64}),
	'T', 'N', sa[2], sa[2],
	sa[1], 1.0, A, max(1, stride(A, 2)),
	B, max(1, stride(B, 2)), - one(Float64), ol,
	max(1, stride(ol, 2)))
	return ol::Matrix{Float64}
end

function cutR(x;
	start = "minimum",
	stop = "maximum",
	nBins = 2,
	returnBreaks = true,
	returnMidPts = false)
	if (start == "minimum") start = minimum(x) end
	if (stop == "maximum") stop = maximum(x) end
	bw=(stop - start) / (nBins - 1)
	midPts = zeros(nBins)
	for i = 1 : nBins
		midPts[i] = start + (i - 1) * bw
	end
	breaks = collect(range(start - (bw / 2); length = nBins + 1, stop = stop + (bw / 2)))
	y = zeros(size(x, 1))
	for j in 1 : size(x, 1)
		for i in 1 : nBins
			if (x[j] >= breaks[i]) && (x[j] < breaks[i + 1])
				y[j] = i
			end
			if i == nBins && x[j] == breaks[i + 1]
				y[j] = i
			end
		end
	end
	if (returnBreaks == true || returnMidPts == true)
		if returnMidPts == false
			return (Int.(y), breaks)
		elseif returnBreaks == false
			return (Int.(y), midPts)
		else
			return (Int.(y), breaks, midPts)
		end
	else
		return Int.(y)
	end
end

function ItemCharFun(pars::DataFrame, theta;
	model = "2PL", #1PL, 2PL, 3PL, grm
	parametrization = "at-b", #"at-b, at-ab, atb, atab"
	D = 1) #true, false
	nItems = size(pars, 1)
	if (model == "1PL")
		(:b in names(pars)) || (:d in names(pars)) || error("discrimination parameter b or d not defined")
		pars.c = zeros(Float64, nItems)
		pars.a = ones(Float64, nItems)
	elseif (model == "2PL")
		(:a in names(pars) || :a in names(pars)) || error("discrimination parameter a1 or a not defined")
		(:b in names(pars)) || (:d in names(pars)) || error("discrimination parameter b or d not defined")
		pars.c = zeros(Float64, nItems)
	elseif (model == "3PL")
		(:a1 in names(pars) || :a in names(pars)) || error("discrimination parameter a1 not defined")
		(:b in names(pars)) || (:d in names(pars)) || error("discrimination parameter b or d not defined")
		(:c in names(pars)) || (:g in names(pars)) || error("discrimination parameter c or g not defined")
	elseif (model == "grm")
		(:a1 in names(pars) || :a in names(pars)) || error("discrimination parameter a1 not defined")
	end

	if :a in names(pars)
		a = pars.a
	else
		a = pars.a1
	end
	if (:c in names(pars))
		c = pars.c
	else
		c = pars.g
	end
	nb = 0
	b = zeros(Float64, size(pars,1), 1)
	for n in names(pars)
		if startswith(string(n), "b")
			nb += 1
			b = hcat(b, pars[!, n])
		end
	end
	nd = 0
	b = zeros(Float64, size(pars,1), 1)
	for n in names(pars)
		if startswith(string(n), "d")
			nd += 1
			d = hcat(d, pars[!, n])
		end
	end
	if nd > 0
		b = mapslices(x -> b .- x, d;  dims = 2)
	end
	I = size(b, 1)
	K = size(theta, 1)
	a2 = ones(I)
	if parametrization == "at-b" #a*theta - b
		b = -b
	elseif parametrization == "at-ab" #a*(theta - b)
		a2 = copy(a)
		b = -b
	elseif parametrization == "atab" #a*theta + b
		a2 = copy(a)
	end
	p = sig_c.([D * (a[i] * theta[k] + a2[i] * b[i, j]) for i = 1 : I, k = 1 : K, j = 1 : nb])
	#pder =  eachslice(((1 .- p) .* p), dims = 1) .* a
	pder =  mapslices(x -> (1 .- x) .* x .* a, p;  dims = 1)
	if model != "grm"
		p = c .+ ((1 .- c) .* p)
	else
		#generalized PCM, does not work
		# p[:, :, 1] = 1 .- p[:, :, 1]
		# pder[:, :, 1] = .- pder[:, :, 1]
		# for k = (nb - 1) :  -1 : 2
		# 	p[:, :, k] = p[:, :, k] - p[:, :, k-1]
		# 	pder[:, :, k] = pder[:, :, k] - pder[:, :, k-1]
		# end
	end
	return p, pder
end

function ItemInfoFun(pars::DataFrame, theta;
	model = "2PL", #1PL, 2PL, 3PL, grm
	parametrization = "at-b", #"at-b, at-ab, atb, atab"
	D = 1) #true, false
	nItems = size(pars, 1)
	if (model == "1PL")
		(:b in names(pars)) || (:d in names(pars)) || error("discrimination parameter b or d not defined")
		pars.c = zeros(Float64, nItems)
		pars.a = ones(Float64, nItems)
	elseif (model == "2PL")
		(:a in names(pars) || :a in names(pars)) || error("discrimination parameter a1 or a not defined")
		(:b in names(pars)) || (:d in names(pars)) || error("discrimination parameter b or d not defined")
		pars.c = zeros(Float64, nItems)
	elseif (model == "3PL")
		(:a1 in names(pars) || :a in names(pars)) || error("discrimination parameter a1 not defined")
		(:b in names(pars)) || (:d in names(pars)) || error("discrimination parameter b or d not defined")
		(:c in names(pars)) || (:g in names(pars)) || error("discrimination parameter c or g not defined")
	elseif (model == "grm")
		(:a1 in names(pars) || :a in names(pars)) || error("discrimination parameter a1 not defined")
	end

	if :a in names(pars)
		a = pars.a
	else
		a = pars.a1
	end
	if (:c in names(pars))
		c = pars.c
	else
		c = pars.g
	end
	nb = 0
	b = zeros(Float64, size(pars,1), 1)
	for n in names(pars)
		if startswith(string(n), "b") || startswith(string(n), "d")
			nb += 1
			b = hcat(b, pars[!, n])
		end
	end
	I = size(b, 1)
	K = size(theta, 1)
	a2 = ones(I)
	if parametrization == "at-b" #a*theta - b
	elseif parametrization == "at-ab" #a*(theta - b)
		a2 = copy(a)
	elseif parametrization == "atab" #a*theta + b
		a2 = copy(a)
	end
	p, pder = ItemCharFun(pars, theta; model = model, parametrization = parametrization)
	nb = size(p, 3)
	if model != "grm"
		i = (a .^ 2) .* ((1 .- p) ./ p) .* ((p .- c) ./ (1 .- c)).^2
	else
		i = pder.^2 ./ p
		i = sum(i, dims = 3)[:, :, 1]
	end
	return i
end

function StudentLikeFun(f::Float64, r::Vector{Float64}, iIndex::Vector{Int64}, a::Vector{Float64}, b::Vector{Float64}, θ::Float64;
	logL = false)
	I = size(a, 1)
	likel = 0
	for i in iIndex
		p = b[i] + a[i] * θ
		likel += r[i] * p - log1p_c(exp_c(p))
	end
	return likel::Float64
end

function genResp(f::Vector{Float64}, pars::DataFrame, θ::Vector{Float64}, design::Matrix{Float64};
	model = "2PL",
	method = "classicUniform"
	)
	if (model == "1PL")
		(:b in names(pars)) || (:d in names(pars)) || error("discrimination parameter b or d not defined")
		pars.g = fill(0, size(pars, 1))
		pars.a1 = fill(1, size(pars, 1))
	end
	if (model == "2PL")
		(:a1 in names(pars) || :a in names(pars)) || error("discrimination parameter a1 not defined")
		(:b in names(pars)) || (:d in names(pars)) || error("discrimination parameter b or d not defined")
		pars.g = fill(0, size(pars, 1))
	end

	if (model == "3PL")
		(:a1 in names(pars) || :a in names(pars)) || error("discrimination parameter a1 not defined")
		(:b in names(pars)) || (:d in names(pars)) || error("discrimination parameter b or d not defined")
		(:c in names(pars)) || (:g in names(pars)) || error("discrimination parameter c or g not defined")
	end
	if :a in names(pars)
		a = pars.a
	else
		a = pars.a1
	end
	if (:c in names(pars))
		c  =pars.c
	else
		c = pars.g
	end
	if (:b in names(pars))
		b = pars.b
	else
		b = pars.d
	end
	I = size(b, 1)
	N = size(θ, 1)
	nindex = Array{Array{Int64, 1}, 1}(undef, I)
	iindex = Array{Array{Int64, 1}, 1}(undef, N)
	for n = 1 : N
		iindex[n] = findall(design[:, n] .== 1.0)
		if n <= I
			nindex[n] = findall(design[n, :] .== 1.0)
		end
	end
	p = Matrix(undef, I, N)
	lp = Matrix(undef, I, N)
	lq = Matrix(undef, I, N)
	for i = 1 : I
		for n = 1 : N
			pr = c[i] + ((1 - c[i]) * (1 / (1 + exp(-(b[i] + a[i] * θ[n])))))
			p[i, n] = pr
			lp[i, n] = log(pr)
			lq[i, n] = log(1 - pr)
		end
	end
	resp = Matrix(undef, I, N)
	if method == "cumulatedPersons"
		@fastmath @inbounds  for i = 1 : I
			gapScore = 3
			while gapScore >= 2
				p2 = p[i, :]
				p2 = hcat((1 .- p2) , p2)#2
				unif = rand(Uniform(0, 1), N)#5
				n = 1#6
				while n <= N#7
					csum = p2[n, 1]#8
					cat = 0#9
					if design[i, n] == 0
						resp[i, n] = 0
						n = n + 1
					else
						while csum < unif[n]
							cat = cat+1
							if (cat == 2) break end
							csum = csum + p2[n, (cat + 1)]
						end
						resp[i, n] = cat
						n = n + 1
					end
				end
				gapScore = abs(sum(resp[i, :]) - sum(p[i, n] for n in nindex[i]))
			end
		end
	end
	if method == "cumulatedItems"
		@fastmath @inbounds  for n = 1 : N#4
			#gapScore=3
			#while gapScore>=2
			p2 = p[:, n]
			p2 = hcat((1 .- p2) , p2)#2
			unif = rand(Uniform(0, 1), I)#5
			samplei = sample(collect(1 : I), I, replace = false)
			i = 1#6
			while i <= I#7
				csum = p2[samplei[i], 1]#8
				cat = 0#9
				if design[samplei[i], n] == 0
					resp[samplei[i], n] = 0#missing
					i = i + 1
				else
					while csum < unif[samplei[i]]
						cat = cat + 1
						if (cat == 2) break end
						csum = csum + p2[samplei[i], cat + 1]
					end
					resp[samplei[i], n]=cat
					i = i + 1
				end
			end
			#gapScore=abs(sum(skipmissing(resp[:,n]))-sum(p[i,n] for i in iindex[n]))
			#if gapScore>=0.5
			#	println("person ",n," gap=",gapScore)
			#end
			#end
		end
	end
	if method == "classicUniform"
		for n = 1 : N#4
			#gapScore=2
			#while gapScore>=0.5
			unif = rand(Uniform(0, 1), I)#5
			samplei = sample(collect(1 : I), I, replace = false)
			i = 1#6
			while i <= I#7
				if design[samplei[i], n] == 0
					resp[samplei[i], n] = 0#missing
					i = i + 1
				else
					if unif[samplei[i]] < p[samplei[i], n]
						resp[samplei[i], n] = 1
					else
						resp[samplei[i], n] = 0
					end
					i = i + 1
				end
			end
			#gapScore=abs(sum(skipmissing(resp[:,n]))-sum(p[i,n] for i in iindex[n]))
			#if gapScore>=0.5
			#	println("person ",n," gap=",gapScore)
			#end
			#end
		end
	end
	if method == "cumItemsPattern"
		for n = 1 : N#4+
			println("start person ", n)
			iIndex = iindex[n]
			IIndex = size(iIndex, 1)
			p2 = p[iIndex, n]
			p2 = hcat((1 .- p2) , p2)#2
			patterns = Vector{Vector{Int64}}(undef, 1)
			nPattern = Vector{Int64}(undef, 1)
			for r = 1:1000
				respn = Vector{Int64}(undef, IIndex)
				unif = rand(Uniform(0, 1), IIndex)#5
				samplei=sample(collect(1 : IIndex), IIndex, replace=false)
				for i = 1 : IIndex
					csum = p2[samplei[i], 1]#8
					cat = 0#9
					while csum < unif[samplei[i]]
						cat = cat + 1
						if (cat == 2) break end
						csum = csum + p2[samplei[i], cat + 1]
					end
					respn[i] = cat
				end
				if r == 1
					patterns[1] = respn
					nPattern[1] = 1
				else
					println(size(patterns, 1))
					corr = 0
					for pat = 1 : size(patterns, 1)
						if patterns[pat] == respn
							corr = pat
						end
					end
					if corr == 0
						push!(patterns, respn)
						nPattern = hcat(nPattern, 1)
					else
						nPattern[corr] = nPattern[corr] + 1
					end
				end

			end
			nPatterns = size(patterns, 1)
			probPatterns = nPattern ./ 1000
			println(probPatterns)
			resp[iIndex, n] .= sample(patterns, pweights(probPatterns), 1)
			println("end person ", n)
		end
	end
	if method == "classicUniformPattern"
		for n = 1 : N#4+
			println("start person ", n)
			iIndex = iindex[n]
			IIndex = size(iIndex, 1)
			patterns = Vector{Vector{Int64}}(undef, 1)
			nPattern = Vector{Int64}(undef, 1)
			for r = 1 : 1000
				respn = Vector{Int64}(undef, IIndex)
				unif = rand(Uniform(0, 1) , IIndex)#5
				samplei = sample(collect(1 : IIndex), IIndex, replace = false)
				for  i = 1 : IIndex#7
					if unif[samplei[i]] < p[samplei[i], n]
						respn[i] = 1
					else
						respn[i] = 0
					end
				end
				if r == 1
					patterns[1] = respn
					nPattern[1] = 1
				else
					println(size(patterns, 1))
					corr=0
					for pat = 1:size(patterns, 1)
						if patterns[pat] == respn
							corr = pat
						end
					end
					if corr == 0
						push!(patterns, respn)
						nPattern = hcat(nPattern, 1)
					else
						nPattern[corr] = nPattern[corr] + 1
					end
				end
			end
			nPatterns = size(patterns, 1)
			probPatterns = nPattern ./ 1000
			resp[iIndex, n] .= sample(patterns, pweights(probPatterns), 1)
			println("end person ", n)
		end
	end
	return resp
end
