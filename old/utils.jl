function mycompare(a, b)::Cint
   return (a < b) ? -1 : ((a > b) ? +1 : 0)
end
mycompare_c = @cfunction(mycompare, Cint, (Ref{Cdouble}, Ref{Cdouble}));
function sort_c(x::Vector{Float64},lx::Int64)
	ccall(:qsort, Cvoid, (Ptr{Cdouble}, Csize_t, Csize_t, Ptr{Cvoid}), x, lx, sizeof(eltype(x)), mycompare_c)
	return x
end
function fast_select_pivot!(v::Vector{Float64}, lo::Int64, hi::Int64)
	@inbounds begin
		mi = (lo+hi)>>>1
		# sort v[mi] <= v[lo] <= v[hi] such that the pivot is immediately in place
		if v[lo]<v[mi]
			v[mi], v[lo] = v[lo], v[mi]
		end
		if v[hi]< v[lo]
			if v[hi]< v[mi]
				v[hi], v[lo], v[mi] = v[lo], v[mi], v[hi]
			else
				v[hi], v[lo] = v[lo], v[hi]
			end
		end
		# return the pivot
		return v[lo]
	end
end
function fast_partition!(v::Vector{Float64}, lo::Int64, hi::Int64)
	pivot = fast_select_pivot!(v, lo, hi)
	# pivot == v[lo], v[hi] > pivot
	i, j = lo, hi
	@inbounds while true
		i += 1; j -= 1
		while (v[i]< pivot); i += 1; end;
		while (pivot< v[j]); j -= 1; end;
		i >= j && break
		v[i], v[j] = v[j], v[i]
	end
	v[j], v[lo] = pivot, v[j]
	return j
end
function fast_sort!(v::Vector{Float64},lo::Int64,hi::Int64,k::Int64)
	@inbounds while lo < hi
		j = fast_partition!(v, lo, hi)
		if j >= k
			# we don't need to sort anything bigger than j
			hi = j-1
		elseif j-lo < hi-j
			# recurse on the smaller chunk
			# this is necessary to preserve O(log(n))
			# stack space in the worst case (rather than O(n))
			lo < (j-1) && fast_sort!(v, lo, j-1, k)
			lo = j+1
		else
			(j+1) < k && fast_sort!(v, j+1, hi, k)
			hi = j-1
		end
	end
	return v[k]
end
function myqle(x::Vector{Float64},lx::Int64,k::Int64)
	y=zeros(lx)
	lo=copy(x[1])
	k2=1
	for i=1:lx-1
		if x[i+1]<=lo
			lo=copy(x[i+1])
			y[k2]=copy(lo)
			println("new lo ",lo)
			k2+=1
		end
	end
	return y
end

# if method=="MIP"
# 	for n=1:N
# 		iIndex=iindex[n]
# 		p=[c[i]+((1-c[i])*(1 / (1 + exp(-a[i]*(θ[n]-b[i]))))) for i in iIndex]
# 		p[p.==0].=0.00001
# 		m=Model(solver=CplexSolver(CPX_PARAM_PREIND=0,CPX_PARAM_MIPEMPHASIS=1))
# 		@variable(m, x[i=1:size(iIndex,1)], Bin)
# 		@objective(m, Max, sum((x[i]*log(p[i]) + ((f[i] - x[i]) * log(1-p[i]))) for  i=1:size(iIndex,1)))
# 		@constraint(m, sum(x[i] for i=1:size(iIndex,1))-sum(p) <=+1)
# 		@constraint(m, sum(x[i] for i=1:size(iIndex,1))-sum(p) >=-1)
#
# 		#@constraint(m, [n=1:(size(nindex[i],1)-1)], x[n] <= x[n+1] )
# 		solve(m)
# 		xopt=getvalue(x)
# 		resp[:,n].=missing
# 		resp[iindex[n],n].=Int.(xopt)
# 		# m=1
# 		# for n=1:N
# 		#  		if design[i,n]==0
# 		#  			resp[i,n]=missing
# 		#  		else
# 		#  			resp[i,n]=Int(xopt[m])
# 		# 			m=m+1
# 		#  		end
# 		#  end
# 	end
# end
