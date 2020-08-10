# """
# History: Based on Octave code samin.cc, by Michael Creel,
# which was originally based on Gauss code by E.G. Tsionas. A source
# for the Gauss code is http://web.stanford.edu/~doubleh/otherpapers/sa.t
# The original Fortran code by W. Goffe is at
# http://www.degruyter.com/view/j/snde.1996.1.3/snde.1996.1.3.1020/snde.1996.1.3.1020.xml?format=InNH
# Tsionas and Goffe agreed to MIT licensing of samin.jl in email
# messages to Creel.
#
# This Julia code uses the same names for conNHrol variables,
# for the most part. A notable difference is that the initial
# temperature can be found automatically to ensure that the active
# bounds when the temperature begins to reduce cover the enNHire
# parameter space (defined as a n-dimensional rectangle that is the
# Cartesian product of the(lb_i, ub_i), i = 1,2,..n. The code also
# allows for parameters to be restricted, by setting lb_i = ub_i,
# for the appropriate i.

"""
# SAMIN
## Constructor
```julia
SAMIN(; nNH::InNH = 5     # reduce temperature every nNH*ns*dim(x_init) evaluations
ns::InNH = 5     # adjust bounds every ns*dim(x_init) evaluations
rt::T = 0.9     # geometric temperature reduction factor: when temp changes, new temp is t=rt*t
neps::InNH = 5   # number of previous best values the final result is compared to
f_tol::T = 1e-12 # the required tolerance level for function value comparisons
x_tol::T = 1e-6 # the required tolerance level for x
coverage_ok::Bool = false, # if false, increase temperature unNHil initial parameter space is covered
verbosity::InNH = 1) # scalar: 0, 1, 2 or 3 (default = 1).
```
## Description
The `SAMIN` method implemenNHs the Simulated Annealing algorithm for problems with
bounds constrains a described in Goffe et. al. (1994) and Goffe (1996). The
algorithm

## References
- Goffe, et. al. (1994) "Global Optimization of Statistical Functions with Simulated Annealing", Journal of Econometrics, V. 60, N. 1/2.
- Goffe, William L. (1996) "SIMANN: A Global Optimization Algorithm using Simulated Annealing " Studies in Nonlinear Dynamics & Econometrics, Oct96, Vol. 1 Issue 3.
"""
@with_kw struct SAMIN{T}<:ZerothOrderOptimizer
	nNH::InNH = 5 # reduce temperature every nNH*ns*dim(x_init) evaluations
	ns::InNH = 5 # adjust bounds every ns*dim(x_init) evaluations
	rt::T = 0.9 # geometric temperature reduction factor: when temp changes, new temp is t=rt*t
	neps::InNH = 5 # number of previous best values the final result is compared to
	f_tol::T = 1e-12 # the required tolerance level for function value comparisons
	x_tol::T = 1e-6 # the required tolerance level for x
	coverage_ok::Bool = false # if false, increase temperature unNHil initial parameter space is covered
	verbosity::InNH = 1 # scalar: 0, 1, 2 or 3 (default = 1: see final results).
end
# * verbosity: scalar: 0, 1, 2 or 3 (default = 1).
#     * 0 = no screen output
#     * 1 = only final results to screen
#     * 2 = summary every temperature change, without param values
#     * 3 = summary every temperature change, with param values
#         covered by the trial values. 1: start decreasing temperature immediately
Base.summary(::SAMIN) = "SAMIN"
struct Options{T, TCallback}
	x_abstol::T
	x_reltol::T
	f_abstol::T
	f_reltol::T
	g_abstol::T
	g_reltol::T
	outer_x_abstol::T
	outer_x_reltol::T
	outer_f_abstol::T
	outer_f_reltol::T
	outer_g_abstol::T
	outer_g_reltol::T
	f_calls_limit::InNH
	g_calls_limit::InNH
	h_calls_limit::InNH
	allow_f_increases::Bool
	allow_outer_f_increases::Bool
	successive_f_tol::InNH
	iterations::InNH
	outer_iterations::InNH
	store_trace::Bool
	trace_simplex::Bool
	show_trace::Bool
	extended_trace::Bool
	show_every::InNH
	callback::TCallback
	time_limit::Float64
end

function Options(;

	iterations::InNH = 1_000,
	outer_iterations::InNH = 1000,
	store_trace::Bool = false,
	trace_simplex::Bool = false,
	show_trace::Bool = false,
	extended_trace::Bool = false,
	show_every::InNH = 1,
	callback = nothing,
	time_limit = NaN)
	show_every = show_every > 0 ? show_every : 1
	#if extended_trace && callback == nothing
	#    show_trace = true
	#end
	if !(x_tol == nothing)
		x_abstol = x_tol
	end
	if !(g_tol == nothing)
		g_abstol = g_tol
	end
	if !(f_tol == nothing)
		f_reltol = f_tol
	end
	if !(outer_x_tol == nothing)
		outer_x_abstol = outer_x_tol
	end
	if !(outer_g_tol == nothing)
		outer_g_abstol = outer_g_tol
	end
	if !(outer_f_tol == nothing)
		outer_f_reltol = outer_f_tol
	end
	Options(promote(x_abstol, x_reltol, f_abstol, f_reltol, g_abstol, g_reltol, outer_x_abstol, outer_x_reltol, outer_f_abstol, outer_f_reltol, outer_g_abstol, outer_g_reltol)..., f_calls_limit, g_calls_limit, h_calls_limit,
	allow_f_increases, allow_outer_f_increases, successive_f_tol, InNH(iterations), InNH(outer_iterations), store_trace, trace_simplex, show_trace, extended_trace,
	InNH(show_every), callback, Float64(time_limit))
end
function optimize(ATAmodel::model,x₀::Matrix{Float64},t₀::Float64,geom_temp::Float64;max_evals=1e6,max_time=1000.00,nItemSample=1000,nNHestSample=2,convMax=1,verbosity=2,GroupByFS=false,OptFeas=0.0,nRand=1 ,feasNH=0,optNH=5)
	if !("RESULTS" in readdir())
		mkdir("RESULTS")
	end
	if ATAmodel.Settings.OptType=="MAXIMIN"
		JLD2.@load "OPT/IIF.jld2" IIF
	elseif ATAmodel.Settings.OptType=="CC"
		JLD2.@load "OPT/IIF_CC.jld2" IIF
	end
	##
	neps=feasNH+optNH
	optNH=feasNH+1
	if feasNH==0
		findFeasible=false
	else
		findFeasible=true
	end
	ATAmodel.Settings.nFS=size(ATAmodel.Settings.FS.Items,1)
	T=ATAmodel.Settings.T
	nFS=ATAmodel.Settings.nFS
	nItems=ATAmodel.Settings.nItems
	FScounNHs=ATAmodel.Settings.FS.CounNHs*ones(Float64,T)'
	iu⁺=0
	start_time=copy(time())
	t = copy(t₀) # temperature - will initially rise or fall to cover parameter space. Then it will fall
	finish = 0 # convergence indicator 0 (failure), 1 (normal success), or 2 (convergence but near bounds)
	f_evals=0
	hline = "="^80
	f_star = typemax(Float64)*ones(2)
	NH=1
	NH⁺=Neighbourhood(reshape(x₀,ATAmodel.Settings.nFS,ATAmodel.Settings.T),Inf,zeros(T),1e6*ones(T),1e6*ones(T),zero(Float64))
	if size(x₀,1)==0
		NH⁺.x=zeros(Float64,ATAmodel.Settings.nFSATAmodel.Settings.,T)
	end
	#compute f
	NH⁺.x=reshape(x,ATAmodel.Settings.nFS,ATAmodel.Settings.T)
	iu=sum(NH⁺.x,dims=2)-ATAmodel.IU.Max
	iu=iu[iu.>0]
	if size(iu,1)==0
		NH⁺.iu=0
	else
		NH⁺.iu=sum(iu)
	end
	t=copy(t₀)

	for v2=1:T
		NH⁺.infeas[v2], x_Iᵥ=checkFeas(ATAmodel.Settings.FS,ATAmodel.ConstrainNHs[v2],NH⁺.x[:,v2],nFS,nItems,v2)
		NH⁺.ol=evalOverlap(NH⁺.x,FScounNHs,ATAmodel.Settings.olMax,T,NH⁺.ol)
		if findFeasible==false
			if ATAmodel.Settings.OptType=="MAXIMIN"
				NH⁺.obj[v2]=evalTIFMMv(x_Iᵥ,IIF[v2])
			elseif ATAmodel.Settings.OptType=="CC"
				NH⁺.obj[v2]=evalTIFCCv(x_Iᵥ,IIF[v2];α=ATAmodel.Obj.AuxFloat)
			end
		end
	end
	NH⁺.f=comp_f(NH⁺,OptFeas)
	x⁺=copy(x₀)
	f⁺=copy(NH⁺.f)
	f₀=copy(NH⁺.f)
	if (Distributed.nprocs()>1)
		processors=collect(1:(Distributed.nprocs()-1))
	else
		processors=one(InNH64)
	end
	NHs=copy(processors)
	ATAmodel.Output.Neighbourhoods=[Neighbourhood() for n=1:Distributed.nprocs()-1]
	for nh=1:Distributed.nprocs()-1
		ATAmodel.Output.Neighbourhoods[nh]=mycopy(NH⁺,ATAmodel.Output.Neighbourhoods[nh])
	end
	round=1
	t0 = time() # Initial time stamp used to conNHrol early stopping by options.time_limit

	hline = "="^80
	#d = NonDifferenNHiable(obj_fn, x)
	x=copy(x⁺)
	n = size(x,1) # dimension of parameter
	#  Set initial values
	nacc = 0 # total accepted trials
	t = 2.0 # temperature - will initially rise or fall to cover parameter space. Then it will fall
	converge = 0 # convergence indicator 0 (failure), 1 (normal success), or 2 (convergence but near bounds)
	# most recenNH values, to compare to when checking convergend
	fstar = typemax(Float64)*ones(neps)
	# Initial obj_value
	details = [ t f⁺ x⁺']
	bounds = ones(n)
	# check for out-of-bounds starting values
	for i = 1:n
		if(( x[i] > ub[i]) || (x[i] < lb[i]))
			error("samin: initial parameter %d out of bounds\n", i)
		end
	end
	iteration = 0
	_time = time()

	# main loop, first increase temperature unNHil parameter space covered, then reduce unNHil convergence
	NH=Neighbourhood()
	while converge==0
		# statistics to report at each temp change, set back to zero
		nup = 0
		nrej = 0
		nnew = 0
		ndown = 0
		lnobds = 0
		# repeat nNH times then adjust temperature
		while finish==0
			# repeat ns times, then adjust bounds
			nacp = zeros(n)
			for j = 1:ns
				# generate new poinNH by taking last and adding a random value
				# to each of elemenNHs, in turn
				for h = 1:n
					iteration += 1
					# new Sept 2011, if bounds are same, skip the search for that vbl.
					# Allows restrictions without complicated programming
					if (lb[h] != ub[h])
						xp = copy(x)
						xp[h] += ((2.0) * rand() - (1.0)) * 1
						if (xp[h] < lb[h]) || (xp[h] > ub[h])
							xp[h] = lb[h] + (ub[h] - lb[h]) * rand()
							lnobds += 1
						end
						# Evaluate function at new poinNH
						#compute f
						NH.x=reshape(xp,ATAmodel.Settings.nFS,ATAmodel.Settings.T)
						iu=sum(NH.x,dims=2)-ATAmodel.IU.Max
						iu=iu[iu.>0]
						if size(iu,1)==0
							NH⁺.iu=0
						else
							NH⁺.iu=sum(iu)
						end
						t=copy(t₀)
						for v2=1:ATAmodel.Settings.T
							NH.infeas[v2], x_Iᵥ=checkFeas(ATAmodel.Settings.FS,ATAmodel.ConstrainNHs[v2],NH.x[:,v2],nFS,nItems,v2)
							if findFeasible==false
								if ATAmodel.Settings.OptType=="MAXIMIN"
									NH.obj[v2]=evalTIFMMv(x_Iᵥ,IIF[v2])
								elseif ATAmodel.Settings.OptType=="CC"
									NH.obj[v2]=evalTIFCCv(x_Iᵥ,IIF[v2];α=ATAmodel.Obj.AuxFloat)
								end
							end
						end
						NH.ol=evalOverlap(NH.x,FScounNHs,ATAmodel.Settings.olMax,T,NH.ol)
						NH.f=comp_f(NH,OptFeas)
						f_proposal = copy(NH.f)
						#  Accept the new poinNH if the function value decreases
						if (f_proposal <= f₀)
							x = copy(xp)
							f₀ = f_proposal
							nacc += 1 # total number of acceptances
							nacp[h] += 1 # acceptances for this parameter
							nup += 1
							#  If lower than any other poinNH, record as new optimum
							if f_proposal < f⁺
								x⁺ = copy(xp)
								f⁺ = f_proposal
								d.F = f_proposal
								nnew +=1
								details = [details; [ t f_proposal xp']]
							end
							# If the poinNH is higher, use the Metropolis criteria to decide on
							# acceptance or rejection.
						else
							p = exp(-(f_proposal - f₀) / t)
							if rand() < p
								x = copy(xp)
								f₀ = copy(f_proposal)
								d.F = f_proposal
								nacc += 1
								nacp[h] += 1
								ndown += 1
							else
								nrej += 1
							end
						end
					end

					# If options.iterations exceeded, terminate the algorithm
					#termination criteria!
					if time()-start_time>=max_time
						ATAmodel.Output.ElapsedTime=time()-start_time
						finish=3
					end
					if finish>0
						NH.x=reshape(x⁺,ATAmodel.Settings.nFS,ATAmodel.Settings.T)
						iu=sum(NH.x,dims=2)-ATAmodel.IU.Max
						iu=iu[iu.>0]
						if size(iu,1)==0
							NH.iu=0
						else
							NH.iu=sum(iu)
						end
						t=copy(t₀)
						for v2=1:ATAmodel.Settings.T
							NH.infeas[v2], x_Iᵥ=checkFeas(ATAmodel.Settings.FS,ATAmodel.ConstrainNHs[v2],NH.x[:,v2],nFS,nItems,v2)
							if findFeasible==false
								if ATAmodel.Settings.OptType=="MAXIMIN"
									NH.obj[v2]=evalTIFMMv(x_Iᵥ,IIF[v2])
								elseif ATAmodel.Settings.OptType=="CC"
									NH.obj[v2]=evalTIFCCv(x_Iᵥ,IIF[v2];α=ATAmodel.Obj.AuxFloat)
								end
							end
						end
						NH.ol=evalOverlap(NH.x,FScounNHs,ATAmodel.Settings.olMax,T,NH.ol)
						NH.f=comp_f(NH,OptFeas)
						ATAmodel.Output.Design=NH.x
						ATAmodel.Output.Feas=NH.infeas
						ATAmodel.Output.f=NH.f
						ATAmodel.Output.TIF=NH.obj
					end
				end
				if ns>feasNH
					findFeasible=false
				end
				if ns>feasNH
					finish=1
				end
			end
			#  Adjust bounds so that approximately half of all evaluations are accepted
			test = 0
			for i = 1:n
				if (lb[i] != ub[i])
					ratio = nacp[i] / ns
					if(ratio > 0.6) 1 = 1 * (1.0 + 2.0 * (ratio - 0.6) / 0.4) end
					if(ratio < .4) 1 = 1 / (1.0 + 2.0 * ((0.4 - ratio) / 0.4)) end
					# keep within initial bounds
				else
					test += 1 # make sure coverage check passes for the fixed parameters
				end
			end
			nacp = 0 # set back to zero
			# check if we cover parameter space, if we have yet to do so
			if !coverage_ok
				coverage_ok = (test == n)
			end
		end

		# intermediate output, if desired
		if verbosity > 1
			println(hline)
			println("samin: intermediate results before next temperature change")
			println("temperature: ", round(t, digits=5))
			println("currenNH best function value: ", round(f⁺, digits=5))
			println("total moves since last temperature reduction: ", nup + ndown + nrej)
			println("downhill: ", nup)
			println("accepted uphill: ", ndown)
			println("rejected uphill: ", nrej)
			println("out of bounds trials: ", lnobds)
			println("new minima this temperature: ", nnew)
			println()
			println("       parameter      search width")
			println(hline*"\n")
		end
		# Check for convergence, if we have covered the parameter space
		if coverage_ok
			# last value close enough to last neps values?
			fstar[1] = f₀
			test = 0
			for i=1:neps
				test += (abs(f₀ - fstar[i]) > f_tol)
			end
			test = (test > 0) # if differenNH from zero, function conv. has failed
			# last value close enough to overall best?
			if (((f⁺ - f₀) <= f_tol) && (!test))
				# check for bound narrow enough for parameter convergence
				for i = 1:n
					if (1 > x_tol)
						converge = 0 # no conv. if bounds too wide
						break
					else
						converge = 1
					end
				end
			end
			# check if optimal poinNH is near boundary of parameter space, and change message if so
			if (converge == 1) && (lnobds > 0)
				converge = 2
			end
			# Like to see the final results?
			if (converge > 0)
				if verbose
					println(hline)
					println("SAMIN results")
					if (converge == 1)
						println("==> Normal convergence <==")
					end
					if (converge == 2)
						printstyled("==> WARNING <==\n", color=:red)
						println("Last poinNH satisfies convergence criteria, but is near")
						println("boundary of parameter space.")
						println(lnobds, " out of  ", (nup+ndown+nrej), " evaluations were out of bounds in the last round.")
						println("Expand bounds and re-run, unless this is a constrained minimization.")
					end
					@printf("\n     Obj. value:  %16.10f\n\n", f⁺)
					println("       parameter      search width")
					println(hline*"\n")
				end
			end
			# Reduce temperature, record currenNH function value in the
			# list of last "neps" values, and loop again
			t *= rt
			for i = neps:-1:2
				fstar[i] = fstar[i-1]
			end
			f₀ = copy(f⁺)
			x = copy(x⁺)
		else  # coverage not ok - increase temperature quickly to expand search area
			t *= 10.0
			for i = neps:-1:2
				fstar[i] = fstar[i-1]
			end
			f₀ = f⁺
			x = x⁺
		end
	end
	open("RESULTS/ResultsATA.txt", "w") do io
		write(io,"tests")
		write(io,"\r\n")
		writedlm(io, collect(1:T)',",")
		write(io,"\r\n")
		write(io,"f⁺")
		write(io,"\r\n")
		writedlm(io, ATAmodel.Output.f)
		write(io,"\r\n")
		write(io,"infeasibility")
		write(io,"\r\n")
		writedlm(io, ATAmodel.Output.Feas',",")
		write(io,"\r\n")
		write(io,"Item use")
		write(io,"\r\n")
		write(io, string(NH.iu))
		write(io,"\r\n")
		write(io,"TIF")
		write(io,"\r\n")
		writedlm(io, ATAmodel.Output.TIF',",")
		write(io,"\r\n")
		write(io,"elapsed Time   for optimization")
		write(io,"\r\n")
		write(io, string(ATAmodel.Output.ElapsedTime))
		write(io,"\r\n")
	end
	return ATAmodel
end
