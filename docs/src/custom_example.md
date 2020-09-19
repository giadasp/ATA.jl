
## Custom objective function

To add a custom objective function modify step 8 in previous examples.
It works only with the siman solver.

8\. Add objective function (Optional)

```julia
@info add_obj_fun!(ATAmodel)[2] 

ATAmodel.obj.type = "custom"

ATAmodel.obj.fun = function (x::Matrix{Float64}, obj_args::NamedTuple)
	IIF = obj_args.IIF
	T = size(IIF, 1)
	TIF = zeros(Float64, T)
	for t = 1:T
        K, I = size(IIF[t])
        # ungroup items
        xₜ = FS_to_items(x[:, t], obj_args.FS_items)
		if K > 1
			TIF[t] = Inf
			for k = 1:K
				TIF[t] = min(TIF, LinearAlgebra.dot(IIF[1, :], xₜ))
			end
		else
			TIF[t] = LinearAlgebra.BLAS.gemv('N', IIF, xₜ)[1]
		end
    end
    min_TIF = minimum(TIF)
    TIF = [min_TIF for t = 1:T] 
    # Must return a vector of length T.
    # The resulting objective function is the minimum of all values in this vector.
	return TIF::Vector{Float64}
end

ATAmodel.obj.args = (IIF = FileIO.load("data/IIF.jld2", "IIF"), FS_items = ATAmodel.settings.FS.items)
```