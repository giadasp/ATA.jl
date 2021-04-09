
function eval_TIFₜ(xₜ::Vector{Float64}, coreₜ::MaximinObjectiveCore) # x = I
    return minimum(LinearAlgebra.BLAS.gemv('N', coreₜ.IIF, xₜ))::Float64
end

function eval_TIFₜ(xₜ::Vector{Float64}, coreₜ::CCMaximinObjectiveCore) # x = I
    IIF = coreₜ.IIF
    α = coreₜ.alpha
    K, I, R = size(IIF)
    alphaR = Int(ceil(R * (α)))
    αQle = Inf
    if K > 1
        for k = 1:K
            αQle = min(αQle, sort(LinearAlgebra.BLAS.gemv('T', IIF[k, :, :], xₜ))[alphaR])
        end
    else
        αQle = sort(_gemvblasT(IIF[1, :, :], xₜ, R))[alphaR]
    end
    return αQle::Float64
end

function eval_TIFₜ(xₜ::Vector{Float64}, coreₜ::MinimaxObjectiveCore) # x = I
    IIF = coreₜ.IIF
    targets = coreₜ.targets
    return -maximum(abs.(LinearAlgebra.BLAS.gemv('N', IIF, xₜ) - targets))::Float64
end


