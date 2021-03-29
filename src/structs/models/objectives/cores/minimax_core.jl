mutable struct MinimaxObjectiveCore
    points::Vector{Float64}
    targets::Vector{Float64}
    IIF::Matrix{Float64}

    MinimaxObjectiveCore() =
        new(Vector{Float64}(), Vector{Float64}(), Matrix{Float64}(undef, 0, 0))
end