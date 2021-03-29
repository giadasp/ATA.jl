mutable struct MaximinObjectiveCore
    points::Vector{Float64}
    IIF::Matrix{Float64}

    MaximinObjectiveCore() = new(Vector{Float64}(), Matrix{Float64}(undef, 0, 0))
end
