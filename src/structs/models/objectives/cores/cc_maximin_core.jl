mutable struct CCMaximinObjectiveCore
    alpha::Float64
    points::Vector{Float64}
    IIF::Array{Float64,3}
    CCMaximinObjectiveCore() =
        new(1, 0.05, Vector{Float64}(), Array{Float64,3}(undef, 0, 0, 0))
end