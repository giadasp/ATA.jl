mutable struct RobustMaximinObjectiveCore
    points::Vector{Float64}
    standard_deviation::Array{Float64,2}
    IIF::Array{Float64,2}
    RobustMaximinObjectiveCore() =
        new(Vector{Float64}(), Array{Float64,2}(undef, 0, 0), Array{Float64,2}(undef, 0, 0))
end
