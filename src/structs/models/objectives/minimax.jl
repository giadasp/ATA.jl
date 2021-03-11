mutable struct MinimaxObjectiveCore
    points::Vector{Float64}
    targets::Vector{Float64}
    IIF::Matrix{Float64}

    MinimaxObjectiveCore() = new(Vector{Float64}(), Vector{Float64}(), Matrix{Float64}(undef,0,0))
end

mutable struct MinimaxObjective <: AbstractObjective
    name::String
    sense::String
    cores::Vector{MinimaxObjectiveCore}

    MinimaxObjective(name, sense, cores) = new("MINIMAX", "min", cores)
    MinimaxObjective() = new("MINIMAX", "min", MinimaxObjectiveCore[])
end