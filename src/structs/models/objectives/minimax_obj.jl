mutable struct MinimaxObjective <: AbstractObjective
    name::String
    sense::String
    cores::Vector{MinimaxObjectiveCore}

    MinimaxObjective(cores) = new("minimax", "min", cores)
    MinimaxObjective() = new("minimax", "min", MinimaxObjectiveCore[])
end
