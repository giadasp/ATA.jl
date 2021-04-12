mutable struct RobustMaximinObjective <: AbstractObjective
    name::String
    sense::String
    R::Int64
    Gamma::Int64
    cores::Vector{RobustMaximinObjectiveCore}

    RobustMaximinObjective(
        R::Int64,
        Gamma::Int64,
        cores::Vector{RobustMaximinObjectiveCore},
    ) = new("robust_maximin", "max", R, Gamma, cores)
    RobustMaximinObjective(cores::Vector{RobustMaximinObjectiveCore}) =
        new("robust_maximin", "max", 1, 1, cores)
    RobustMaximinObjective() =
        new("robust_maximin", "max", 1, 1, RobustMaximinObjectiveCore[])
end
