mutable struct DeJongMaximinObjective <: AbstractObjective
    name::String
    sense::String
    R::Int64
    cores::Vector{MaximinObjectiveCore}

    DeJongMaximinObjective(R::Int64, cores::Vector{MaximinObjectiveCore}) =
        new("de_jong_maximin", "max", R, cores)
    DeJongMaximinObjective(cores::Vector{MaximinObjectiveCore}) =
        new("de_jong_maximin", "max", 1, cores)
    DeJongMaximinObjective() = new("de_jong_maximin", "max", 1, MaximinObjectiveCore[])
end
