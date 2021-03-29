mutable struct DeJongMaximinObjective <: AbstractObjective
    name::String
    sense::String
    cores::Vector{MaximinObjectiveCore}

    DeJongMaximinObjective(cores::Vector{MaximinObjectiveCore}) = new("de_jong_maximin", "max", cores)
    DeJongMaximinObjective() = new("de_jong_maximin", "max", MaximinObjectiveCore[])
end
