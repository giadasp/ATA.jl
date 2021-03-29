mutable struct MaximinObjective <: AbstractObjective
    name::String
    sense::String
    cores::Vector{MaximinObjectiveCore}

    MaximinObjective(cores::Vector{MaximinObjectiveCore}) = new("maximin", "max", cores)
    MaximinObjective() = new("maximin", "max", MaximinObjectiveCore[])
end
