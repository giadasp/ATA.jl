mutable struct SoysterMaximinObjective <: AbstractObjective
    name::String
    sense::String
    cores::Vector{MaximinObjectiveCore}

    SoysterMaximinObjective(cores::Vector{MaximinObjectiveCore}) = new("soyster_maximin", "max", cores)
    SoysterMaximinObjective() = new("soyster_maximin", "max", MaximinObjectiveCore[])
end
