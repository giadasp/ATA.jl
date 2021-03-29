mutable struct SoysterMaximinObjective <: AbstractObjective
    name::String
    sense::String
    R::Int64
    cores::Vector{MaximinObjectiveCore}

    SoysterMaximinObjective(R::Int64, cores::Vector{MaximinObjectiveCore}) =
        new("soyster_maximin", "max", R, cores)
    SoysterMaximinObjective(cores::Vector{MaximinObjectiveCore}) =
        new("soyster_maximin", "max", 1, cores)
    SoysterMaximinObjective() = new("soyster_maximin", "max", 1, MaximinObjectiveCore[])
end
