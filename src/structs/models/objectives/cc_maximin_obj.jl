mutable struct CCMaximinObjective <: AbstractObjective
    name::String
    sense::String
    R::Int64
    cores::Vector{CCMaximinObjectiveCore}

    CCMaximinObjective(R::Int64, cores::Vector{CCMaximinObjectiveCore}) = new("cc_maximin", "max", R, cores)
    CCMaximinObjective(cores::Vector{CCMaximinObjectiveCore}) = new("cc_maximin", "max", 1, cores)
    CCMaximinObjective() = new("cc_maximin", "max", 1, CCMaximinObjectiveCore[])
end
