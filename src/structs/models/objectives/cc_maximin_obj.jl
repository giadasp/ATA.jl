mutable struct CCMaximinObjective <: AbstractObjective
    name::String
    sense::String
    cores::Vector{CCMaximinObjectiveCore}
    CCMaximinObjective(cores) = new("cc_maximin", "max", cores)
    CCMaximinObjective() = new("cc_maximin", "max", CCMaximinObjectiveCore[])
end
