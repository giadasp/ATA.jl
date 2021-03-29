mutable struct MaximinObjectiveCore
    points::Vector{Float64}
    IIF::Matrix{Float64}
    MaximinObjectiveCore() = new(Vector{Float64}(), Matrix{Float64}(undef, 0, 0))
end

mutable struct MaximinObjective <: AbstractObjective
    name::String
    sense::String
    cores::Vector{MaximinObjectiveCore}

    MaximinObjective(cores::Vector{MaximinObjectiveCore}) = new("maximin", "max", cores)
    MaximinObjective() = new("maximin", "max", MaximinObjectiveCore[])
end
