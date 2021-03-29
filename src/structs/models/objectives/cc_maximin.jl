mutable struct CcMaximinObjectiveCore
    R::Int64
    alpha::Float64
    points::Vector{Float64}
    IIF::Array{Float64,3}
    CcMaximinObjectiveCore() =
        new(1, 0.05, Vector{Float64}(), Array{Float64,3}(undef, 0, 0, 0))
end

mutable struct CcMaximinObjective <: AbstractObjective
    name::String
    sense::String
    cores::Vector{CcMaximinObjectiveCore}
    CcMaximinObjective(cores) = new("cc_maximin", "max", cores)
    CcMaximinObjective() = new("cc_maximin", "max", CcMaximinObjectiveCore[])
end
