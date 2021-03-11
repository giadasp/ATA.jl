mutable struct Neighbourhood
    x::Matrix{Float64}
    f::Float64
    obj::Vector{Float64}
    infeas::Vector{Float64}
    ol::Vector{Float64}
    iu::Float64
    Neighbourhood(x, f, obj, infeas, ol, iu) = new(x, f, obj, infeas, ol, iu)
    Neighbourhood() = new(
        Matrix{Float64}(undef, 0, 0),
        Inf,
        Float64[],
        Float64[],
        Float64[],
        zero(Float64),
    )
end
