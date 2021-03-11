mutable struct ExpectedScore
    var::Symbol
    val::Matrix{Float64}
    min::Vector{Float64}
    max::Vector{Float64}
    pts::Vector{Float64}
    ExpectedScore(var, val, min, max, pts) = new(var, val, min, max, pts)
    ExpectedScore() = new(Symbol(""), zeros(Float64, 0, 0), Float64[], Float64[], Float64[])
end
