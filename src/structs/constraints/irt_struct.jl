mutable struct IRT
    model::String
    parameters::DataFrames.DataFrame
    parametrization::String
    D::Float64
    metric::Vector{Float64} # [meanTarget,stdTarget]
    X::Vector{Float64}
    W::Vector{Float64}
    IRT() = new("2PL", DataFrames.DataFrame(), "at-b", 1.0, [0.0, 1.0], zeros(1), zeros(1))
    IRT(model, parameters, parametrization, D, metric, X, W) =
        new(model, parameters, parametrization, D, metric, X, W)
end
