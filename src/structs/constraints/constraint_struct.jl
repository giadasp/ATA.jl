include("expected_score_struct.jl")
mutable struct Constraint
    length_min::Int64
    length_max::Int64
    expected_score::ExpectedScore
    mean_vars::Vector{Symbol} # constrain the mean to be
    mean_vars_min::Vector{Float64}
    mean_vars_max::Vector{Float64}
    sum_vars::Vector{Symbol}# constrain the sum to be
    sum_vars_min::Vector{Float64}
    sum_vars_max::Vector{Float64}
    constr_A::Matrix{Float64}
    constr_b::Vector{Float64}
    # ol_max::Matrix{Int64}
    Constraint(
        length_min,
        length_max,
        expected_score,
        mean_vars,
        mean_vars_min,
        mean_vars_max,
        sum_vars,
        sum_vars_min,
        sum_vars_max,
        constr_A,
        constr_b,
    ) = new(
        length_min,
        length_max,
        expected_score,
        mean_vars,
        mean_vars_min,
        mean_vars_max,
        sum_vars,
        sum_vars_min,
        sum_vars_max,
        constr_A,
        constr_b,
    )
    Constraint() = new(
        zero(Int64),
        10000000,
        ExpectedScore(),
        Vector{Symbol}(undef, 0),
        Vector{Float64}(undef, 0),
        Vector{Float64}(undef, 0),
        Vector{Symbol}(undef, 0),
        Vector{Float64}(undef, 0),
        Vector{Float64}(undef, 0),
        Matrix{Float64}(undef, 0, 0),
        Vector{Float64}(undef, 0),
    )
end
