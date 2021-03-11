mutable struct Output
    categories::Vector{Symbol}
    quantitative::Vector{Symbol}
    summ_quan::Vector{Function}
    design::Matrix{Float64}
    f::Float64
    feas::Vector{Float64}
    TIF::Vector{Float64}
    elapsed_time::Float64
    neighbourhoods::Vector{Neighbourhood}
    infos::Vector{Vector{String}}
    Output(
        categories,
        quantitative,
        summ_quan,
        design,
        f,
        feas,
        TIF,
        elapsed_time,
        neighbourhoods,
        infos,
    ) = new(
        categories,
        quantitative,
        summ_quan,
        design,
        f,
        feas,
        TIF,
        elapsed_time,
        neighbourhoods,
        infos,
    )
    Output() = new(
        Vector{Vector{Symbol}}(undef, 0),
        Vector{Vector{Symbol}}(undef, 0),
        Function[],
        zeros(Float64, 0, 0),
        zero(Float64),
        Float64[],
        Float64[],
        zero(Float64),
        Neighbourhood[],
        Vector{Vector{String}}(undef, 0),
    )
end
