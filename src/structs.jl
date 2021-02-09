mutable struct InputSettings
    T::Vector{Int64}
    n_items::Int64
    n_groups::Int64
    groups::Vector{String}
    IRT_model::String
    IRT_parameters::Vector{String}
    IRT_parametrization::String
    IRT_D::Float64
    enemy_sets_var::Vector{String}
    friend_sets_var::Vector{String}
    item_use_min::Vector{Int64}
    item_use_max::Vector{Int64}
    length_min::Vector{Int64}
    length_max::Vector{Int64}
    length_weight::Vector{Float64}
    expected_score_var::Vector{String}
    expected_score_pts::Vector{Vector{Float64}}
    expected_score_min::Vector{Vector{Float64}}
    expected_score_max::Vector{Vector{Float64}}
    mean_vars::Vector{Vector{String}}
    mean_vars_min::Vector{Vector{Float64}}
    mean_vars_max::Vector{Vector{Float64}}
    sum_vars::Vector{Vector{String}}
    sum_vars_min::Vector{Vector{Float64}}
    sum_vars_max::Vector{Vector{Float64}}
    obj_type::String
    obj_points::Vector{Vector{Float64}}
    obj_targets::Vector{Vector{Float64}}
    obj_aux_int::Int64
    obj_aux_float::Float64
    categories::Vector{String}
    InputSettings(
        T,
        n_items,
        n_groups,
        groups,
        IRT_model,
        IRT_parameters,
        IRT_parametrization,
        IRT_D,
        enemy_sets_var,
        friend_sets_var,
        item_use_min,
        item_use_max,
        length_min,
        length_max,
        length_weight,
        expected_score_var,
        expected_score_pts,
        expected_score_min,
        expected_score_max,
        mean_vars,
        mean_vars_min,
        mean_vars_max,
        sum_vars,
        sum_vars_min,
        sum_vars_max,
        obj_type,
        obj_points,
        obj_targets,
        obj_aux_int,
        obj_aux_float,
        categories,
    ) = new(
        T,
        n_items,
        n_groups,
        groups,
        IRT_model,
        IRT_parameters,
        IRT_parametrization,
        IRT_D,
        enemy_sets_var,
        friend_sets_var,
        item_use_min,
        item_use_max,
        length_min,
        length_max,
        length_weight,
        expected_score_var,
        expected_score_pts,
        expected_score_min,
        expected_score_max,
        mean_vars,
        mean_vars_min,
        mean_vars_max,
        sum_vars,
        sum_vars_min,
        sum_vars_max,
        obj_type,
        obj_points,
        obj_targets,
        obj_aux_int,
        obj_aux_float,
        categories,
    )
    InputSettings() = new(
        Int64[],
        zero(Int64),
        zero(Int64),
        String[],
        "",
        String[],
        "at-b",
        1.0,
        String[],
        String[],
        Int64[],
        Int64[],
        Int64[],
        Int64[],
        Float64[],
        String[],
        Float64[],
        Float64[],
        Float64[],
        Vector{Vector{String}}(undef, 0),
        Vector{Vector{Float64}}(undef, 0),
        Vector{Vector{Float64}}(undef, 0),
        Vector{Vector{String}}(undef, 0),
        Vector{Vector{Float64}}(undef, 0),
        Vector{Vector{Float64}}(undef, 0),
        Vector{Vector{Float64}}(undef, 0),
        "",
        Vector{Vector{Float64}}(undef, 0),
        zero(Int64),
        zero(Float64),
        String[],
    )
end

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

mutable struct FS
    var::Vector{Symbol}
    sets::Vector{String}
    counts::Vector{Int64}
    items::Vector{Vector{Int64}}
    FS(var, sets, counts, items) = new(var, sets, counts, items)
    FS() = new(Symbol[], String[], Int64[], Vector{Vector{Int64}}(undef, 0))
end

mutable struct ES
    var::Vector{Symbol}
    names::Vector{String}
    sets::Vector{Vector{Int64}}
   ES(var, names, sets) = new(var, names, sets)
   ES() = new(Symbol[], String[], Vector{Vector{Int64}}(undef, 0))
end

mutable struct ExpectedScore
    var::Symbol
    val::Matrix{Float64}
    min::Vector{Float64}
    max::Vector{Float64}
    pts::Vector{Float64}
    ExpectedScore(var, val, min, max, pts) = new(var, val, min, max, pts)
    ExpectedScore() = new(Symbol(""), zeros(Float64, 0, 0), Float64[], Float64[], Float64[])
end

mutable struct IU
    min::Vector{Int64}
    max::Vector{Int64}
    IU(min, max) = new(min, max)
    IU() = new(Int64[], Int64[])
end

mutable struct Settings
    n_items::Int64
    n_fs::Int64
    bank::DataFrames.DataFrame
    IRT::IRT
    theta_bounds::Vector{Vector{Float64}}
    forced0::Vector{Vector{Bool}}
    n_groups::Int64
    T::Int64
    Tg::Vector{Int64}
    fs::FS # friend Sets
    es::ES # enemy Sets
    iu::IU
    ol_max::Matrix{Float64}
    Settings(
        n_items,
        n_fs,
        bank,
        IRT,
        theta_bounds,
        forced0,
        n_groups,
        T,
        Tg,
        fs,
        es,
        iu,
        ol_max,
    ) = new(
        n_items,
        n_fs,
        bank,
        IRT,
        theta_bounds,
        forced0,
        n_groups,
        T,
        Tg,
        fs,
        es,
        iu,
        ol_max,
    ) # no pattern mode
    Settings() = new(
        0,
        0,
        DataFrames.DataFrame(),
        IRT(),
        [[-6.0, 6.0]],
        Vector{Vector{Bool}}(undef, 0),
        1,
        1,
        [1],
        FS(),
        ES(),
        IU(),
        zeros(Int64, 0, 0),
    )
end

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

mutable struct Obj
    sense::String
    type::String
    points::Vector{Vector{Float64}}
    targets::Vector{Vector{Float64}}
    aux_int::Int64
    aux_float::Float64
    # must be test separable and accept x_Iₜ (design of test v, item level, not grouped by friend sets) as first argument and NamedTuple as second argument, it must return a Vector of Float64 with T elements.
    fun::Function
    args::NamedTuple # ex: (a = [0, 0, 0], b = "hello") 
    Obj(sense, type, points, targets, aux_int, aux_float, fun, args) =
        new(sense, type, points, targets, aux_int, aux_float, fun, args)
    Obj() = new(
        "max",
        "",
        Vector{Vector{Float64}}(undef, 0),
        Vector{Vector{Float64}}(undef, 0),
        zero(Int64),
        zero(Float64),
        () -> nothing,
        (default_arg_1 = 0, default_arg_2 = 0),
    )
end

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

mutable struct Model
    settings::Settings
    constraints::Vector{Constraint}
    obj::Obj
    output::Output
    Model(settings, constraints, obj, output) =
        new(settings, constraints, obj, output)
    Model() = new(Settings(), Constraint[], Obj(), Output())
end

mutable struct Opt
    x::Matrix{Float64}
    f::Float64
    obj::Vector{Float64}
    infeas::Vector{Float64}
    ol::Vector{Float64}
    iu::Vector{Float64}
    Opt(x, f, obj, infeas, ol, iu) = new(x, f, obj, infeas, ol, iu)
    Opt() =
        new(Matrix{Float64}(undef, 0, 0), Inf, Float64[], Float64[], Float64[], Float64[])
end
