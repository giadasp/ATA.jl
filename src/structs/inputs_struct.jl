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
    
    function InputSettings(;
        T = Int64[],
        n_items = zero(Int64),
        n_groups = zero(Int64),
        groups = String[],
        IRT_model = "",
        IRT_parameters = String[],
        IRT_parametrization = "at-b",
        IRT_D = 1.0,
        enemy_sets_var = String[],
        friend_sets_var = String[],
        item_use_min = Int64[],
        item_use_max = Int64[],
        length_min = Int64[],
        length_max = Int64[],
        length_weight = Float64[],
        expected_score_var = String[],
        expected_score_pts = Float64[],
        expected_score_min = Float64[],
        expected_score_max = Float64[],
        mean_vars = Vector{Vector{String}}(undef, 0),
        mean_vars_min = Vector{Vector{Float64}}(undef, 0),
        mean_vars_max = Vector{Vector{Float64}}(undef, 0),
        sum_vars = Vector{Vector{String}}(undef, 0),
        sum_vars_min = Vector{Vector{Float64}}(undef, 0),
        sum_vars_max = Vector{Vector{Float64}}(undef, 0),
        obj_type = "",
        obj_points = Vector{Vector{Float64}}(undef, 0),
        obj_targets = Vector{Vector{Float64}}(undef, 0),
        obj_aux_int = zero(Int64),
        obj_aux_float = zero(Float64),
        categories = String[],
    ) 
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
    )
    end
end
