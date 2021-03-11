Inputs = InputSettings(

    # 1. T
    [8, 8, 8],
    # 2. n_items
    Int(366),
    # 3. n_groups
    Int(3),
    # 4. groups
    ["D", "C", "B"],

    ####################################################################
    ########################## IRT #####################################
    ####################################################################

    # 5. IRT_model
    "1PL",
    # 6. IRT_parameters
    ["b"],
    # 7. IRT_parametrization
    "at-ab",
    # 8. IRT_D
    1,

    ####################################################################
    ######################FRIENDS AND ENEMIes###########################
    ####################################################################

    # 9. enemy_sets_var
    ["ENEMY_SET"],
    # 10. friend_sets_var
    ["UNIT"],

    ####################################################################
    ########################### ITEM USE  ##############################
    ####################################################################

    # 11. item_use_min
    fill(0, 366),
    # 12. item_use_max
    fill(7, 366),

    ####################################################################
    ########################## TesT LENGHT #############################
    ####################################################################

    # 13. lenght_min
    [36, 36, 36],
    # 14. lenght_max
    [40, 40, 40],
    # 15. lenght_weight
    [1.0, 1.0, 1.0],

    ####################################################################
    ####################### EXPECTED SCORE #############################
    ####################################################################

    # 16. expected_score_var
    ["PROP_CORR", "PROP_CORR", "PROP_CORR"],
    # 17. expected_score_pts
    [zeros(Float64, 1), zeros(Float64, 1), zeros(Float64, 1)],
    # 18. expected_score_min
    [[0.50], [0.50], [0.50]],
    # 19. expected_score_max
    [[0.57], [0.57], [0.57]],

    ####################################################################
    ###################### GENERIC CONSTRAINTS#################
    ####################################################################

    # 20. mean_vars
    Vector{Vector{String}}(undef, 0),  # (future)
    # 21. mean_vars_min
    Vector{Vector{Float64}}(undef, 0), # (future)
    # 22. mean_vars_max
    Vector{Vector{Float64}}(undef, 0), # (future)
    # 23. sum_vars
    Vector{Vector{String}}(undef, 0),
    # 24. sum_vars_min
    Vector{Vector{Float64}}(undef, 0),
    # 25. sum_vars_max
    Vector{Vector{Float64}}(undef, 0),

    ####################################################################
    ######################### OBJECTIVE ################################
    ####################################################################

    # 26. obj_type
    "custom", # "MAXIMIN", "CCMAXIMIN", "", "custom", "MINIMAX"
    # 27. obj_points (required in MAXIMIN, CC, MINIMAX)
    [[-0.60], [0.30], [0.60]],
    # 28. obj_targets (required in MINIMAX)
    Vector{Vector{Float64}}(),
    # 29. obj_aux_int
    zero(Float64),
    # 30. obj_aux_float
    0.05,

    ####################################################################
    ######################### OUTPUT ###################################
    ####################################################################
    # 31. categories
    ["UNIT", "CAT_1", "CAT_2", "CAT_3", "CAT_4", "CAT_5_6", "CAT_5", "CAT_6"],
);
