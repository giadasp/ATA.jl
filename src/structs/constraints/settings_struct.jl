include("friend_sets_struct.jl")
include("enemy_sets_struct.jl")
include("item_use_struct.jl")
include("irt_struct.jl")

mutable struct Settings
    n_items::Int64
    n_fs::Int64
    bank::DataFrames.DataFrame
    irt::IRT
    theta_bounds::Vector{Vector{Float64}}
    forced0::Vector{Vector{Bool}}
    n_groups::Int64
    T::Int64
    Tg::Vector{Int64}
    fs::FriendSets # friend Sets
    es::EnemySets # enemy Sets
    iu::ItemUse
    ol_max::Matrix{Float64}
    to_apply::Vector{Bool} #[iu_max?, iu_min?, ol?]
    Settings(
        n_items,
        n_fs,
        bank,
        irt,
        theta_bounds,
        forced0,
        n_groups,
        T,
        Tg,
        fs,
        es,
        iu,
        ol_max,
        to_apply,
    ) = new(
        n_items,
        n_fs,
        bank,
        irt,
        theta_bounds,
        forced0,
        n_groups,
        T,
        Tg,
        fs,
        es,
        iu,
        ol_max,
        to_apply,
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
        FriendSets(),
        EnemySets(),
        ItemUse(),
        zeros(Int64, 0, 0),
        [false, false, false],
    )
end
