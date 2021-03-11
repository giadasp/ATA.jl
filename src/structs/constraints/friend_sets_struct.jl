mutable struct FriendSets
    var::Vector{Symbol}
    sets::Vector{String}
    counts::Vector{Int64}
    items::Vector{Vector{Int64}}
    FriendSets(var, sets, counts, items) = new(var, sets, counts, items)
    FriendSets() = new(Symbol[], String[], Int64[], Vector{Vector{Int64}}(undef, 0))
end
