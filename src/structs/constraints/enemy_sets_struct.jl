mutable struct EnemySets
    var::Vector{Symbol}
    names::Vector{String}
    sets::Vector{Vector{Int64}}
    EnemySets(var, names, sets) = new(var, names, sets)
    EnemySets() = new(Symbol[], String[], Vector{Vector{Int64}}(undef, 0))
end
