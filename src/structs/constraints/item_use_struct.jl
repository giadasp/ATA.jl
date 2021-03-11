mutable struct ItemUse
    min::Vector{Int64}
    max::Vector{Int64}
    ItemUse(min, max) = new(min, max)
    ItemUse() = new(Int64[], Int64[])
end
