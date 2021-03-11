struct NoObjective <: AbstractObjective
    name::String
    NoObjective() = new("feasibility model")
end
