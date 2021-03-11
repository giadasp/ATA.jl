mutable struct NoObjModel <: AbstractModel
    settings::Settings
    constraints::Vector{Constraint}
    obj::NoObjective
    output::Output

    NoObjModel(settings, constraints, obj, output) = new(settings, constraints, obj, output)
    NoObjModel() = new(Settings(), Constraint[], NoObjective(), Output())
end
