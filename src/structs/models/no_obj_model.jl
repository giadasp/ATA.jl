mutable struct NoObjModel <: AbstractModel
    settings::Settings
    constraints::Vector{Constraint}
    output::Output
    NoObjModel(settings, constraints, output) = new(settings, constraints, output)
    NoObjModel() = new(Settings(), Constraint[], Output())
end
