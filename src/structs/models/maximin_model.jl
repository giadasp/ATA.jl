mutable struct MaximinModel <: AbstractModel
    settings::Settings
    constraints::Vector{Constraint}
    obj::MaximinObjective
    output::Output
    MaximinModel(settings, constraints, obj, output) =
        new(settings, constraints, obj, output)
    MaximinModel() = new(Settings(), Constraint[], MaximinObjective(), Output())
end
