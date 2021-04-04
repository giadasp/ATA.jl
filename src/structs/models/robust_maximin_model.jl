mutable struct RobustMaximinModel <: AbstractModel
    settings::Settings
    constraints::Vector{Constraint}
    obj::SoysterMaximinObjective
    output::Output
    RobustMaximinModel(settings, constraints, obj, output) =
        new(settings, constraints, obj, output)
        RobustMaximinModel() =
        new(Settings(), Constraint[], RobustMaximinObjective(), Output())
end
