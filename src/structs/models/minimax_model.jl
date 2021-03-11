mutable struct MinimaxModel <: AbstractModel
    settings::Settings
    constraints::Vector{Constraint}
    obj::MinimaxObjective
    output::Output
    MinimaxModel(settings, constraints, obj, output) =
        new(settings, constraints, obj, output)
    MinimaxModel() = new(Settings(), Constraint[], MinimaxObjective(), Output())
end
