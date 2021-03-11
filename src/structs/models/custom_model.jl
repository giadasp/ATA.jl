mutable struct CustomModel <: AbstractModel
    settings::Settings
    constraints::Vector{Constraint}
    obj::CustomObjective
    output::Output
    CustomModel(settings, constraints, obj, output) =
        new(settings, constraints, obj, output)
    CustomModel() = new(Settings(), Constraint[], CustomObjective(), Output())
end
