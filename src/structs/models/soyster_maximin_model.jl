mutable struct SoysterMaximinModel <: AbstractModel
    settings::Settings
    constraints::Vector{Constraint}
    obj::MaximinObjective
    output::Output
    SoysterMaximinModel(settings, constraints, obj, output) =
        new(settings, constraints, obj, output)
    SoysterMaximinModel() = new(Settings(), Constraint[], MaximinObjective(), Output())
end
