mutable struct DeJongMaximinModel <: AbstractModel
    settings::Settings
    constraints::Vector{Constraint}
    obj::MaximinObjective
    output::Output
    DeJongMaximinModel(settings, constraints, obj, output) =
        new(settings, constraints, obj, output)
    DeJongMaximinModel() = new(Settings(), Constraint[], MaximinObjective(), Output())
end
