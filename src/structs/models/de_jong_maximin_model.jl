mutable struct DeJongMaximinModel <: AbstractModel
    settings::Settings
    constraints::Vector{Constraint}
    obj::DeJongMaximinObjective
    output::Output
    DeJongMaximinModel(settings, constraints, obj, output) =
        new(settings, constraints, obj, output)
    DeJongMaximinModel() = new(Settings(), Constraint[], DeJongMaximinObjective(), Output())
end
