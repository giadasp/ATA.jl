mutable struct CCMaximinModel <: AbstractModel
    settings::Settings
    constraints::Vector{Constraint}
    obj::CCMaximinObjective
    output::Output
    CCMaximinModel(settings, constraints, obj, output) =
        new(settings, constraints, obj, output)
    CCMaximinModel() = new(Settings(), Constraint[], CCMaximinObjective(), Output())
end
