mutable struct CCMaximinModel <: AbstractModel
    settings::Settings
    constraints::Vector{Constraint}
    obj::CcMaximinObjective
    output::Output
    CCMaximinModel(settings, constraints, obj, output) =
        new(settings, constraints, obj, output)
    CCMaximinModel() = new(Settings(), Constraint[], CcMaximinObjective(), Output())
end
