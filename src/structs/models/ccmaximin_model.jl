mutable struct CcMaximinModel <: AbstractModel
    settings::Settings
    constraints::Vector{Constraint}
    obj::CcMaximinObjective
    output::Output
    CcMaximinModel(settings, constraints, obj, output) =
        new(settings, constraints, obj, output)
    CcMaximinModel() = new(Settings(), Constraint[], CcMaximinObjective(), Output())
end
