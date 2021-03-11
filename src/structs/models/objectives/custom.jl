mutable struct CustomObjective <: AbstractObjective
    name::String
    sense::String
    fun::Function
    args::NamedTuple # ex: (a = [0, 0, 0], b = "hello") 
    CustomObjective(name, sense, fun, args) = new(name, sense, fun, args)
    CustomObjective() =
        new("nothing fun", "max", () -> nothing, (default_arg_1 = 0, default_arg_2 = 0))
end
