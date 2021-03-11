abstract type AbstractModel end

include("objectives/obj.jl")
include("maximin_model.jl")
include("minimax_model.jl")
include("custom_model.jl")
include("ccmaximin_model.jl")
include("no_obj_model.jl")
