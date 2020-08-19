include("simAn.jl")
include("jump.jl")
function assemble!(ATAmodel::Model; solver="jumpATA", starting_design = Matrix{Float64}(undef,0,0), start_temp = 0.1, geom_temp = 0.1, max_time = 1000.00, results_folder = "RESULTS", n_item_sample = 1, n_test_sample = 1, conv_max = 2, verbosity = 2,  opt_feas = 0.0, n_fill = 1 , feas_nh = 0, opt_nh = 5, optimizer_constructor = "GLPK", optimizer_attributes = [("tm_lim",1000)])
    if solver=="siman"
        siman!(ATAmodel; starting_design = starting_design, start_temp = start_temp, geom_temp = geom_temp, max_time = max_time, results_folder = results_folder, n_item_sample = n_item_sample, n_test_sample = n_test_sample, conv_max = conv_max, verbosity = verbosity, opt_feas = opt_feas , n_fill = n_fill , feas_nh = feas_nh, opt_nh = opt_nh)
    elseif solver=="jumpATA"
        jumpATA!(ATAmodel; starting_design = starting_design, optimizer_constructor = optimizer_constructor, optimizer_attributes = optimizer_attributes, results_folder = results_folder)
    end
end
