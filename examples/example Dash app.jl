# cd("folder in which the package is saved")
# using Pkg
# Pkg.activate(".")  # required
# Pkg.instantiate()
# cd("where your input files are")
using ATA
using Dash
using DashCoreComponents
using DashHtmlComponents

# Before running the app, if you want to use a MILP solver, remember to load it
# (ex: using Cbc; run_app!()).
run_app!()
# Navigate with the browser to localhost:8080
