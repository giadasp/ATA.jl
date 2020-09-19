function run_app!()
  global ATAmodel = Model()
  js_files = [
  Dict(
  "src"=> "https://code.jquery.com/jquery-3.5.1.slim.min.js",
  "integrity" => "sha384-DfXdz2htPH0lsSSs5nCTpuj/zy4C+OGpamoFVy38MVBnE+IbbVYUew+OrCXaRkfj",
  "crossorigin" => "anonymous",
  ),
  Dict(
  "src"=> "https://stackpath.bootstrapcdn.com/bootstrap/4.5.0/js/bootstrap.min.js",
  "integrity" => "sha384-OgVRvuATP1z7JjHLkuOU7Xw704+h835Lr+6QL9UvYjZE3Ipu6Tp75j7Bh/kR0JKI",
  "crossorigin" => "anonymous"
  )
  ]
  css_files = [
  Dict(
  "href"=> "https://stackpath.bootstrapcdn.com/bootstrap/4.5.0/css/bootstrap.min.css",
  "integrity" => "sha384-9aIt2nRpC12Uk9gS9baDl411NQApFmC26EwAOH8WgZl5MYYxFfc+NcPb1dKGj7Sk",
  "crossorigin" => "anonymous",
  "rel" => "stylesheet"
  )
  ,
  Dict(
  "href"=> "https://codepen.io/giadasp/pen/MWKbOJv.css",
  "rel" => "stylesheet"
  )
  ]

  app = Dash.dash(external_scripts=js_files, external_stylesheets=css_files, assets_folder="css")
  app.layout =
  DashCoreComponents.dcc_tabs([
  DashCoreComponents.dcc_tab(label = "Create model", children = [
  DashHtmlComponents.html_div( [
  DashHtmlComponents.html_div( className="modal-dialog", role="document",[
  DashHtmlComponents.html_div( className="modal-content", [
  DashHtmlComponents.html_div( className="modal-body", [
  DashHtmlComponents.html_div( className="spinner-border", role="status", [
  DashHtmlComponents.html_span("Loading...", className="sr-only")
  ])
  ])
  ])
  ])
  ];
  id="loader", className="modal fade", tabIndex="-1", role="dialog", (Symbol("aria-labelledby")=>"exampleModalLabel") ,(Symbol("aria-hidden")=>"true")),
  DashHtmlComponents.html_div(className = "container-fluid", [
  DashHtmlComponents.html_div(className = "row-md-12 justify-content-center m-3", [
  DashHtmlComponents.html_div(className = "col-md-12 ", [

  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6 ", [
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_label("folder:", className="mr-1"),
  DashCoreComponents.dcc_input(id="folder_txt", value=pwd(), type = "text", className="form-control")
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1 text-center", [
  DashHtmlComponents.html_button("Start ATA"; id = "start_btn", n_clicks=0, className = "btn btn-primary")
  ])
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_div(id = "restart_card", className = "card bg-light h-100", [
  DashHtmlComponents.html_div(className = "card-body", [
  DashHtmlComponents.html_p(children=[""], id = "start_lbl"),
  ])

  ])
  ])
  ])
  ])
  ]),


  DashHtmlComponents.html_div(className = "row d-md-flex justify-content-md-center m-3", [
  #DashHtmlComponents.html_p("hidden",id ="settings_loader", style="display: none;"),
  DashHtmlComponents.html_div(className = "col-md-12", [

  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("file name:"),
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(id="settings_txt", value="SettingsATA.jl", className="form-control", type = "text")
  ])
  ]),
  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("item pool file:", id="pool_lbl"),
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(id="bank_txt", value="bank.csv", className="form-control", type = "text"),
  ])
  ]),
  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("bank delimitator:", id="delim_bank_lbl"),
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(id="bank_delim_txt", value=";", className="form-control", type = "text"),
  ])
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center text-center", [
  DashHtmlComponents.html_button("Load Settings"; id="settings_btn", n_clicks=0, className = "btn btn-primary")
  ]),
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_div(id = "settings_card", className = "card h-100 bg-light", [
  DashHtmlComponents.html_div(className = "card-body", [
  DashHtmlComponents.html_p(children=[""], id = "settings_lbl")
  ])

  ])
  ])
  ])
  ])
  ]),


  DashHtmlComponents.html_div(className = "row justify-content-center m-3", [
  DashHtmlComponents.html_div(className = "col-md-12", [
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6 align-self-center text-center", [
  DashHtmlComponents.html_button("Add Friend sets"; id="addFS_btn", n_clicks=0, className = "btn btn-primary")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_div(id = "addFS_card", className = "card h-100 bg-light", [
  DashHtmlComponents.html_div(className = "card-body", [
  DashHtmlComponents.html_p(children=[""], id = "addFS_lbl")
  ])
  ])

  ])
  ])
  ])
  ]),


  DashHtmlComponents.html_div(className = "row justify-content-center m-3", [
  DashHtmlComponents.html_div(className = "col-md-12", [

  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6 align-self-center text-center", [
  DashHtmlComponents.html_button("Add Enemy sets"; id="addES_btn", n_clicks=0, className = "btn btn-primary")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_div(id = "addES_card", className = "card h-100 bg-light", [
  DashHtmlComponents.html_div(className = "card-body", [
  DashHtmlComponents.html_p(children=[""], id = "addES_lbl")
  ])

  ])
  ])
  ])
  ])
  ]),


  DashHtmlComponents.html_div(className = "row justify-content-center m-3", [
  DashHtmlComponents.html_div(className = "col-md-12", [

  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("file name:"),
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(id="constraints_file_txt", value="Constraints.csv", className="form-control", type = "text")
  ])
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("delimitator:", id="constraints_delim_lbl"),
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(id="constraints_delim_txt", value=";", className="form-control", type = "text"),
  ])
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1 text-center", [
  DashHtmlComponents.html_button("Add constraints"; id="addConstraints_btn", n_clicks=0, className = "btn btn-primary")
  ]),
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_div(id = "addConstraints_card", className = "card h-100 bg-light", [
  DashHtmlComponents.html_div(className = "card-body", [
  DashHtmlComponents.html_p(children=[""], id = "addConstraints_lbl"),
  ])
  ])

  ])
  ])
  ])
  ]),


  DashHtmlComponents.html_div(className = "row justify-content-center m-3", [
  DashHtmlComponents.html_div(className = "col-md-12", [

  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("file name:"),
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(id="ol_txt", value="OverlapMatrix.csv", className="form-control", type = "text"),
  ])
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("matrix delimitator:", id="ol_delim_lbl"),
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(id="ol_delim_txt", value=";", className="form-control", type = "text"),
  ])
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center text-center m-1", [
  DashHtmlComponents.html_button("Add Overlap Matrix"; id="addOl_btn", n_clicks=0, className = "btn btn-primary")
  ]),
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_div(id = "ol_card",className = "card h-100 bg-light", [
  DashHtmlComponents.html_div(className = "card-body", [
  DashHtmlComponents.html_p(children=[""], id = "ol_lbl"),
  ])

  ])
  ])
  ])
  ])
  ]),

  DashHtmlComponents.html_div(className = "row justify-content-center m-3", [
  DashHtmlComponents.html_div(className = "col-md-12 align-self-center", [

  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6 align-self-center text-center", [
  DashHtmlComponents.html_button("Add Expected Score"; id="addExpScore_btn", n_clicks=0, className = "btn btn-primary")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6 align-self-center", [
  DashHtmlComponents.html_div(id = "addExpScore_card", className = "card h-100 bg-light", [
  DashHtmlComponents.html_div(className = "card-body", [
  DashHtmlComponents.html_p(children=[""], id = "addExpScore_lbl"),
  ])
  ])
  ])

  ])
  ])
  ]),


  DashHtmlComponents.html_div(className = "row justify-content-center m-3", [
  DashHtmlComponents.html_div(className = "col-md-12 align-self-center", [

  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6 align-self-center text-center", [
  DashHtmlComponents.html_button("Group By Friend sets"; id="group_btn", n_clicks=0, className = "btn btn-primary")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6 align-self-center", [
  DashHtmlComponents.html_div(id = "group_card", className = "card h-100 bg-light", [
  DashHtmlComponents.html_div(className = "card-body", [
  DashHtmlComponents.html_p(children=[""], id = "group_lbl"),
  ])
  ])
  ])

  ])
  ])
  ]),


  DashHtmlComponents.html_div(className = "row justify-content-center m-3", [
  DashHtmlComponents.html_div(className = "col-md-12 align-self-center ", [

  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6 align-self-center text-center", [
  DashHtmlComponents.html_button("Add objective function"; id="addObj_btn", n_clicks=0, className = "btn btn-primary")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6 align-self-center", [
  DashHtmlComponents.html_div(id = "addObj_card", className = "card h-100 bg-light", [
  DashHtmlComponents.html_div(className = "card-body", [
  DashHtmlComponents.html_p(children=[""], id = "addObj_lbl"),
  ])
  ])
  ])

  ])
  ])
  ])

  ])
  ]),
  DashCoreComponents.dcc_tab(label = "assemble!", children = [
  DashHtmlComponents.html_div(className = "container-fluid", [
  DashHtmlComponents.html_div(className = "row m-3", [
  DashHtmlComponents.html_div(className = "col-md-12 align-self-center", [
  DashHtmlComponents.html_div(className = "card h-100", [
  DashHtmlComponents.html_div(className = "card-body", [
  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_h3("SimAn solver", className="text-center"),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Maximum time (in seconds):")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value=100,className="form-control", type="number", min=0, max = 10000, step=10, id="max_time_txt")
  ]),
  ]),

  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Starting temperature:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value=0.1,className="form-control", type="number", min=0, max = 1000, step=0.00001, id="start_temp_txt")
  ]),
  ]),
  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Geometric temp decreasing:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value=0.1, className="form-control", type="number", min=0, max = 1000, step=0.00001, id="geom_temp_txt")
  ]),
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Items to sample:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value=1,className="form-control", type="number", min=1, max = 1000, step=1, id="n_item_StatsBase.sample_txt")

  ]),
  ]),
  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Tests to sample:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value=1,className="form-control", type="number", min=1, max = 1000, step=1, id="n_test_StatsBase.sample_txt")
  ]),
  ]),
  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Maximum convergence:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value=10,className="form-control", type="number", min=2, max=1000, step=1, id="max_conv_txt")

  ]),
  ]),
  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Verbosity:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_dropdown(value=1, options=[(label="minimal",value=1),(label="detailed",value=2)], id="verbosity_txt")
  ]),
  ]),
  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Optimality/Feasibility balancer:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value=0.1, className="form-control", type="number", min=0, max=1, step=0.01, id="opt_feas_txt")
  ]),
  ]),
  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("#Fill-up phases:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value=1,className="form-control", type="number", min=0, max=2000, step=1, id="n_fill_txt")
  ]),
  ]),
  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("#Feasibility neighbourhoods:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value=1, className="form-control", type="number", min=0, max=2000, step=1, id="feas_nh_txt")
  ]),
  ]),
  DashHtmlComponents.html_div(className = "row m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("#Optimality neighbourhoods:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value=100, className="form-control", type="number", min=0, max=2000, step=1, id="opt_nh_txt")
  ]),
  ]),
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_h3("JuMP solver", className = "text-center"),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Optimizer constructor:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_dropdown(
  value="GLPK",
  id="optimizer_constructor_txt",
  options = [
  (label="GLPK", value="GLPK"),
  (label="CPLEX", value="CPLEX"),
  (label="Gurobi", value="Gurobi"),
  (label="KNITRO", value="KNITRO"),
  (label="Cbc", value="Cbc"),
  (label="Xpress", value="Xpress"),
  (label="Juniper", value="Juniper"),
  (label="CPLEX", value="CPLEX"),
  (label="MosekTools", value="MosekTools"),
  (label="SCIP", value="SCIP")
  ]
  )
  ]),
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Optimizer attributes:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value="[(\"tm_lim\",1000)]",className="form-control", type="text", id="optimizer_attributes_txt")
  ]),
  ]),

  ]),
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Starting design file:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value="design.csv",className="form-control", type="text",id="starting_design_txt")
  ]),
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Results folder name:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_input(value="RESULTS",className="form-control", type="text", id="results_folder_txt")
  ]),
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Print plots:")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_dropdown(
  value="no",
  id="plots_out_txt",
  options = [
  (label="yes", value="yes"),
  (label="no", value="no")
  ]
  )
  ]),
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashHtmlComponents.html_label("Solver (siman or jumpATA):")
  ]),
  DashHtmlComponents.html_div(className = "col-md-6", [
  DashCoreComponents.dcc_dropdown(
  value="siman",
  id="solver_txt",
  options = [
  (label="siman", value="siman"),
  (label="jumpATA", value="jumpATA")
  ]
  )
  ]),
  ]),
  DashHtmlComponents.html_div(className = "row justify-content-center m-1", [
  DashHtmlComponents.html_button("assemble!"; id="assemble_btn", n_clicks=0, className = "btn btn-success")
  ])
  ])
  ])
  ])
  ])
  ])
  ]),

  DashCoreComponents.dcc_tab(label = "Results", children = [
  DashHtmlComponents.html_div(className = "container-fluid", [
  DashHtmlComponents.html_div(className = "row-md-12 justify-content-center m-3", [
  DashHtmlComponents.html_div(className = "col-md-12", [
  DashHtmlComponents.html_div(id = "results_card", className = "card", [
  DashHtmlComponents.html_div(className = "card-body", [
  DashHtmlComponents.html_p(children=[""], id = "results_lbl"),
  ])
  ])
  ])
  ])
  ])
  ])
  # ,
  # DashCoreComponents.dcc_tab(label = "Assembled Tests", children = [
  # "hi3"
  # ]),
  ])


  Dash.callback!(app, [Dash.output("results_lbl","children"), Dash.output("results_card","className")], Dash.Input("assemble_btn", "n_clicks"), [Dash.State("start_temp_txt","value"),
  Dash.State("geom_temp_txt","value"),
  Dash.State("max_time_txt","value"),
  Dash.State("n_item_StatsBase.sample_txt","value"),
  Dash.State("n_test_StatsBase.sample_txt","value"),
  Dash.State("max_conv_txt","value"),
  Dash.State("verbosity_txt","value"),
  Dash.State("opt_feas_txt","value"),
  Dash.State("n_fill_txt","value"),
  Dash.State("feas_nh_txt","value"),
  Dash.State("opt_nh_txt","value"),
  Dash.State("optimizer_constructor_txt","value"),
  Dash.State("optimizer_attributes_txt","value"),
  Dash.State("starting_design_txt","value"),
  Dash.State("solver_txt","value"),
  Dash.State("plots_out_txt","value"),
  Dash.State("results_folder_txt","value")]) do n_clicks, start_temp, geom_temp, max_time, n_item_sample, n_test_sample, max_conv, verbosity, opt_feas, n_fill, feas_nh, opt_nh, optimizer_constructor, optimizer_attributes, starting_design, solver, plots_out, results_folder
    if n_clicks > 0
      output=""
      try
        results_folder="RESULTS"
        if isfile(starting_design)
          starting_design=Matrix{Float64}(CSV.read(starting_design, delim = ";", header = false))
        else
          starting_design=Matrix{Float64}(undef,0,0)
        end
        println(" "
        ,start_temp," "
        ,geom_temp," "
        ,max_time," "
        ,results_folder," "
        ,n_item_sample," "
        ,n_test_sample," "
        ,max_conv," "
        ,verbosity," "
        ,opt_feas," "
        ,opt_feas," "
        ,n_fill," "
        ,feas_nh," "
        ,opt_nh," ")
        message = assemble!(ATAmodel, solver = solver,
        starting_design = starting_design,
        start_temp = Float64(start_temp),
        geom_temp = Float64(geom_temp),
        max_time =  Float64(max_time),
        results_folder = results_folder,
        n_item_sample =  Int64(n_item_sample),
        n_test_sample = Int64(n_test_sample),
        max_conv = Int64(max_conv),
        verbosity = Int64(verbosity),
        opt_feas = Float64(opt_feas) ,
        n_fill = Int64(n_fill),
        feas_nh = Int64(feas_nh),
        opt_nh = Int64(opt_nh),
        optimizer_attributes = eval(Meta.parse(optimizer_attributes)),
        optimizer_constructor = optimizer_constructor #eval(Meta.parse(string(optimizer_constructor,".Optimizer")))
        )
        plots_out = (plots_out == "yes") ? true : false
        if ATAmodel.settings.n_FS!=ATAmodel.settings.n_items
          print_results!(ATAmodel, group_by_fs = true, plots_out = plots_out, results_folder = results_folder)
        else
          print_results!(ATAmodel, results_folder = results_folder, plots_out = plots_out)
        end
        output = read(string(results_folder,"/Results.txt"), String)
      catch e
        println(showerror,e)
      end
      return output, "card h-100 bg-light"
    else
      return "", "card h-100 bg-light"
    end
  end





  Dash.callback!(app, [Dash.output("settings_lbl","children"), Dash.output("settings_card","className")], [Dash.Input("settings_btn", "n_clicks")], [Dash.State("settings_txt","value"), Dash.State("bank_txt","value"), Dash.State("bank_delim_txt","value")]) do n_clicks, settings_file, bank_file, bank_delim
    if n_clicks>0
      message=["success",""]
      try
        message = load_settings!(ATAmodel; settings_file = settings_file, bank_file = bank_file, bank_delim = bank_delim)
      catch e
        message[1] = "danger"
        message[2] = sprint(showerror,e)
      end
      return message[2], string("card h-100 text-white bg-", message[1])
    else
      return "", "card h-100 bg-light"
    end
  end



  Dash.callback!(app, [Dash.output("addFS_lbl","children"), Dash.output("addFS_card","className")], [Dash.Input("addFS_btn", "n_clicks")] , Dash.State("addFS_btn", "n_clicks")) do n_clicks, n_click_state
    if n_click_state > 0
      message=["success",""]
      try
        message = add_friends!(ATAmodel)
      catch e
        message[1] = "danger"
        message[2] = sprint(showerror,e)
      end
      return message[2], string("card h-100 text-white bg-", message[1])
    else
      return "", "card h-100 bg-light"
    end
  end

  Dash.callback!(app, [Dash.output("addES_lbl","children"), Dash.output("addES_card","className")], [Dash.Input("addES_btn", "n_clicks")], Dash.State("addES_btn", "n_clicks")) do n_clicks, n_click_state
    if n_click_state > 0
      message=["success",""]
      try
        message = add_enemies!(ATAmodel)
      catch e
        message[1] = "danger"
        message[2] = sprint(showerror,e)
      end
      return message[2], string("card h-100 text-white bg-", message[1])
    else
      return "", "card h-100 bg-light"
    end
  end

  Dash.callback!(app, [Dash.output("addConstraints_lbl","children"), Dash.output("addConstraints_card","className")], [Dash.Input("addConstraints_btn", "n_clicks")], [Dash.State("constraints_file_txt","value"), Dash.State("constraints_delim_txt","value")]) do n_clicks, constraints_file, constraints_delim
    if n_clicks > 0
      message=["success",""]
      try
        message = add_constraints!(ATAmodel; constraints_file = constraints_file, constraints_delim = constraints_delim)
      catch e
        message[1] = "danger"
        message[2] = sprint(showerror,e)
      end
      return message[2], string("card h-100 text-white bg-", message[1])
    else
      return "", "card h-100 bg-light"
    end
  end


  Dash.callback!(app, [Dash.output("ol_lbl","children"), Dash.output("ol_card","className")], [Dash.Input("addOl_btn", "n_clicks")], [ Dash.State("ol_txt","value"), Dash.State("ol_delim_txt","value")]) do n_clicks, ol_file, ol_delim
    if n_clicks > 0
      message=["success",""]
      try
        message = add_overlap!(ATAmodel; overlap_file = ol_file, overlap_delim = ol_delim)
      catch e
        message[1] = "danger"
        message[2] = sprint(showerror,e)
      end
      return message[2], string("card h-100 text-white bg-", message[1])
    else
      return "", "card h-100 bg-light"
    end
  end


  Dash.callback!(app, [Dash.output("addExpScore_lbl","children"), Dash.output("addExpScore_card","className")], [Dash.Input("addExpScore_btn", "n_clicks")]) do n_clicks
    if n_clicks > 0
      message=["success",""]
      try
        message = add_exp_score!(ATAmodel)
      catch e
        message[1] = "danger"
        message[2] = sprint(showerror,e)
      end
      return message[2], string("card h-100 text-white bg-", message[1])
    else
      return "", "card h-100 bg-light"
    end
  end

  Dash.callback!(app, [Dash.output("group_lbl","children"), Dash.output("group_card","className")], [Dash.Input("group_btn", "n_clicks")]) do n_clicks
    if n_clicks > 0
      message=["success",""]
      try
        message = group_by_friends!(ATAmodel)
      catch e
        message[1] = "danger"
        message[2] = sprint(showerror,e)
      end
      return message[2], string("card h-100 text-white bg-", message[1])
    else
      return "", "card h-100 bg-light"
    end
  end

  Dash.callback!(app, [Dash.output("addObj_lbl","children"), Dash.output("addObj_card","className")], [Dash.Input("addObj_btn", "n_clicks")]) do n_clicks
    if n_clicks > 0
      message=["success",""]
      try
        message = add_obj_fun!(ATAmodel)
      catch e
        message[1] = "danger"
        message[2] = sprint(showerror,e)
      end
      return message[2], string("card h-100 text-white bg-", message[1])
    else
      return "", "card h-100 bg-light"
    end
  end

  Dash.callback!(app, [Dash.output("start_lbl","children"), Dash.output("restart_card" ,"className")], Dash.Input("start_btn", "n_clicks"), [Dash.State("folder_txt","value")]) do n_clicks, folder
    if n_clicks > 0
      message = ""
      try
        cd(folder)
        global ATAmodel = start_ATA()
      catch e
        message = sprint(showerror,e)
      end
      if message==""
        return "ATA model successfully initialized", "card h-100 text-white bg-success"
      else
        return message, "card h-100 text-white bg-danger"
      end
    else
      return "", "card h-100 bg-light"
    end
  end

  Dash.callback!(app, [Dash.output("n_item_StatsBase.sample_txt","value"), Dash.output("n_test_StatsBase.sample_txt","value")], Dash.Input("settings_card", "className")) do success
    if success == "card h-100 text-white bg-success"
      return maximum([ATAmodel.constraints[t].length_max for t=1:ATAmodel.settings.T]), ATAmodel.settings.T
    else
      return throw(Dash.PreventUpdate())
    end
  end
  #
  Dash.callback!(app, Dash.output("settings_btn", "n_clicks"), Dash.Input("start_lbl","children")) do starter_lbl
    if starter_lbl=="ATA model successfully initialized"
      return 0
    else
      return Dash.no_update()
    end
  end

  Dash.callback!(app, Dash.output("addFS_btn", "n_clicks"), Dash.Input("settings_lbl","children")) do starter_lbl
    if starter_lbl==""
      return 0
    else
      return Dash.no_update()
    end
  end

  Dash.callback!(app, Dash.output("addES_btn", "n_clicks"), Dash.Input("addFS_lbl","children")) do starter_lbl
    if starter_lbl==""
      return 0
    else
      return Dash.no_update()
    end
  end

  Dash.callback!(app, Dash.output("addConstraints_btn", "n_clicks"), Dash.Input("addES_lbl","children")) do starter_lbl
    if starter_lbl==""
      return 0
    else
      return Dash.no_update()
    end
  end

  Dash.callback!(app, Dash.output("addOl_btn", "n_clicks"), Dash.Input("addConstraints_lbl","children")) do starter_lbl
    if starter_lbl==""
      return 0
    else
      return Dash.no_update()
    end
  end

  Dash.callback!(app, Dash.output("addExpScore_btn", "n_clicks"), Dash.Input("ol_lbl","children")) do starter_lbl
    if starter_lbl==""
      return 0
    else
      return Dash.no_update()
    end
  end

  Dash.callback!(app, Dash.output("group_btn", "n_clicks"), Dash.Input("addExpScore_lbl","children")) do starter_lbl
    if starter_lbl==""
      return 0
    else
      return Dash.no_update()
    end
  end

  Dash.callback!(app, Dash.output("addObj_btn", "n_clicks"), Dash.Input("group_lbl","children")) do starter_lbl
    if starter_lbl==""
      return 0
    else
      return Dash.no_update()
    end
  end

  Dash.callback!(app, Dash.output("assemble_btn", "n_clicks"), Dash.Input("addObj_lbl","children")) do starter_lbl
    if starter_lbl==""
      return 0
    else
      return Dash.no_update()
    end
  end

  println("Started at localhost:8080")
  Dash.run_server(app, "0.0.0.0", 8080)
end
