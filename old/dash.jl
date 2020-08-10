using HTTP
using Dash
using DashHtmlComponents
using DashCoreComponents

external_stylesheets = ["https://stackpath.bootstrapcdn.com/bootstrap/4.5.0/css/bootstrap.min.css"]
external_scripts = [
Dict(
"src"=> "jquery/jquery-3.5.1.min.js"),
Dict(
"src"=> "jquery/js/bootstrap.min.js.map"),
Dict(
"src"=> "jquery/js/bootstrap.min.js"),
]
app = Dash.dash() #assets_folder="assets",external_scripts = external_scripts)#external_stylesheets = external_stylesheets, external_scripts = external_scripts)
#
app.layout =
# [html_nav([
#   html_div(className = "nav nav-tabs", id="nav-tab", role="tablist", [
#     html_a(;className="nav-item nav-link active", id="nav-home-tab", (Symbol("data-toggle")=>"tab"), href="#nav-home", role="tab", (Symbol("aria-controls")=>"nav-home"), (Symbol("aria-selected")=>"true")),
#     html_a(;className="nav-item nav-link", id="nav-profile-tab", (Symbol("data-toggle")=>"tab"), href="#nav-profile", role="tab", (Symbol("aria-controls")=>"nav-profile"), (Symbol("aria-selected")=>"false")),
#     html_a(;className="nav-item nav-link", id="nav-contact-tab", (Symbol("data-toggle")=>"tab"), href="#nav-contact", role="tab",(Symbol("aria-controls")=>"nav-contact"), (Symbol("aria-selected")=>"false"))
#   ])
# ])
# html_div(className = "tab-content", id="nav-tabContent", [
#   html_div(;className="tab-pane fade show active", id="nav-home", role="tabpanel", (Symbol("aria-labelledby")=>"nav-home-tab")),
#   html_div(;className="tab-pane fade show", id="nav-profile", role="tabpanel", (Symbol("aria-labelledby")=>"nav-profile-tab")),
#   html_div(;className="tab-pane fade show", id="nav-contact", role="tabpanel", (Symbol("aria-labelledby")=>"nav-contact-tab"))
# ])
# ])
#
dcc_tabs([

dcc_tab(label = "Create model", children = [
html_div( [
html_div( className="modal-dialog", role="document",[
html_div( className="modal-content", [
html_div( className="modal-body", [
html_div( className="spinner-border", role="status", [
html_span("Loading...", className="sr-only")
])
])
])
])
];
id="loader", className="modal fade", tabIndex="-1", role="dialog", (Symbol("aria-labelledby")=>"exampleModalLabel") ,(Symbol("aria-hidden")=>"true")),
html_div(className = "container-fluid", [
html_div(className = "row-md-12 justify-content-center m-3", [
html_div(className = "col-md-12 ", [
html_div(className = "card m-1 p-1", [
html_div(className = "card-body", [
html_div(className = "row justify-content-center m-1", [
html_div(className = "col-md-6 ", [
html_div(className = "row justify-content-center m-1", [
html_label("folder:", className="mr-1"),
dcc_input(id="folder_txt", value=pwd(), type = "text", className="form-control")
]),
html_div(className = "row justify-content-center m-1", [
html_button("Restart ATA"; id = "start_btn", n_clicks=0, className = "btn btn-primary")
])
]),
html_div(className = "col-md-6", [
html_div(id = "restart_card", className = "card p-2 m-1", [
html_div(className = "card-body", [
html_p("", id = "start_lbl"),
])
])
])
])
])
])
])
]),


html_div(className = "row d-md-flex justify-content-md-center m-3", [
#html_p("hidden",id ="settings_loader", style="display: none;"),
html_div(className = "col-md-12", [
html_div(className = "card m-1 p-1", [
html_div(className = "card-body", [
html_div(className = "row m-1", [
html_div(className = "col-md-6", [
html_div(className = "row m-1", [
html_div(className = "col-md-6", [
html_label("file name:"),
]),
html_div(className = "col-md-6", [
dcc_input(id="settings_txt", value="settingsATA.jl", className="form-control", type = "text")
])
]),
html_div(className = "row m-1", [
html_div(className = "col-md-6", [
html_label("item pool file:", id="pool_lbl"),
]),
html_div(className = "col-md-6", [
dcc_input(id="bank_txt", value="bank.csv", className="form-control", type = "text"),
])
]),
html_div(className = "row m-1", [
html_div(className = "col-md-6", [
html_label("bank delimitator:", id="delim_bank_lbl"),
]),
html_div(className = "col-md-6", [
dcc_input(id="bank_delim_txt", value=";", className="form-control", type = "text"),
])
]),
html_div(className = "row justify-content-center m-1", [
html_button("Load Settings"; id="settings_btn", n_clicks=0, className = "btn btn-primary")
]),
]),
html_div(className = "col-md-6", [
html_div(id = "settings_card", className = "card", [
html_div(className = "card-body", [
html_p("", id = "settings_lbl"),
])
])
])
])
])
])
])
]),


html_div(className = "row justify-content-center m-3", [
html_div(className = "col-md-12", [
html_div(className = "card m-1 p-1", [
html_div(className = "card-body", [
html_div(className = "row justify-content-center m-3", [
html_div(className = "col-md-6  align-self-center text-center", [
html_button("Add Friend Sets"; id="addFS_btn", n_clicks=0, className = "btn btn-primary")
]),
html_div(className = "col-md-6", [
html_div(id = "addFS_card", className = "card", [
html_div(className = "card-body", [
html_p("", id = "addFS_lbl"),
])
])
])
])
])
])
])
]),


html_div(className = "row justify-content-center m-3", [
html_div(className = "col-md-12", [
html_div(className = "card m-1 p-1", [
html_div(className = "card-body", [
html_div(className = "row justify-content-center m-1", [
html_div(className = "col-md-6 align-self-center text-center", [
html_button("Add Enemy Sets"; id="addES_btn", n_clicks=0, className = "btn btn-primary")
]),
html_div(className = "col-md-6", [
html_div(id = "addES_card", className = "card", [
html_div(className = "card-body", [
html_p("", id = "addES_lbl"),
])
])
])
])
])
])
])
]),


html_div(className = "row justify-content-center m-3", [
html_div(className = "col-md-12", [
html_div(className = "card m-1 p-1", [
html_div(className = "card-body", [
html_div(className = "row justify-content-center m-1", [
html_div(className = "col-md-6", [
html_div(className = "row justify-content-center m-1", [
html_div(className = "col-md-6", [
html_label("file name:"),
]),
html_div(className = "col-md-6", [
dcc_input(id="constraints_file_txt", value="Constraints.csv", className="form-control", type = "text")
])
]),
html_div(className = "row justify-content-center m-1", [
html_div(className = "col-md-6", [
html_label("delimitator:", id="constraints_delim_lbl"),
]),
html_div(className = "col-md-6", [
dcc_input(id="constraints_delim_txt", value=";", className="form-control", type = "text"),
])
]),
html_div(className = "row justify-content-center m-1", [
html_button("Add Constraints"; id="addConstraints_btn", n_clicks=0, className = "btn btn-primary")
]),
]),
html_div(className = "col-md-6", [
html_div(id = "addConstraints_card", className = "card", [
html_div(className = "card-body", [
html_p("", id = "addConstraints_lbl"),
])
])
])
])
])
])
])
]),


html_div(className = "row justify-content-center m-3", [
html_div(className = "col-md-12", [
html_div(className = "card m-1 p-1", [
html_div(className = "card-body", [
html_div(className = "row justify-content-center m-1", [
html_div(className = "col-md-6", [
html_div(className = "row justify-content-center m-1", [
html_div(className = "col-md-6", [
html_label("file name:"),
]),
html_div(className = "col-md-6", [
dcc_input(id="ol_txt", value="Overlap Matrix.csv", className="form-control", type = "text"),
])
]),
html_div(className = "row justify-content-center m-1", [
html_div(className = "col-md-6", [
html_label("matrix delimitator:", id="ol_delim_lbl"),
]),
html_div(className = "col-md-6", [
dcc_input(id="ol_delim_txt", value=";", className="form-control", type = "text"),
])
]),
html_div(className = "row justify-content-center m-1", [
html_button("Add Overlap Matrix"; id="addOl_btn", n_clicks=0, className = "btn btn-primary")
]),
]),
html_div(className = "col-md-6", [
html_div(id = "ol_card",className = "card m-1 p-1", [
html_div(className = "card-body", [
html_p("", id = "ol_lbl"),
])
])
])
])
])
])
])
]),

html_div(className = "row justify-content-center m-3", [
    html_div(className = "col-md-12 align-self-center text-center", [
        html_div(className = "card m-1 p-1", [
            html_div(className = "card-body", [
                html_div(className = "row justify-content-center m-1", [
                    html_div(className = "col-md-6 align-self-center text-center", [
                        html_button("Add Expected Score"; id="addExpScore_btn", n_clicks=0, className = "btn btn-primary")
                        ]),
                        html_div(className = "col-md-6 align-self-center text-center", [
                            html_div(id = "addExpScore_card", className = "card m-1 p-1", [
                                html_div(className = "card-body", [
                                    html_p("", id = "addExpScore_lbl"),
                            ])
                        ])
                    ])
                ])
            ])
        ])
    ])
]),


html_div(className = "row justify-content-center m-3", [
    html_div(className = "col-md-12 align-self-center text-center", [
        html_div(className = "card m-1 p-1", [
            html_div(className = "card-body", [
                html_div(className = "row justify-content-center m-1", [
                    html_div(className = "col-md-6 align-self-center text-center", [
                        html_button("Group By Friend Sets"; id="group_btn", n_clicks=0, className = "btn btn-primary")
                    ]),
                    html_div(className = "col-md-6 align-self-center text-center", [
                        html_div(id = "group_card", className = "card m-1 p-1", [
                            html_div(className = "card-body", [
                                html_p("", id = "group_lbl"),
                            ])
                        ])
                    ])
                ])
            ])
        ])
    ])
]),


html_div(className = "row justify-content-center m-3", [
    html_div(className = "col-md-12 align-self-center text-center", [
        html_div(className = "card m-1 p-1", [
            html_div(className = "card-body", [
                html_div(className = "row justify-content-center m-1", [
                    html_div(className = "col-md-6 align-self-center text-center", [
                        html_button("Add objective function"; id="addObj_btn", n_clicks=0, className = "btn btn-primary")
                    ]),
                    html_div(className = "col-md-6 align-self-center text-center", [
                        html_div(id = "addObj_card", className = "card m-1 p-1", [
                            html_div(className = "card-body", [
                                html_p("", id = "addObj_lbl"),
                            ])
                        ])
                    ])
                ])
            ])
        ])
    ])
])

])
]),
dcc_tab(label = "Assemble!", children = [
    html_div(className = "container-fluid", [
        html_div(className = "row justify-content-center m-3", [
            html_div(className = "col-md-12 align-self-center text-center", [
                html_div(className = "card m-1 p-1", [
                    html_div(className = "card-body", [
                        html_div(className = "row justify-content-center m-1", [
                            html_div(className = "col-md-6", [
                            html_h1("SimAn solver"),
                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                html_label("Maximum time (in seconds):")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="100",className="form-control", type="text", id="max_time_txt")
                                ]),
                            ]),

                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                html_label("Starting temperature:")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="0.1",className="form-control", type="text", id="start_temp_txt")
                                ]),
                            ]),
                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                    html_label("Geometric temp decreasing:")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="0.1",className="form-control", type="text", id="geom_temp_txt")
                                ]),
                            ]),
                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                    html_label("Items to sample:")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="1",className="form-control", type="text", id="n_item_sample_txt")

                                ]),
                            ]),
                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                html_label("Tests to sample:")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="1",className="form-control", type="text", id="n_test_sample_txt")
                                ]),
                            ]),
                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                html_label("Minimum convergence:")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="10",className="form-control", type="text", id="conv_max_txt")

                                ]),
                            ]),
                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                html_label("Verbosity (1=minimal, 2=detailed):")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="2",className="form-control", type="text", id="verbosity_txt")
                                ]),
                            ]),
                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                html_label("Optimality/Feasibility balancer:")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="0.1", className="form-control", type="text", id="opt_feas_txt")
                                ]),
                            ]),
                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                html_label("#Fill-up phases:")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="1",className="form-control", type="text", id="n_fill_txt")
                                ]),
                            ]),
                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                html_label("#Feasibility neighbourhoods:")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="0",className="form-control", type="text", id="feas_nh_txt")
                                ]),
                            ]),
                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                html_label("#Optimality neighbourhoods:")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="100",className="form-control", type="text", id="opt_nh_txt")
                                ]),
                            ]),
                            ]),
                            html_div(className = "col-md-6", [
                            html_h1("JuMP solver"),
                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                html_label("Optimizer constructor:")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="CPLEX.Optimizer",className="form-control", type="text", id="optimizer_constructor_txt")
                                ]),
                            ]),
                            html_div(className = "row justify-content-center m-1", [
                                html_div(className = "col-md-6", [
                                html_label("Optimizer attributes:")
                                ]),
                                html_div(className = "col-md-6", [
                                dcc_input(value="[(\"tm_lim\",1000)]",className="form-control", type="text", id="optimizer_attributes_txt")
                                ]),
                            ]),

                            ]),
                        ]),
                        html_div(className = "row justify-content-center m-1", [
                            html_div(className = "col-md-6", [
                            html_label("Starting design file:")
                            ]),
                            html_div(className = "col-md-6", [
                            dcc_input(value="Design.csv",className="form-control", type="text",id="starting_design_txt")
                            ]),
                        ]),
                        html_div(className = "row justify-content-center m-1", [
                            html_div(className = "col-md-6", [
                            html_label("Results folder name:")
                            ]),
                            html_div(className = "col-md-6", [
                            dcc_input(value="RESULTS",className="form-control", type="text", id="results_folder_txt")
                            ]),
                        ]),
                        html_div(className = "row justify-content-center m-1", [
                            html_div(className = "col-md-6", [
                            html_label("Print plots (yes or no):")
                            ]),
                            html_div(className = "col-md-6", [
                            dcc_input(value="yes",className="form-control", type="text", id="plots_out_txt")
                            ]),
                        ]),
                        html_div(className = "row justify-content-center m-1", [
                            html_div(className = "col-md-6", [
                            html_label("Solver (siman or jumpATA):")
                            ]),
                            html_div(className = "col-md-6", [
                            dcc_input(value="siman",className="form-control", type="text", id="solver_txt")
                            ]),
                        ]),
                        html_div(className = "row justify-content-center m-1", [
                            html_button("Assemble!"; id="assemble_btn", n_clicks=0, className = "btn btn-success")
                        ])
                    ])
                ])
            ])
        ])
    ])
]),

dcc_tab(label = "Results", children = [
html_div(className = "container-fluid", [
html_div(className = "row-md-12 justify-content-center m-3", [
html_div(className = "col-md-12", [
    html_div(id = "results_card", className = "card m-1 p-1", [
        html_div(className = "card-body", [
            html_p("", id = "results_lbl"),
        ])
    ])
])
])
])
])
# ,
# dcc_tab(label = "Assembled Tests", children = [
# "hi3"
# ]),
])


callback!(app, [Dash.Output("results_lbl","children"), Dash.Output("results_card","className")], Dash.Input("assemble_btn", "n_clicks"), [Dash.State("start_temp_txt","value"),
Dash.State("geom_temp_txt","value"),
Dash.State("max_time_txt","value"),
Dash.State("n_item_sample_txt","value"),
Dash.State("n_test_sample_txt","value"),
Dash.State("conv_max_txt","value"),
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
Dash.State("results_folder_txt","value")]) do n_clicks, start_temp, geom_temp, max_time, n_item_sample, n_test_sample, conv_max, verbosity, opt_feas, n_fill, feas_nh, opt_nh, optimizer_constructor, optimizer_attributes, starting_design, solver, plots_out, results_folder
    if n_clicks > 0
        output=""
        try
            results_folder="RESULTS"
            if isfile(starting_design)
                starting_design=Matrix{Float64}(CSV.read(starting_design, delim = ";", header = false))
            else
                starting_design=Matrix{Float64}(undef,0,0)
            end
            message = Assemble!(ATAmodel, solver = solver,
            starting_design = starting_design,
            start_temp = parse(Float64, start_temp),
            geom_temp = parse(Float64, geom_temp),
            max_time = parse(Int, max_time),
            results_folder = results_folder,
            n_item_sample = parse(Int, n_item_sample),
            n_test_sample = parse(Int, n_test_sample),
            conv_max = parse(Int, conv_max),
            verbosity = parse(Int, verbosity),
            opt_feas = parse(Float64, opt_feas) ,
            n_fill = parse(Int, n_fill) ,
            feas_nh = parse(Int, feas_nh),
            opt_nh = parse(Int, opt_nh),
            optimizer_attributes =  eval(Meta.parse(optimizer_attributes)),
            optimizer_constructor = ( () ->eval(Meta.parse(optimizer_constructor)))
            )
            plots_out = (plots_out == "yes") ? true : false
            if ATAmodel.Settings.nFS!=ATAmodel.Settings.nItems
                PrintResults(ATAmodel, GroupByFS= true, plots_out = plots_out)
            else
                PrintResults(ATAmodel, results_folder = results_folder, plots_out = plots_out)
            end
            output = read(string(results_folder,"/Results.txt"), String)
        catch e
            println(showerror,e)
        end
            return output, string("card")
    else
        throw(PreventUpdate())
    end
end

global ATAmodel = StartATA()

callback!(app, [Dash.Output("start_lbl","children"), Dash.Output("restart_card" ,"className")],  Dash.Input("start_btn", "n_clicks"), [Dash.State("start_btn","n_clicks"),  Dash.State("folder_txt","value")]) do trigger, n_clicks, folder
    if n_clicks > 0
        message = ""
        try
            cd(folder)
            global ATAmodel = StartATA()
        catch e
            message = sprint(showerror,e)
        end

        if message==""
            return "ATA model successfully initialized", "card text-white bg-success"
        else
            return message, "card text-white bg-danger"
        end
    else
        throw(PreventUpdate())
    end
end
#


callback!(app, [Dash.Output("settings_lbl","children"), Dash.Output("settings_card","className")], [Dash.Input("settings_btn", "n_clicks"),Dash.Input("start_lbl","children")], [Dash.State("settings_txt","value"), Dash.State("bank_txt","value"), Dash.State("bank_delim_txt","value"),Dash.State("start_lbl","children")]) do n_clicks, restarted_trigger,  settings_file, bank_file, bank_delim, restarted
        if n_clicks>0
            message=["success",""]
            try
                message = LoadSettings!(ATAmodel; settings_file = settings_file, bank_file = bank_file, bank_delim = bank_delim)
            catch e
                message[1] = "danger"
                message[2] = sprint(showerror,e)
            end
            return message[2], string("card text-white bg-", message[1])
        else
            throw(PreventUpdate())
        end
end


callback!(app, [Dash.Output("addFS_lbl","children"), Dash.Output("addFS_card","className")], Dash.Input("addFS_btn", "n_clicks")) do  n_clicks
    if n_clicks > 0
        message=["success",""]
        try
            message = AddFriendSets!(ATAmodel)
        catch e
            message[1] = "danger"
            message[2] = sprint(showerror,e)
        end
            return message[2], string("card text-white bg-", message[1])
        else
            throw(PreventUpdate())
        end
end

callback!(app, [Dash.Output("addES_lbl","children"), Dash.Output("addES_card","className")], Dash.Input("addES_btn", "n_clicks")) do  n_clicks
    if n_clicks > 0
        message=["success",""]
        try
            message = AddEnemySets!(ATAmodel)
        catch e
            message[1] = "danger"
            message[2] = sprint(showerror,e)
        end
            return message[2], string("card text-white bg-", message[1])
    else
        throw(PreventUpdate())
    end
end

callback!(app, [Dash.Output("addConstraints_lbl","children"), Dash.Output("addConstraints_card","className")], Dash.Input("addConstraints_btn", "n_clicks"),  [Dash.State("constraints_file_txt","value"), Dash.State("constraints_delim_txt","value")]) do n_clicks, constraints_file, constraints_delim
    if n_clicks > 0
        message=["success",""]
        try
            message = AddConstr!(ATAmodel; constraints_file = constraints_file, constraints_delim = constraints_delim)
        catch e
            message[1] = "danger"
            message[2] = sprint(showerror,e)
        end
            return message[2], string("card text-white bg-", message[1])
    else
        throw(PreventUpdate())
    end
end


callback!(app, [Dash.Output("ol_lbl","children"), Dash.Output("ol_card","className")], Dash.Input("addOl_btn", "n_clicks"),  [ Dash.State("ol_txt","value"), Dash.State("ol_delim_txt","value")]) do n_clicks, ol_file, ol_delim
    if n_clicks > 0
        message=["success",""]
        try
            message = AddOverlaps!(ATAmodel; overlap_file = ol_file, overlap_delim = ol_delim)
        catch e
            message[1] = "danger"
            message[2] = sprint(showerror,e)
        end
            return message[2], string("card text-white bg-", message[1])
    else
        throw(PreventUpdate())
    end
end


callback!(app, [Dash.Output("addExpScore_lbl","children"), Dash.Output("addExpScore_card","className")], Dash.Input("addExpScore_btn", "n_clicks")) do  n_clicks
    if n_clicks > 0
        message=["success",""]
        try
            message = AddExpScore!(ATAmodel)
        catch e
            message[1] = "danger"
            message[2] = sprint(showerror,e)
        end
            return message[2], string("card text-white bg-", message[1])
    else
        throw(PreventUpdate())
    end
end

callback!(app, [Dash.Output("group_lbl","children"), Dash.Output("group_card","className")], Dash.Input("group_btn", "n_clicks")) do  n_clicks
    if n_clicks > 0
        message=["success",""]
        try
            message = GroupByFriendSet!(ATAmodel)
        catch e
            message[1] = "danger"
            message[2] = sprint(showerror,e)
        end
            return message[2], string("card text-white bg-", message[1])
    else
        throw(PreventUpdate())
    end
end

callback!(app, [Dash.Output("addObj_lbl","children"), Dash.Output("addObj_card","className")], Dash.Input("addObj_btn", "n_clicks")) do n_clicks
    if n_clicks > 0
        message=["success",""]
        try
            message = AddObjFun!(ATAmodel)
        catch e
            message[1] = "danger"
            message[2] = sprint(showerror,e)
        end
            return message[2], string("card text-white bg-", message[1])
    else
        throw(PreventUpdate())
    end
end


handler = Dash.make_handler(app)
println("Started at localhost:8080")
HTTP.serve(handler, HTTP.Sockets.localhost, 8080)
