"""
    plot_results(
        ata_model::AbstractModel;
        plots_folder = "plots",
    )

# Description

Plot the ICFs and TIFs of the assembled tests.

# Arguments

- **`ata_model::AbstractModel`** : Required. The model built with `ATA` fuctions, `ata_model.design` matrix must be `IxT` or `nfsxT` if the items are grouped by friend sets. 
- **`plots_folder`** : Optional. Default: "plots". The folder in which the output is stored.
"""
function plot_results(ata_model::AbstractModel; plots_folder = "plots")
    if !isdir(plots_folder)
        mkdir(plots_folder)
    else
        success!(
            ata_model,
            string(
                "There is already a folder with this name, files in ",
                plots_folder,
                " will be overwritten.",
            ),
        )
    end
    T = ata_model.settings.T
    if size(ata_model.output.design, 1) > 0
        n_items = ata_model.settings.n_items
        length_max = [ata_model.constraints[t].length_max for t = 1:T]
        n_fs = size(ata_model.settings.fs.items, 1)
        if n_fs < ata_model.settings.n_items
            n = n_fs * T
            design = reshape(ata_model.output.design, n_fs, T)
            new_design = zeros(Float64, n_items, T)
            for t = 1:T
                new_design[
                    vcat(
                        ata_model.settings.fs.items[findall(
                            design[:, t] .== one(Float64),
                        )]...,
                    ),
                    t,
                ] .= one(Float64)
            end
            design = new_design
            DelimitedFiles.writedlm(string(plots_folder, "/designItemLevel.csv"), design)
        else
            n = n_items * T
            design = reshape(ata_model.output.design, n_items, T)
        end

        if ata_model.obj.name == "cc_maximin"
            JLD2.@load "opt/IIF_CC.jld2" IIF
            JLD2.@load "opt/ICF_CC.jld2" ICF
            IIF_CC = copy(IIF)
            ICF_CC = copy(ICF)
        end

        if ata_model.obj.name in
           ["maximin", "soyster_maximin", "de_jong_maximin", "cc_maximin", "minimax"]
            JLD2.@load "opt/IIF.jld2" IIF
            JLD2.@load "opt/ICF.jld2" ICF
        end


        #TIF e ICF
        if ata_model.obj.name in
           ["maximin", "soyster_maximin", "de_jong_maximin", "cc_maximin", "minimax"]
            if isfile("sim_pool.csv")
                sim_pool = CSV.read("sim_pool.csv", DataFrames.DataFrame)
            else
                sim_pool = Float64[]
            end

            ThetasPlot = collect(range(-4, stop = 4, length = 101)) #nqp values in interval/r/n",
            IIFplot = item_info(
                ata_model.settings.irt.parameters,
                ThetasPlot,
                model = ata_model.settings.irt.model,
                parametrization = ata_model.settings.irt.parametrization,
                D = ata_model.settings.irt.D,
            )
            ICFplot = item_char(
                ata_model.settings.irt.parameters,
                ThetasPlot,
                model = ata_model.settings.irt.model,
                parametrization = ata_model.settings.irt.parametrization,
                D = ata_model.settings.irt.D,
            )[1][
                :,
                :,
                1,
            ]
            IIFf = Array{Float64,2}(undef, T, 101)
            ICFf = Array{Float64,2}(undef, T, 101)

            for t = 1:T
                IIFf[t, :] = [sum(IIFplot[findall(design[:, t] .> 0), i]) for i = 1:101]
                ICFf[t, :] = [sum(ICFplot[findall(design[:, t] .> 0), i]) for i = 1:101]
            end

            ATAPlot.plot_ATA(
                ata_model,
                IIFf,
                ICFf,
                design;
                sim_pool = sim_pool,
                plots_folder = plots_folder,
            )

            if ata_model.obj.name == "cc_maximin"
                ATAPlot.plot_ATA_CC(
                    ata_model,
                    IIFf,
                    ICFf,
                    design;
                    sim_pool = sim_pool,
                    plots_folder = plots_folder,
                )
            end

        end

    else
        println("No solution found")
    end

end
