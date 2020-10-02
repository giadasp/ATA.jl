"""
	plot_results(ATAmodel; group_by_fs = false, results_folder = "PLOTS")

# Description

Plot the ICFs and TIFs of the assembled tests.

# Arguments

- **`ATAmodel::Model`** : Required. The model built with `ATA` fuctions, `ATAmodel.design` matrix must be `IxT` or `nFSxT` if the items are grouped by friend sets. 
- **`group_by_fs`** : Optional. Default: `false`. Set to `true` if items have been grouped by friend sets by [`group_by_friends!`](#ATA.group_by_friends!-Tuple{ATA.Model}).
- **`results_folder`** : Optional. Default: "PLOTS". The folder in which the output is stored.
"""
function plot_results(ATAmodel; group_by_fs = false, results_folder = "PLOTS")

    if size(ATAmodel.output.design, 1) > 0
        n_items = ATAmodel.settings.n_items
        length_max = [ATAmodel.constraints[t].length_max for t = 1:T]
        n_FS = size(ATAmodel.settings.FS.items, 1)

        if n_FS < ATAmodel.settings.n_items
            n = n_FS * T
        else
            n = n_items * T
        end

        if ATAmodel.obj.type == "CC"
            JLD2.@load "OPT/IIF_CC.jld2" IIF
            JLD2.@load "OPT/ICF_CC.jld2" ICF
            IIF_CC = copy(IIF)
            ICF_CC = copy(ICF)
        end

        if ATAmodel.obj.type == "MAXIMIN" || ATAmodel.obj.type == "CC"
            JLD2.@load "OPT/IIF.jld2" IIF
            JLD2.@load "OPT/ICF.jld2" ICF
        end

        if group_by_fs == true
            design = reshape(ATAmodel.output.design, n_FS, T)
            new_design = zeros(Float64, n_items, T)

            for t = 1:T
                new_design[
                    vcat(ATAmodel.settings.FS.items[findall(
                        design[:, t] .== one(Float64),
                    )]...),
                    t,
                ] .= one(Float64)
            end

            design = new_design
            DelimitedFiles.writedlm(string(results_folder, "/designItemLevel.csv"), design)
        else
            design = reshape(ATAmodel.output.design, n_items, T)
        end

        #TIF e ICF
        if ATAmodel.obj.type == "MAXIMIN" || ATAmodel.obj.type == "CC"
            if isfile("simPool.csv")
                simPool = CSV.read("simPool.csv")
            else
                simPool = Float64[]
            end

            if plots_out == true
                ThetasPlot = collect(range(-4, stop = 4, length = 101)) #nqp values in interval/r/n",
                IIFplot = item_info(
                    ATAmodel.settings.IRT.parameters,
                    ThetasPlot,
                    model = ATAmodel.settings.IRT.model,
                )
                ICFplot = item_char(
                    ATAmodel.settings.IRT.parameters,
                    ThetasPlot,
                    model = ATAmodel.settings.IRT.model,
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

                ATAPlot.PlotATA(
                    ATAmodel,
                    IIFf,
                    ICFf,
                    design;
                    simPool = simPool,
                    results_folder = results_folder,
                )

                if ATAmodel.obj.type == "CC"
                    ATAPlot.PlotATA_CC(
                        ATAmodel,
                        IIFf,
                        ICFf,
                        design;
                        simPool = simPool,
                        results_folder = results_folder,
                    )
                end

            end

        end

    else
        println("No solution found")
    end

end
