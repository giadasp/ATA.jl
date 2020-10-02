function jumpATA!(
    ATAmodel::Model;
    starting_design = Matrix{Float64}(undef, 0, 0),
    results_folder = "RESULTS",
    optimizer_constructor = "GLPK",
    optimizer_attributes = [("tm_lim", 500000), ("msg_lev", 3)],
)
    message = ""
    if !(results_folder in readdir())
        mkdir(results_folder)
    else
        message *= string(
            "You have already a folder with this name, files in ",
            results_folder,
            " will be overwritten.\n",
        )
    end

    n_items = ATAmodel.settings.n_items
    if ATAmodel.obj.type == "MAXIMIN"
        if isfile("OPT/IIF.jld2")
            JLD2.@load "OPT/IIF.jld2" IIF
            message *= "- Assembling tests with MAXIMIN..."
        else
            message *= "No IIF.jld2 file in OPT folder, Run add_obj_fun!() first!"
            return message
        end
    elseif ATAmodel.obj.type == "CC"
        message *= "You must use the Simulated Annealing algorithm to assemble tests with CC objective function."
        return message
    elseif ATAmodel.obj.type == ""
        IIF = []
        message *= "Assembling tests with NO objective function..."
    end

    if ATAmodel.settings.n_FS == 0
        ATAmodel.settings.n_FS = ATAmodel.settings.n_items
    end

    ################################################################################
    #                                overlap
    ################################################################################	
    opMatrix = ATAmodel.settings.ol_max
    nPairs = 0
    if size(opMatrix, 1) > 0
        Pairs_t = _combinations(ATAmodel.settings.T)
        nPairs_t = size(Pairs_t, 1)
        Pairs = _combinations(ATAmodel.settings.T)
        nPairs = size(Pairs, 1)
        ol_max = Array{Int64,1}(undef, nPairs)
        fInd = [Pairs[pair][1] for pair = 1:nPairs]
        fIndFirst = [Pairs[pair][2] for pair = 1:nPairs]
        for pair = 1:nPairs
            ol_max[pair] = opMatrix[fInd[pair], fIndFirst[pair]]
        end
    end
    #Friend sets groups
    IIF_new = [zeros(Float64, 0, 0) for t = 1:ATAmodel.settings.T]
    ICF_new = [zeros(Float64, 0, 0) for t = 1:ATAmodel.settings.T]
    if ATAmodel.settings.n_FS != ATAmodel.settings.n_items
        #group IIFs
        for t = 1:ATAmodel.settings.T
            if size(IIF, 1) > 0
                IIF_new[t] = Matrix{Float64}(undef, size(IIF[t], 1), ATAmodel.settings.n_FS)
            end
            if size(ATAmodel.constraints[t].expected_score.val, 1) > 0
                ICF_new[t] = Matrix{Float64}(
                    undef,
                    size(ATAmodel.constraints[t].expected_score.val, 2),
                    ATAmodel.settings.n_FS,
                )
            end
        end
        for fs = 1:ATAmodel.settings.n_FS
            for t = 1:ATAmodel.settings.T
                if size(IIF[t], 1) > 0
                    IIF_new[t][:, fs] =
                        sum(IIF[t][:, ATAmodel.settings.FS.items[fs]], dims = 2)
                end
                if size(ATAmodel.constraints[t].expected_score.val, 1) > 0
                    ICF_new[t][:, fs] = sum(
                        ATAmodel.constraints[t].expected_score.val[
                            ATAmodel.settings.FS.items[fs],
                            :,
                        ],
                        dims = 1,
                    )
                end
            end
        end
    else
        IIF_new = IIF
        ICF_new = [ATAmodel.constraints[t].expected_score.val for t = 1:ATAmodel.settings.T]
    end
    ################################################################################
    #                                count constraints
    ################################################################################
    c = 1
    #length
    for t = 1:ATAmodel.settings.T
        c += ATAmodel.constraints[t].length_max .> 0 ? 1 : 0
        c += ATAmodel.constraints[t].length_min .> 0 ? 1 : 0
    end

    # item use
    c += size(ATAmodel.IU.max, 1)
    c += sum(ATAmodel.IU.min .> 0)

    #overlap old
    c += nPairs

    #expected score
    for t = 1:ATAmodel.settings.T
        if size(ICF_new[t], 1) > 0
            for k = 1:size(ICF_new[t], 1)
                if ATAmodel.constraints[t].expected_score.min[k] > 0
                    c += 1
                end
                if ATAmodel.constraints[t].expected_score.max[k] < 1
                    c += 1
                end
            end
        end
    end

    #constraints
    for t = 1:ATAmodel.settings.T
        if size(ATAmodel.constraints[t].constr_A, 1) > 0
            c += size(ATAmodel.constraints[t].constr_A, 1)
        end
    end
    ncons = copy(c) - 1
    ################################################################################
    #                                JuMP model
    ################################################################################

    m = JuMP.Model()
    if optimizer_constructor == "GLPK"
        add_GLPK!(m)
    elseif optimizer_constructor == "CPLEX"
        add_CPLEX!(m)
    elseif optimizer_constructor == "Gurobi"
        add_Gurobi!(m)
    elseif optimizer_constructor == "KNITRO"
        add_KNITRO!(m)
    elseif optimizer_constructor == "Cbc"
        add_Cbc!(m)
    elseif optimizer_constructor == "Xpress"
        add_Xpress!(m)
    elseif optimizer_constructor == "Juniper"
        add_Juniper!(m)
    elseif optimizer_constructor == "MosekTools"
        add_MosekTools!(m)
    elseif optimizer_constructor == "SCIP"
        add_SCIP!(m)
    end


    for (name, value) in optimizer_attributes
        JuMP.set_optimizer_attribute(m, name, value)
    end

    #starting design check
    if size(starting_design, 1) > 0
        if ATAmodel.settings.n_FS != ATAmodel.settings.n_items
            if (
                size(starting_design, 1) != ATAmodel.settings.n_FS ||
                size(starting_design, 2) != ATAmodel.settings.T
            )
                message *= "- Starting design must be of size: (n_FS x T).\n"
                return message
            end
        else
            if (
                size(starting_design, 1) != ATAmodel.settings.n_items ||
                size(starting_design, 2) != ATAmodel.settings.T
            )
                message *= "- Starting design must be of size: (n_items x T).\n"
                return message
            end
        end
        if (any(starting_design != 0 && starting_design != 1))
            message *= "- Starting design must contain only 1 or 0.\n"
            return message
        end
    else
        starting_design = zeros(ATAmodel.settings.n_FS, ATAmodel.settings.T)
    end

    #decision variables
    JuMP.@variable(
        m,
        x[i = 1:ATAmodel.settings.n_FS, t = 1:ATAmodel.settings.T] <=
        Int(ATAmodel.settings.forced0[t][i]),
        Bin,
        start = starting_design[i, t]
    )

    if size(opMatrix, 1) > 0
        if maximum(opMatrix) > 0
            #Overlap Vars new
            JuMP.@variable(m, y[i = 1:ATAmodel.settings.n_FS, p = 1:nPairs], Bin)
        end
    end

    # Item Use
    for i = 1:ATAmodel.settings.n_FS
        if ATAmodel.IU.min[i] .> 0
            JuMP.@constraint(
                m,
                ATAmodel.IU.min[i] - sum(x[i, t] for t = 1:ATAmodel.settings.T) <= 0
            ) # z[c])
            c += 1
        end
        JuMP.@constraint(
            m,
            sum(x[i, t] for t = 1:ATAmodel.settings.T) - ATAmodel.IU.max[i] <= 0
        ) # z[c])
        c += 1
    end
    if size(opMatrix, 1) > 0
        if maximum(opMatrix) > 0
            #overlap classic
            for p = 1:nPairs
                JuMP.@constraint(
                    m,
                    sum(
                        y[i, p] * ATAmodel.settings.FS.counts[i]
                        for i = 1:ATAmodel.settings.n_FS
                    ) <= ol_max[p]
                )
                JuMP.@constraint(
                    m,
                    [i = 1:ATAmodel.settings.n_FS],
                    2 * y[i, p] <= x[i, fInd[p]] + x[i, fIndFirst[p]]
                )
                JuMP.@constraint(
                    m,
                    [i = 1:ATAmodel.settings.n_FS],
                    y[i, p] >= x[i, fInd[p]] + x[i, fIndFirst[p]] - 1
                )
            end
        else
            #no overlap
            for p = 1:nPairs
                for i = 1:ATAmodel.settings.n_FS
                    JuMP.@constraint(m, x[i, fInd[p]] + x[i, fIndFirst[p]] <= 1)
                end
            end
        end
    end

    #expected score
    for t = 1:ATAmodel.settings.T
        if size(ICF_new[t], 1) > 0
            for k = 1:size(ICF_new, 1)
                if ATAmodel.constraints[t].expected_score.min[k] > 0
                    JuMP.@constraint(
                        m,
                        sum(
                            x[i, t] * round(ICF_new[t][k, i]; digits = 3)
                            for i = 1:ATAmodel.settings.n_FS
                        ) >= round(
                            ATAmodel.constraints[t].expected_score.min[k] *
                            ATAmodel.constraints[t].length_min;
                            digits = 3,
                        )
                    )
                    c += 1
                end
                if ATAmodel.constraints[t].expected_score.max[k] < 1
                    JuMP.@constraint(
                        m,
                        sum(
                            x[i, t] * round(ICF_new[t][k, i]; digits = 3)
                            for i = 1:ATAmodel.settings.n_FS
                        ) <= round(
                            ATAmodel.constraints[t].expected_score.max[k] *
                            ATAmodel.constraints[t].length_max;
                            digits = 3,
                        )
                    )
                    c += 1
                end
            end
        end
    end

    #constraints
    for t = 1:ATAmodel.settings.T
        if size(ATAmodel.constraints[t].constr_A, 1) > 0
            for constr = 1:size(ATAmodel.constraints[t].constr_A, 1)
                JuMP.@constraint(
                    m,
                    sum(
                        x[i, t] * ATAmodel.constraints[t].constr_A[constr, i]
                        for i = 1:ATAmodel.settings.n_FS
                    ) <= ATAmodel.constraints[t].constr_b[constr] + 0
                ) # z[c])
                c += 1
            end
        end
    end


    ncons = copy(c) - 1

    if ATAmodel.obj.type == "MAXIMIN"
        #Objective bound
        JuMP.@variable(m, w >= 0)
        for t = 1:ATAmodel.settings.T
            for k = 1:size(ATAmodel.obj.points[t], 1)
                JuMP.@constraint(
                    m,
                    sum(
                        round(IIF_new[t][k, i]; digits = 4) * x[i, t]
                        for i = 1:ATAmodel.settings.n_FS
                    ) >= w
                )
            end
        end
        JuMP.@objective(m, Min, (-w))
    end
    JuMP.optimize!(m)
    message *= string("The model has termination status:", JuMP.termination_status(m))
    println(string("The model has termination status:", JuMP.termination_status(m)))
    design = abs.(round.(JuMP.value.(x)))
    DelimitedFiles.writedlm(string(results_folder, "/design.csv"), design)
    ATAmodel.output.design = design
    JLD2.@save string(results_folder, "/ATAmodel.jld2") ATAmodel
    return message
end
