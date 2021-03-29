function jump!(
    ata_model::AbstractModel;
    starting_design = Matrix{Float64}(undef, 0, 0),
    results_folder = "results",
    optimizer_constructor = "GLPK",
    optimizer_attributes = [("tm_lim", 500000), ("msg_lev", 3)],
    kwargs...,
)
    message = ["", ""]
    if !isdir(results_folder)
        mkdir(results_folder)
    else
        message[2] =
            message[2] * string(
                "- There is already a folder with this name, files in ",
                results_folder,
                " will be overwritten.\n",
            )
    end
    n_items = ata_model.settings.n_items
    if ata_model.obj.name in ["maximin", "minimax", "soyster_maximin", "de_jong_maximin"]
        if any([size(ata_model.obj.cores[t].IIF, 1) > 0 for t = 1:ata_model.settings.T])
            IIF = [ata_model.obj.cores[t].IIF for t = 1:ata_model.settings.T]
            message[2] =
                message[2] *
                string("- Assembling tests with ", ata_model.obj.name, " objective...\n")
        else
            message[1] = "danger"
            message[2] =
                message[2] * "- IIFs have not been computed. Run add_obj_fun!() first.\n"
            push!(ata_model.output.infos, message)
            return nothing
        end
    elseif ata_model.obj.name == "cc_maximin"
        message[1] = "danger"
        message[2] =
            message[2] *
            "- You must use the siman solver to assemble tests with chance-constrained MAXIMIN objective function.\n"
        push!(ata_model.output.infos, message)
        return nothing
    elseif ata_model.obj.name == ""
        IIF = []
        message[2] = message[2] * "- Assembling tests with NO objective function...\n"
    end

    if ata_model.settings.n_fs == 0
        ata_model.settings.n_fs = ata_model.settings.n_items
    end

    ################################################################################
    #                                overlap
    ################################################################################	
    overlap_matrix = ata_model.settings.ol_max
    nPairs = 0
    if ata_model.settings.to_apply[3]
        if size(overlap_matrix, 1) > 0
            Pairs_t = _combinations(ata_model.settings.T)
            nPairs_t = size(Pairs_t, 1)
            Pairs = _combinations(ata_model.settings.T)
            nPairs = size(Pairs, 1)
            ol_max = Array{Int64,1}(undef, nPairs)
            fInd = [Pairs[pair][1] for pair = 1:nPairs]
            fIndFirst = [Pairs[pair][2] for pair = 1:nPairs]
            for pair = 1:nPairs
                ol_max[pair] = overlap_matrix[fInd[pair], fIndFirst[pair]]
            end
        end
    end
    #Group IIFs by friend set
    if ata_model.obj.name in ["maximin", "soyster_maximin", "de_jong_maximin", "minimax"]
        IIF_new = [zeros(Float64, 0, 0) for t = 1:ata_model.settings.T]
        if ata_model.settings.n_fs != ata_model.settings.n_items
            #group IIFs
            for t = 1:ata_model.settings.T
                if size(IIF, 1) > 0
                    IIF_new[t] =
                        Matrix{Float64}(undef, size(IIF[t], 1), ata_model.settings.n_fs)
                end
            end
            for fs = 1:ata_model.settings.n_fs
                for t = 1:ata_model.settings.T
                    if size(IIF[t], 1) > 0
                        IIF_new[t][:, fs] =
                            sum(IIF[t][:, ata_model.settings.fs.items[fs]], dims = 2)
                    end
                end
            end
        else
            IIF_new = IIF
            ICF_new =
                [ata_model.constraints[t].expected_score.val for t = 1:ata_model.settings.T]
        end
    end
    # group expected score by friend set
    ICF_new = [zeros(Float64, 0, 0) for t = 1:ata_model.settings.T]
    if ata_model.settings.n_fs != ata_model.settings.n_items
        #group IIFs
        for t = 1:ata_model.settings.T
            if size(ata_model.constraints[t].expected_score.val, 1) > 0
                ICF_new[t] = zeros(
                    Float64,
                    size(ata_model.constraints[t].expected_score.val, 1),
                    ata_model.settings.n_fs,
                )
            end
        end
        for fs = 1:ata_model.settings.n_fs
            for t = 1:ata_model.settings.T
                if size(ata_model.constraints[t].expected_score.val, 1) > 0
                    ICF_new[t][:, fs] = sum(
                        ata_model.constraints[t].expected_score.val[
                            :,
                            ata_model.settings.fs.items[fs],
                        ],
                        dims = 2,
                    )
                end
            end
        end
    else
        ICF_new =
            [ata_model.constraints[t].expected_score.val for t = 1:ata_model.settings.T]
    end
    ################################################################################
    #                                count constraints
    ################################################################################
    c = 1
    #length
    for t = 1:ata_model.settings.T
        c += ata_model.constraints[t].length_max .> 0 ? 1 : 0
        c += ata_model.constraints[t].length_min .> 0 ? 1 : 0
    end

    # item use
    if ata_model.settings.to_apply[1]
        c += size(ata_model.settings.iu.max, 1)
    end
    if ata_model.settings.to_apply[2]
        c += sum(ata_model.settings.iu.min .> 0)
    end
    #overlap old
    if ata_model.settings.to_apply[3]
        c += nPairs
    end
    #expected score
    for t = 1:ata_model.settings.T
        if size(ICF_new[t], 1) > 0
            for k = 1:size(ICF_new[t], 1)
                if ata_model.constraints[t].expected_score.min[k] > 0
                    c += 1
                end
                if ata_model.constraints[t].expected_score.max[k] < 1
                    c += 1
                end
            end
        end
    end

    #constraints
    for t = 1:ata_model.settings.T
        if size(ata_model.constraints[t].constr_A, 1) > 0
            c += size(ata_model.constraints[t].constr_A, 1)
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
        if ata_model.settings.n_fs != ata_model.settings.n_items
            if (
                size(starting_design, 1) != ata_model.settings.n_fs ||
                size(starting_design, 2) != ata_model.settings.T
            )
                message[1] = "danger"
                message[2] = message[2] * "- Starting design must be of size: (n_fs x T).\n"
                push!(ata_model.output.infos, message)
                return nothing
            end
        else
            if (
                size(starting_design, 1) != ata_model.settings.n_items ||
                size(starting_design, 2) != ata_model.settings.T
            )
                message[1] = "danger"
                message[2] = message[2] * "- Starting design must be of size: (I x T).\n"
                push!(ata_model.output.infos, message)
                return nothing
            end
        end
        if (any(starting_design != 0 && starting_design != 1))
            message[1] = "danger"
            message[2] = message[2] * "- Starting design must contain only 1s or 0s.\n"
            push!(ata_model.output.infos, message)
            return nothing
        end
    else
        starting_design = zeros(ata_model.settings.n_fs, ata_model.settings.T)
    end

    #decision variables
    JuMP.@variable(
        m,
        x[i = 1:ata_model.settings.n_fs, t = 1:ata_model.settings.T] <=
        Int(ata_model.settings.forced0[t][i]),
        Bin,
        start = starting_design[i, t]
    )

    if ata_model.settings.to_apply[3]
        if size(overlap_matrix, 1) > 0
            if maximum(overlap_matrix) > 0
                #Overlap Vars new
                JuMP.@variable(m, y[i = 1:ata_model.settings.n_fs, p = 1:nPairs], Bin)
            end
        end
    end

    # Item Use
    for i = 1:ata_model.settings.n_fs
        if ata_model.settings.to_apply[2]
            if ata_model.settings.iu.min[i] .> 0
                JuMP.@constraint(
                    m,
                    ata_model.settings.iu.min[i] -
                    sum(x[i, t] for t = 1:ata_model.settings.T) <= 0
                ) # z[c])
                c += 1
            end
        end
        if ata_model.settings.to_apply[1]
            JuMP.@constraint(
                m,
                sum(x[i, t] for t = 1:ata_model.settings.T) -
                ata_model.settings.iu.max[i] <= 0
            ) # z[c])
            c += 1
        end
    end
    if ata_model.settings.to_apply[3]
        if size(overlap_matrix, 1) > 0
            if maximum(overlap_matrix) > 0
                #overlap classic
                for p = 1:nPairs
                    JuMP.@constraint(
                        m,
                        sum(
                            y[i, p] * ata_model.settings.fs.counts[i] for
                            i = 1:ata_model.settings.n_fs
                        ) <= ol_max[p]
                    )
                    JuMP.@constraint(
                        m,
                        [i = 1:ata_model.settings.n_fs],
                        2 * y[i, p] <= x[i, fInd[p]] + x[i, fIndFirst[p]]
                    )
                    JuMP.@constraint(
                        m,
                        [i = 1:ata_model.settings.n_fs],
                        y[i, p] >= x[i, fInd[p]] + x[i, fIndFirst[p]] - 1
                    )
                end
            else
                #no overlap
                for p = 1:nPairs
                    for i = 1:ata_model.settings.n_fs
                        JuMP.@constraint(m, x[i, fInd[p]] + x[i, fIndFirst[p]] <= 1)
                    end
                end
            end
        end
    end
    #expected score
    for t = 1:ata_model.settings.T
        if size(ICF_new[t], 1) > 0
            for k = 1:size(ICF_new[t], 1)
                if ata_model.constraints[t].expected_score.min[k] > 0
                    JuMP.@constraint(
                        m,
                        sum(
                            x[i, t] * round(ICF_new[t][k, i]; digits = 3) for
                            i = 1:ata_model.settings.n_fs
                        ) >= round(
                            ata_model.constraints[t].expected_score.min[k] *
                            ata_model.constraints[t].length_min;
                            digits = 3,
                        )
                    )
                    c += 1
                end
                if ata_model.constraints[t].expected_score.max[k] < 1
                    JuMP.@constraint(
                        m,
                        sum(
                            x[i, t] * round(ICF_new[t][k, i]; digits = 3) for
                            i = 1:ata_model.settings.n_fs
                        ) <= round(
                            ata_model.constraints[t].expected_score.max[k] *
                            ata_model.constraints[t].length_max;
                            digits = 3,
                        )
                    )
                    c += 1
                end
            end
        end
    end

    #constraints
    for t = 1:ata_model.settings.T
        if size(ata_model.constraints[t].constr_A, 1) > 0
            for constr = 1:size(ata_model.constraints[t].constr_A, 1)
                JuMP.@constraint(
                    m,
                    sum(
                        x[i, t] * ata_model.constraints[t].constr_A[constr, i] for
                        i = 1:ata_model.settings.n_fs
                    ) <= ata_model.constraints[t].constr_b[constr] + 0
                ) # z[c])
                c += 1
            end
        end
    end


    ncons = copy(c) - 1

    if ata_model.obj.name in ["maximin", "soyster_maximin", "de_jong_maximin"]
        #Objective bound
        JuMP.@variable(m, w >= 0)
        for t = 1:ata_model.settings.T
            for k = 1:size(ata_model.obj.cores[t].points, 1)
                JuMP.@constraint(
                    m,
                    sum(
                        round(IIF_new[t][k, i]; digits = 4) * x[i, t] for
                        i = 1:ata_model.settings.n_fs
                    ) >= w
                )
            end
        end
        JuMP.@objective(m, Min, (-w))
    elseif ata_model.obj.name == "minimax"
        #Objective bound
        JuMP.@variable(m, epsilon >= 0)
        for t = 1:ata_model.settings.T
            for k = 1:size(ata_model.obj.cores[t].points, 1)
                JuMP.@constraint(
                    m,
                    sum(
                        round(IIF_new[t][k, i]; digits = 4) * x[i, t] for
                        i = 1:ata_model.settings.n_fs
                    ) - ata_model.obj.cores[t].targets[k] <= epsilon
                )
                JuMP.@constraint(
                    m,
                    sum(
                        round(IIF_new[t][k, i]; digits = 4) * x[i, t] for
                        i = 1:ata_model.settings.n_fs
                    ) - ata_model.obj.cores[t].targets[k] >= -epsilon
                )
            end
        end
        JuMP.@objective(m, Min, epsilon)
    end
    JuMP.optimize!(m)
    message[2] =
        message[2] *
        string("- The model has termination status:", JuMP.termination_status(m), "\n.")
    #println(string("The model has termination status:", JuMP.termination_status(m)))
    design = abs.(round.(JuMP.value.(x)))
    DelimitedFiles.writedlm(string(results_folder, "/design.csv"), design)
    ata_model.output.design = design
    JLD2.@save string(results_folder, "/ata_model.jld2") ata_model
    push!(ata_model.output.infos, ["success", message])
end
