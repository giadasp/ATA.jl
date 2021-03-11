

"""
group_by_friends!(ATAmodel::AbstractModel)

# Description

Group the items by friend sets once the friend sets have been added to the `ATA.AbstractModel` by [`add_friends!`](#ATA.add_friends!-Tuple{ATA.AbstractModel}) function and as **last operation** of the assembly.

# Arguments

- **`ATAmodel::AbstractModel`** : Required. The model built with `start_ATA()`, settings loaded by [`start_ATA`](#ATA.start_ATA) function and friend sets added by [`add_friends!`](#ATA.add_friends!-Tuple{ATA.AbstractModel}) function.
"""
function group_by_friends!(ATAmodel::AbstractModel) #last
    n_items = ATAmodel.settings.n_items
    #only works for categorical variables and item use, all the other contraitns need expansion by fs_items
    if ATAmodel.settings.n_fs == 0
        push!(
            ATAmodel.output.infos,
            ["danger", "No friend sets found, run add_friends!(model) before.\n"],
        )
        return nothing
    end
    if size(ATAmodel.constraints[1].constr_b, 1) == 0
        push!(
            ATAmodel.output.infos,
            ["danger", "No constraints to group, run add_constraints!(model) before.\n"],
        )
        return nothing
    else
        A = Vector{Matrix{Float64}}(undef, ATAmodel.settings.T)
        b = Vector{Vector{Float64}}(undef, ATAmodel.settings.T)
        A_new = Vector{Matrix{Float64}}(undef, ATAmodel.settings.T)
        for t = 1:ATAmodel.settings.T
            for t = 1:ATAmodel.settings.T
                A[t] = ATAmodel.constraints[t].constr_A
                b[t] = ATAmodel.constraints[t].constr_b
            end
        end
        n_fs = ATAmodel.settings.n_fs
        for t = 1:ATAmodel.settings.T
            A_new[t] = Matrix{Float64}(undef, size(A[t], 1), n_fs)
        end
        for fs = 1:n_fs
            for t = 1:ATAmodel.settings.T
                A_new[t][:, fs] = sum(A[t][:, ATAmodel.settings.fs.items[fs]], dims = 2)
            end
        end
        for t = 1:ATAmodel.settings.T
            DelimitedFiles.writedlm("OPT/A_$t.csv", A_new[t])
            DelimitedFiles.writedlm("OPT/b_$t.csv", b[t])
            ATAmodel.constraints[t].constr_A = A_new[t]
        end
        open("OPT/fs_items.jl", "w") do f
            write(f, "fs_items = Vector{Vector{Int64}}(undef, $n_fs)\n")
            for fs = 1:n_fs
                write(f, string("fs_items[$fs] =", ATAmodel.settings.fs.items[fs], "\n"))
            end
        end
        open("OPT/Settings.jl", "a") do f
            write(f, "group_by_fs = true\n\n")
        end
        #transform forced0
        x_forced0 = ATAmodel.settings.forced0
        x_forced0_new = Vector{Vector{Bool}}(undef, ATAmodel.settings.T)
        for v = 1:ATAmodel.settings.T
            x_forced0_new[v] = fill(true, n_fs)
            for i = 1:ATAmodel.settings.n_items
                if x_forced0[v][i] == false
                    for fs = 1:n_fs
                        if any(i .== ATAmodel.settings.fs.items[fs])
                            x_forced0_new[v][fs] = false
                        end
                    end
                end
            end
        end
        #update ATAmodel.forced0
        ATAmodel.settings.forced0 = x_forced0_new
        DelimitedFiles.writedlm("OPT/x_forced0.txt", x_forced0_new)
        #item use
        item_use_min = ATAmodel.settings.iu.min
        item_use_max = ATAmodel.settings.iu.max
        item_use_min_new = zeros(Int, n_fs)
        item_use_max_new = zeros(Int, n_fs)
        for fs = 1:n_fs
            item_use_min_new[fs] =
                Int(maximum(ATAmodel.settings.iu.min[ATAmodel.settings.fs.items[fs]]))
            item_use_max_new[fs] =
                Int(minimum(ATAmodel.settings.iu.max[ATAmodel.settings.fs.items[fs]]))
        end
        #enemy sets

        #update model
        ATAmodel.settings.iu.max = item_use_max_new
        ATAmodel.settings.iu.min = item_use_min_new
        open("OPT/Settings.jl", "a") do f
            write(f, "item_use_min = $item_use_min_new\n\n")
            write(f, "item_use_max = $item_use_max_new\n\n")
        end
        push!(
            ATAmodel.output.infos,
            ["success", string("- Grouped in ", n_fs, " Friend sets.\n")],
        )
    end
    return nothing
end
