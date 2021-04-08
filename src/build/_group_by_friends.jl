

"""
    _group_by_friends!(ata_model::AbstractModel)

# Description

Group the items by friend sets once the friend sets have been added to the `ATA.AbstractModel` by [`add_friends!`](#ATA.add_friends!-Tuple{ATA.AbstractModel}) function and as **last operation** of the assembly.

# Arguments

- **`ata_model::AbstractModel`** : Required. The model built with `start_ata()`, settings loaded by [`start_ata`](#ATA.start_ata) function and friend sets added by [`add_friends!`](#ATA.add_friends!-Tuple{ATA.AbstractModel}) function.
"""
function _group_by_friends!(ata_model::AbstractModel) #last
    message = ["", ""]
    try
        n_items = ata_model.settings.n_items
        #only works for categorical variables and item use, all the other contraitns need expansion by fs_items
        if ata_model.settings.n_fs == 0
            push!(
                ata_model.output.infos,
                ["danger", "No friend sets found, run start_ata() before.\n"],
            )
            return nothing
        end
        if any([
            size(ata_model.constraints[t].constr_b, 1) > 0 for t = 1:ata_model.settings.T
        ])
            A = Vector{Matrix{Float64}}(undef, ata_model.settings.T)
            b = Vector{Vector{Float64}}(undef, ata_model.settings.T)
            A_new = Vector{Matrix{Float64}}(undef, ata_model.settings.T)
            for t = 1:ata_model.settings.T
                A[t] = ata_model.constraints[t].constr_A
                b[t] = ata_model.constraints[t].constr_b
            end
            n_fs = ata_model.settings.n_fs
            for t = 1:ata_model.settings.T
                A_new[t] = Matrix{Float64}(undef, size(A[t], 1), n_fs)
            end
            for fs = 1:n_fs
                for t = 1:ata_model.settings.T
                    A_new[t][:, fs] =
                        sum(A[t][:, ata_model.settings.fs.items[fs]], dims = 2)
                end
            end
            for t = 1:ata_model.settings.T
                DelimitedFiles.writedlm("opt/A_$t.csv", A_new[t])
                DelimitedFiles.writedlm("opt/b_$t.csv", b[t])
                ata_model.constraints[t].constr_A = A_new[t]
            end
            open("opt/fs_items.jl", "w") do f
                write(f, "fs_items = Vector{Vector{Int64}}(undef, $n_fs)\n")
                for fs = 1:n_fs
                    write(
                        f,
                        string("fs_items[$fs] =", ata_model.settings.fs.items[fs], "\n"),
                    )
                end
            end
            # open("opt/Settings.jl", "a") do f
            #     write(f, "group_by_fs = true\n\n")
            # end
            #transform forced0
            x_forced0 = ata_model.settings.forced0
            x_forced0_new = Vector{Vector{Bool}}(undef, ata_model.settings.T)
            for v = 1:ata_model.settings.T
                x_forced0_new[v] = fill(true, n_fs)
                for i = 1:ata_model.settings.n_items
                    if x_forced0[v][i] == false
                        for fs = 1:n_fs
                            if any(i .== ata_model.settings.fs.items[fs])
                                x_forced0_new[v][fs] = false
                            end
                        end
                    end
                end
            end
            #update ata_model.forced0
            ata_model.settings.forced0 = x_forced0_new
            DelimitedFiles.writedlm("opt/x_forced0.txt", x_forced0_new)
            message[2] =
                message[2] * "- Categorical and linear quantitative constraints grouped by friend sets.\n"
        end
        #item use
        if ata_model.settings.to_apply[1]
            item_use_max = ata_model.settings.iu.max
            item_use_max_new = zeros(Int, n_fs)
            for fs = 1:n_fs
                item_use_max_new[fs] =
                    Int(minimum(ata_model.settings.iu.max[ata_model.settings.fs.items[fs]]))
            end
            ata_model.settings.iu.max = item_use_max_new
            open("opt/Settings.jl", "a") do f
                write(f, "item_use_max = $item_use_max_new\n\n")
            end
            message[2] = message[2] * "- Maximum item use grouped.\n"
        end
        if ata_model.settings.to_apply[2]
            item_use_min = ata_model.settings.iu.min
            item_use_min_new = zeros(Int, n_fs)
            for fs = 1:n_fs
                item_use_min_new[fs] =
                    Int(maximum(ata_model.settings.iu.min[ata_model.settings.fs.items[fs]]))
            end
            ata_model.settings.iu.min = item_use_min_new
            open("opt/Settings.jl", "a") do f
                write(f, "item_use_min = $item_use_min_new\n\n")
            end
            message[2] = message[2] * "- Minimum item use grouped.\n"
        end
        message[1] = "success"
        push!(message)
    catch e
        message[1] = "danger"
        message[2] = message[2] * string("- ", sprint(showerror, e), "\n")
        push!(ata_model.output.infos, message)
    end
    return nothing
end
