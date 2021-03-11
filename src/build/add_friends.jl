
"""
add_friends!(ATAmodel::AbstractModel)

# Description

Add friend sets to the `ATA.AbstractModel` as specified in the `settings_file` loaded by [`start_ATA`](#ATA.start_ATA) function.
It requires the [`start_ATA`](#ATA.startATA) step.

# Arguments

- **`ATAmodel::AbstractModel`** : Required. The model built with `start_ATA()` function.
"""
function add_friends!(ATAmodel::AbstractModel)
    message = ["", ""]
    if !isfile("OPT/Settings.jl")
        push!(ATAmodel.output.infos, ["danger", "Run start_ATA before!\n"])
        return nothing
    else
        n_items = ATAmodel.settings.n_items
        FriendSets =
            string.(
                unique(
                    vcat(
                        [
                            unique(
                                skipmissing(
                                    ATAmodel.settings.bank[
                                        !,
                                        (ATAmodel.settings.fs.var[isv]),
                                    ],
                                ),
                            ) for isv = 1:(size(ATAmodel.settings.fs.var, 1))
                        ]...,
                    ),
                ),
            )
        n_fs = size(FriendSets, 1)
        fs_counts = zeros(Int, n_fs)
        fs_items = [zeros(Int, 0) for i = 1:n_fs]
        single_items = Vector{Union{Missing,String}}([missing for i = 1:n_items])
        for i = 1:n_items
            units = [
                ATAmodel.settings.bank[i, ATAmodel.settings.fs.var[isv]] for
                isv = 1:(size(ATAmodel.settings.fs.var, 1))
            ]
            if all(ismissing.(units))
                single_items[i] = string(i)
            else
                nonmissing = units[.!ismissing.(units)]
                fs = Int64[]
                f_i = 1
                for f in FriendSets
                    if f in nonmissing
                        push!(fs, f_i)
                    end
                    f_i += 1
                end
                for f in fs
                    fs_items[f] = vcat(fs_items[f], i)
                    fs_counts[f] += 1
                end
            end
        end
        push!(ATAmodel.settings.fs.var, :SINGLE_fs)
        DataFrames.DataFrames.insertcols!(
            ATAmodel.settings.bank,
            size(ATAmodel.settings.bank, 2),
            :SINGLE_fs => single_items,
        )
        items_single = findall(.!ismissing.(ATAmodel.settings.bank[!, :SINGLE_fs]))
        fs_items = vcat(fs_items, [[i] for i in items_single])
        FriendSets = vcat(FriendSets, ATAmodel.settings.bank[items_single, :SINGLE_fs])
        fs_counts = vcat(fs_counts, ones(Int64, size(items_single, 1)))
        n_fs = size(fs_counts, 1)
        DelimitedFiles.writedlm("OPT/FriendSets.csv", FriendSets)

        #update model
        ATAmodel.settings.n_fs = n_fs
        ATAmodel.settings.fs.sets = FriendSets
        ATAmodel.settings.fs.items = fs_items
        ATAmodel.settings.fs.counts = fs_counts
        JLD2.@save "OPT/ATAmodel.jld2" AbstractModel
        push!(
            ATAmodel.output.infos,
            ["success", string("- ", n_fs, " friend sets added.\n")],
        )
    end
    return nothing
end
