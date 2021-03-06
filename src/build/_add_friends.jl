
"""
    _add_friends!(ata_model::AbstractModel)
"""
function _add_friends!(ata_model::AbstractModel)

    n_items = ata_model.settings.n_items
    FriendSets =
        string.(
            unique(
                vcat(
                    [
                        unique(
                            skipmissing(
                                ata_model.settings.bank[
                                    !,
                                    (ata_model.settings.fs.var[isv]),
                                ],
                            ),
                        ) for isv = 1:(size(ata_model.settings.fs.var, 1))
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
            ata_model.settings.bank[i, ata_model.settings.fs.var[isv]] for
            isv = 1:(size(ata_model.settings.fs.var, 1))
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
    push!(ata_model.settings.fs.var, :SINGLE_fs)
    DataFrames.insertcols!(
        ata_model.settings.bank,
        size(ata_model.settings.bank, 2),
        :SINGLE_fs => single_items,
    )
    items_single = findall(.!ismissing.(ata_model.settings.bank[!, :SINGLE_fs]))
    fs_items = vcat(fs_items, [[i] for i in items_single])
    FriendSets = vcat(FriendSets, ata_model.settings.bank[items_single, :SINGLE_fs])
    fs_counts = vcat(fs_counts, ones(Int64, size(items_single, 1)))
    n_fs = size(fs_counts, 1)
    DelimitedFiles.writedlm("opt/FriendSets.csv", FriendSets)

    #update model
    ata_model.settings.n_fs = n_fs
    ata_model.settings.fs.sets = FriendSets
    ata_model.settings.fs.items = fs_items
    ata_model.settings.fs.counts = fs_counts
    JLD2.@save "opt/ata_model.jld2" AbstractModel
    success!(ata_model, string(n_fs, " friend sets added."))

    return nothing
end
