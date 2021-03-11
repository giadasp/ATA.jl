
"""
add_enemies!(ata_model::AbstractModel)

# Description

Add enemy sets to the `ATA.AbstractModel` as specified in the `settings_file` loaded by [`start_ATA`](#ATA.start_ATA) function.
It requires the [`start_ATA`](#ATA.start_ATA) step.

# Arguments

- **`ata_model::AbstractModel`** : Required. The model built with `start_ATA()` and with settings loaded by [`start_ATA`](#ATA.start_ATA) function.
"""
function add_enemies!(ata_model::AbstractModel)
    bank = ata_model.settings.bank
    n_items = ata_model.settings.n_items
    if size(ata_model.constraints[1].constr_b, 1) > 0
        A = Vector{Matrix{Float64}}(undef, ata_model.settings.T)
        b = Vector{Vector{Float64}}(undef, ata_model.settings.T)
        for t = 1:ata_model.settings.T
            A[t] = ata_model.constraints[t].constr_A
            b[t] = ata_model.constraints[t].constr_b
        end

    else
        A = Vector{Matrix{Float64}}(undef, ata_model.settings.T)
        b = Vector{Vector{Float64}}(undef, ata_model.settings.T)
        for t = 1:ata_model.settings.T
            A[t] = Matrix{Float64}(undef, 0, ata_model.settings.n_items)
            b[t] = Vector{Float64}(undef, 0)
        end
    end

    if size(ata_model.settings.es.var, 1) > 0
        enemy_sets_var = ata_model.settings.es.var
        EnemySets = unique(
            vcat(
                [
                    unique(skipmissing(bank[!, (enemy_sets_var[isv])])) for
                    isv = 1:(size(enemy_sets_var, 1))
                ]...,
            ),
        )
        DelimitedFiles.writedlm("EnemySets.csv", EnemySets)
        nes = size(EnemySets, 1)
        sets = Vector{Vector{Int64}}(undef, nes)
        for es = 1:nes
            sets[es] = Int64[]
            for t = 1:ata_model.settings.T
                A[t] = vcat(A[t], zeros(ata_model.settings.n_items)')
            end
            for i = 1:n_items
                if any(
                    skipmissing(
                        vcat(
                            [
                                (
                                    (bank[!, (enemy_sets_var[isv])][i] == EnemySets[es]) ==
                                    true
                                ) for isv = 1:(size(enemy_sets_var, 1))
                            ]...,
                        ),
                    ),
                )
                    #esMaata_model.settings.Tg[n_items, es] = 1
                    for t = 1:ata_model.settings.T
                        A[t][end, i] = 1
                    end
                    sets[es] = vcat(sets[es], i)
                end
            end
            for t = 1:ata_model.settings.T
                push!(b[t], 1)
            end
        end
        for t = 1:ata_model.settings.T
            DelimitedFiles.writedlm("OPT/A_$t.csv", A[t])
            DelimitedFiles.writedlm("OPT/b_$t.csv", b[t])
        end
        #update model
        ata_model.settings.es.names = EnemySets
        ata_model.settings.es.sets = sets
        for t = 1:ata_model.settings.T
            ata_model.constraints[t].constr_b = b[t]
            ata_model.constraints[t].constr_A = A[t]
        end
        JLD2.@save "OPT/ata_model.jld2" AbstractModel
        push!(ata_model.output.infos, ["success", string("- ", nes, " enemy sets added. ")])
    end
    return nothing
end
