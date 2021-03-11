
"""
add_enemies!(ATAmodel::AbstractModel)

# Description

Add enemy sets to the `ATA.AbstractModel` as specified in the `settings_file` loaded by [`start_ATA`](#ATA.start_ATA) function.
It requires the [`start_ATA`](#ATA.start_ATA) step.

# Arguments

- **`ATAmodel::AbstractModel`** : Required. The model built with `start_ATA()` and with settings loaded by [`start_ATA`](#ATA.start_ATA) function.
"""
function add_enemies!(ATAmodel::AbstractModel)
    bank = ATAmodel.settings.bank
    n_items = ATAmodel.settings.n_items
    if size(ATAmodel.constraints[1].constr_b, 1) > 0
        A = Vector{Matrix{Float64}}(undef, ATAmodel.settings.T)
        b = Vector{Vector{Float64}}(undef, ATAmodel.settings.T)
        for t = 1:ATAmodel.settings.T
            A[t] = ATAmodel.constraints[t].constr_A
            b[t] = ATAmodel.constraints[t].constr_b
        end

    else
        A = Vector{Matrix{Float64}}(undef, ATAmodel.settings.T)
        b = Vector{Vector{Float64}}(undef, ATAmodel.settings.T)
        for t = 1:ATAmodel.settings.T
            A[t] = Matrix{Float64}(undef, 0, ATAmodel.settings.n_items)
            b[t] = Vector{Float64}(undef, 0)
        end
    end

    if size(ATAmodel.settings.es.var, 1) > 0
        enemy_sets_var = ATAmodel.settings.es.var
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
            for t = 1:ATAmodel.settings.T
                A[t] = vcat(A[t], zeros(ATAmodel.settings.n_items)')
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
                    #esMaATAmodel.settings.Tg[n_items, es] = 1
                    for t = 1:ATAmodel.settings.T
                        A[t][end, i] = 1
                    end
                    sets[es] = vcat(sets[es], i)
                end
            end
            for t = 1:ATAmodel.settings.T
                push!(b[t], 1)
            end
        end
        for t = 1:ATAmodel.settings.T
            DelimitedFiles.writedlm("OPT/A_$t.csv", A[t])
            DelimitedFiles.writedlm("OPT/b_$t.csv", b[t])
        end
        #update model
        ATAmodel.settings.es.names = EnemySets
        ATAmodel.settings.es.sets = sets
        for t = 1:ATAmodel.settings.T
            ATAmodel.constraints[t].constr_b = b[t]
            ATAmodel.constraints[t].constr_A = A[t]
        end
        JLD2.@save "OPT/ATAmodel.jld2" AbstractModel
        push!(ATAmodel.output.infos, ["success", string("- ", nes, " enemy sets added. ")])
    end
    return nothing
end
