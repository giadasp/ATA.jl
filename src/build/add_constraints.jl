

"""
add_constraints!(ata_model::AbstractModel; constraints_file = "constraints.csv", constraints_delim = ";")

# Description

Add categorical and sum constraints to the `ATA.AbstractModel` as specified in the `constraints_file`.
It requires the [`start_ATA`](#ATA.start_ATA) step.

# Arguments

- **`ata_model::AbstractModel`** : Required. The model built with `start_ATA()` and with settings loaded by [`start_ATA`](#ATA.start_ATA) function.
- **`constraints_file`** : Optional. Default: "constraints.csv". The path of the file containing the categorical and sum constraints in the form of custom-separated values.
- **`constraints_delim`** : Optional. Default: ";". The custom-separator for the `constraints_file`.

"""
function add_constraints!(
    ata_model::AbstractModel;
    constraints_file = "constraints.csv",
    constraints_delim = ";",
)
    n_items = ata_model.settings.n_items
    bank = ata_model.settings.bank
    A = Vector{Matrix{Float64}}(undef, ata_model.settings.T)
    b = Vector{Vector{Float64}}(undef, ata_model.settings.T)
    if size(ata_model.constraints[1].constr_b, 1) > 0
        for t = 1:ata_model.settings.T
            A[t] = copy(ata_model.constraints[t].constr_A)
            b[t] = copy(ata_model.constraints[t].constr_b)
        end
    else
        for t = 1:ata_model.settings.T
            A[t] = zeros(Float64, 0, ata_model.settings.n_items)
            b[t] = Vector{Float64}(undef, 0)
        end

    end
    if !isfile(constraints_file)
        push!(
            ata_model.output.infos,
            ["danger", string(constraints_file, " doesn't exist.")],
        )
        return nothing
    else
        message = ["", ""]
        x_forced0 = ata_model.settings.forced0
        Categoricalconsts =
            CSV.read(constraints_file, DataFrames.DataFrame, delim = constraints_delim)
        if size(Categoricalconsts, 1) == 0
            message[2] =
                message[2] *
                string("No lines in file ", constraints_file, " nothing added.\n")
        else
            Categoricalconsts.var = Symbol.(Categoricalconsts.var)
            CatCons = Vector{Vector{Float64}}(undef, ata_model.settings.T)
            for g = 1:ata_model.settings.n_groups
                if g > 1
                    tests = collect(
                        (sum(
                            ata_model.settings.Tg[1:(g-1)],
                        )+1):(sum(
                            ata_model.settings.Tg[1:(g-1)],
                        )+ata_model.settings.Tg[g]),
                    )
                else
                    tests = collect(1:ata_model.settings.Tg[1])
                end
                GroupCatCons = findall((Categoricalconsts.group .== g))
                for con in GroupCatCons
                    for t in tests
                        if !ismissing(Categoricalconsts[con, :min])
                            indices =
                                string.(bank[!, Symbol.(Categoricalconsts[con, :var])]) .==
                                string(Categoricalconsts[con, :value])
                            indices[ismissing.(indices)] .= false
                            indices = findall(indices .== true)
                            A[t] = vcat(A[t], zeros(Float64, 1, n_items))
                            if ismissing(Categoricalconsts[con, :weight])
                                A[t][end, indices] .= -1.0
                                push!(b[t], -Categoricalconsts[con, :min])
                            else
                                A[t][end, indices] .= -(Categoricalconsts[con, :weight])
                                push!(
                                    b[t],
                                    -Categoricalconsts[con, :weight] *
                                    Categoricalconsts[con, :min],
                                )
                            end
                        end
                        if !ismissing(Categoricalconsts[con, :max])
                            indices =
                                string.(
                                    ata_model.settings.bank[
                                        !,
                                        Symbol.(Categoricalconsts[con, :var]),
                                    ],
                                ) .== string(Categoricalconsts[con, :value])
                            indices[ismissing.(indices)] .= false
                            indices = findall(indices .== true)
                            if Categoricalconsts[con, :max] == 0
                                x_forced0[t][indices] .= false
                            else
                                A[t] = vcat(A[t], zeros(Float64, 1, n_items))
                                if ismissing(Categoricalconsts[con, :weight])
                                    A[t][end, indices] .= 1.0
                                    push!(b[t], Categoricalconsts[con, :max])
                                else
                                    A[t][end, indices] .= (Categoricalconsts[con, :weight])
                                    push!(
                                        b[t],
                                        Categoricalconsts[con, :weight] *
                                        Categoricalconsts[con, :max],
                                    )
                                end
                            end
                        end
                    end
                end
            end
            #add quantitative constraints
            for t = 1:ata_model.settings.T
                for var = 1:size(ata_model.constraints[t].sum_vars, 1)
                    vals = copy(
                        ata_model.settings.bank[!, ata_model.constraints[t].sum_vars[var]],
                    )
                    vals[findall(ismissing.(vals))] .= 0.0
                    if !ismissing(ata_model.constraints[t].sum_vars_min[var])
                        A[t] = vcat(A[t], .-vals')
                        push!(b[t], -ata_model.constraints[t].sum_vars_min[var])
                    end
                    if !ismissing(ata_model.constraints[t].sum_vars_max[var])
                        A[t] = vcat(A[t], vals')
                        push!(b[t], ata_model.constraints[t].sum_vars_max[var])
                    end
                end
            end
            message[2] = message[2] * "- Sum variables constraints added. \n"
            for t = 1:ata_model.settings.T
                DelimitedFiles.writedlm("OPT/A_$t.csv", A[t])
                DelimitedFiles.writedlm("OPT/b_$t.csv", b[t])
            end
            DelimitedFiles.writedlm("OPT/x_forced0.txt", x_forced0)
        end
        #update model
        ata_model.settings.forced0 = x_forced0
        for t = 1:ata_model.settings.T
            ata_model.constraints[t].constr_b = b[t]
            ata_model.constraints[t].constr_A = A[t]
        end
        JLD2.@save "OPT/ata_model.jld2" AbstractModel
        message[1] = "success"
        message[2] = message[2] * string("- Constraints in file ", constraints_file ," added. \n")
        push!(ata_model.output.infos, message)
    end
    return nothing
end
