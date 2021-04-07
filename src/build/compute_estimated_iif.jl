function compute_estimated_iif!(ata_model)
    T = ata_model.settings.T
    n_items = ata_model.settings.n_items
    irt_parameters = ata_model.settings.irt.parameters
    irt_model = ata_model.settings.irt.model
    irt_D = ata_model.settings.irt.D
    irt_parametrization = ata_model.settings.irt.parametrization
    IIF = Vector{Array{Float64,2}}(undef, T)
    K = zeros(Int, T)
    for t = 1:T
        K[t] = size(ata_model.obj.cores[t].points, 1)
        IIF[t] = zeros(K[t], n_items)
        for k = 1:K[t]
            IIF[t][k, :] = item_info(
                irt_parameters,
                ata_model.obj.cores[t].points[k],
                model = irt_model,
                parametrization = irt_parametrization,
                D = irt_D,
            )# K[t] x I
        end
        ata_model.obj.cores[t].IIF = IIF[t]
    end
    JLD2.@save "opt/IIF.jld2" IIF
    return nothing
end