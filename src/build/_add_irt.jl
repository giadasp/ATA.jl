function _add_irt!(ata_model::AbstractModel)
    try
        if ata_model.settings.irt.model == "1PL"
            ata_model.settings.irt.parameters = DataFrames.DataFrame(ata_model.settings.bank[!, [:b]])#nqp values in interval\r\n",
        elseif ata_model.settings.irt.model == "2PL"
            ata_model.settings.irt.parameters = DataFrames.DataFrame(ata_model.settings.bank[!, [:a, :b]]) #nqp values in interval\r\n",
        elseif ata_model.settings.irt.model == "3PL"
            ata_model.settings.irt.parameters = DataFrames.DataFrame(ata_model.settings.bank[!, [:a, :b, :c]]) #nqp values in interval\r\n",
        else
            error!(ata_model, "Only 1PL, 2PL and 3PL IRT models are allowed.")
        end
        CSV.write("opt/irt_parameters.csv", ata_model.settings.irt.parameters)
    catch e

    end
    success!(ata_model, "IRT item parameters loaded.")
    return nothing
end
