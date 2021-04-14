function success!(ata_model::AbstractModel, message::String)
    push!(ata_model.output.infos, ["success", string("- ", message, "\n")])
    return nothing
end

function error!(ata_model::AbstractModel, message::String)
    push!(ata_model.output.infos, ["danger", string("- ", message, "\n")])
    return nothing
end
