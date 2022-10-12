using Wflow

mutable struct ModelHandler
    model::Union{Nothing, Wflow.Model}
end

function initialize(config_file::AbstractString)
    handler.model = Wflow.BMI.initialize(Wflow.Model, config_file)
end

function update()
    handler.model = Wflow.BMI.update()
end

handler = ModelHandler(nothing)