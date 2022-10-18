using Wflow

mutable struct ModelHandler
    model::Wflow.Model
    ModelHandler() = new()
end

handler = ModelHandler()

 function initialize(config_file::AbstractString)
    try
        handler.model = Wflow.BMI.initialize(Wflow.Model, config_file)
    catch error
        ex = ModelException()
        ex.message = string(error)
        throw(ex)
    end
end

function update()
    try
        handler.model = Wflow.BMI.update(handler.model)
        println(handler.model)
    catch error
        ex = ModelException()
        ex.message = string(error)
        throw(ex)
    end
end