using Wflow

mutable struct ModelHandler
    model::Union{Nothing, Wflow.Model}
end

function initialize(config_file::AbstractString)
    try
        handler.model = Wflow.BMI.initialize(Wflow.Model, config_file)
    catch error
        println(error)
        ex = ModelException()
        ex.message = "testing"
        throw(ex)
    end
end

function update()
    try
        handler.model = Wflow.BMI.update()
    catch error
        ex = ModelException()
        ex.message = string(error)
        throw(ex)
    end
end

handler = ModelHandler(nothing)