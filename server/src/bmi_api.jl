
struct Initialize
    config_file::String
    fn::String
end

struct GetStartTime
    fn::String
end

StructTypes.StructType(::Type{Initialize}) = StructTypes.Struct()

function wflow_bmi(m::Initialize, model::Union{Wflow.Model,Nothing})
    getfield(Wflow.BMI, Symbol(m.fn))(Wflow.Model, m.config_file)
end

function wflow_bmi(m::GetStartTime, model::Wflow.Model)
    start_time = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return(Dict("start_time" => start_time))
end
