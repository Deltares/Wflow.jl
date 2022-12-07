
struct Initialize
    config_file::String
    fn::String
end

struct GetStartTime
    fn::String
end

struct GetEndTime
    fn::String
end

struct UpdateUntil
    fn::String
    time::Float64
end

struct Finalize
    fn::String
end

function wflow_bmi(m::Initialize, model::Union{Wflow.Model,Nothing})
    model = getfield(Wflow.BMI, Symbol(m.fn))(Wflow.Model, m.config_file)
    return model
end

function wflow_bmi(m::GetStartTime, model::Wflow.Model)
    start_time = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("start_time" => start_time)
end

function wflow_bmi(m::GetEndTime, model::Wflow.Model)
    end_time = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("end_time" => end_time)
end

function wflow_bmi(m::UpdateUntil, model::Wflow.Model)
    model = getfield(Wflow.BMI, Symbol(m.fn))(model, m.time)
    return model
end

function wflow_bmi(m::Finalize, model::Wflow.Model)
    getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("status" => "OK")
end