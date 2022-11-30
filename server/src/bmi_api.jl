
struct Initialize
    config_file::String
    fn::String
end
StructTypes.StructType(::Type{Initialize}) = StructTypes.Struct()

function wflow_bmi(m::Initialize)
    getfield(Wflow.BMI, Symbol(m.fn))(Wflow.Model, m.config_file)
end
