
struct Initialize
    path::String
    fn::String
end
StructTypes.StructType(::Type{Initialize}) = StructTypes.Struct()

function wflow_bmi(m::Initialize)
    try
        getfield(Wflow.BMI, Symbol(m.fn))(Wflow.Model, m.path)
    catch e
        @error "Wflow BMI $(m.fn) failed" exception = (e, catch_backtrace())
    end
end