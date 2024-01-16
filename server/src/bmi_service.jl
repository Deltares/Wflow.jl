
struct Initialize
    config_file::String
    fn::String
end

struct GetComponentName
    fn::String
end

struct GetInputItemCount
    fn::String
end

struct GetOutputItemCount
    fn::String
end

struct GetStartTime
    fn::String
end

struct GetEndTime
    fn::String
end

struct GetTimeStep
    fn::String
end

struct GetTimeUnits
    fn::String
end

struct GetCurrentTime
    fn::String
end

struct Update
    fn::String
end

struct UpdateUntil
    fn::String
    time::Float64
end

struct GetInputVarNames
    fn::String
end

struct GetOutputVarNames
    fn::String
end

struct GetVarItemSize
    fn::String
    name::String
end

struct GetVarType
    fn::String
    name::String
end

struct GetVarUnits
    fn::String
    name::String
end

struct GetVarNbytes
    fn::String
    name::String
end

struct GetVarLocation
    fn::String
    name::String
end

struct GetValue
    fn::String
    name::String
end

struct GetValuePtr
    fn::String
    name::String
end

struct GetValueAtIndices
    fn::String
    inds::Vector{Int}
    name::String
end

struct SetValue
    fn::String
    name::String
    src::Vector{Wflow.Float}
end

struct SetValueAtIndices
    fn::String
    inds::Vector{Int}
    name::String
    src::Vector{Wflow.Float}
end

struct GetGridType
    fn::String
    grid::Int
end

struct GetGridShape
    fn::String
    grid::Int
end

struct GetGridRank
    fn::String
    grid::Int
end

struct GetGridSize
    fn::String
    grid::Int
end

struct GetGridX
    fn::String
    grid::Int
end

struct GetGridY
    fn::String
    grid::Int
end

struct GetGridNodeCount
    fn::String
    grid::Int
end

struct GetVarGrid
    fn::String
    name::String
end

struct Finalize
    fn::String
end

function wflow_bmi(m::Initialize, model::Union{Wflow.Model,Nothing})
    model = getfield(Wflow.BMI, Symbol(m.fn))(Wflow.Model, m.config_file)
    return model
end

function wflow_bmi(m::GetComponentName, model::Wflow.Model)
    component_name = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("component_name" => component_name)
end

function wflow_bmi(m::GetInputItemCount, model::Wflow.Model)
    input_item_count = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("input_item_count" => input_item_count)
end

function wflow_bmi(m::GetOutputItemCount, model::Wflow.Model)
    output_item_count = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("output_item_count" => output_item_count)
end

function wflow_bmi(m::GetStartTime, model::Wflow.Model)
    start_time = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("start_time" => start_time)
end

function wflow_bmi(m::GetEndTime, model::Wflow.Model)
    end_time = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("end_time" => end_time)
end

function wflow_bmi(m::GetTimeStep, model::Wflow.Model)
    time_step = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("time_step" => time_step)
end

function wflow_bmi(m::GetTimeUnits, model::Wflow.Model)
    time_units = getfield(Wflow.BMI, Symbol(m.fn))(Wflow.Model)
    return Dict("time_units" => time_units)
end

function wflow_bmi(m::GetCurrentTime, model::Wflow.Model)
    current_time = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("current_time" => current_time)
end

function wflow_bmi(m::UpdateUntil, model::Wflow.Model)
    model = getfield(Wflow.BMI, Symbol(m.fn))(model, m.time)
    return model
end

function wflow_bmi(m::Update, model::Wflow.Model)
    model = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return model
end

function wflow_bmi(m::GetInputVarNames, model::Wflow.Model)
    input_var_names = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("input_var_names" => input_var_names)
end

function wflow_bmi(m::GetOutputVarNames, model::Wflow.Model)
    output_var_names = getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("output_var_names" => output_var_names)
end

function wflow_bmi(m::GetVarItemSize, model::Wflow.Model)
    var_itemsize = getfield(Wflow.BMI, Symbol(m.fn))(model, m.name)
    return Dict("var_itemsize" => var_itemsize)
end

function wflow_bmi(m::GetVarType, model::Wflow.Model)
    var_type = getfield(Wflow.BMI, Symbol(m.fn))(model, m.name)
    return Dict("var_type" => var_type)
end

function wflow_bmi(m::GetVarUnits, model::Wflow.Model)
    var_units = getfield(Wflow.BMI, Symbol(m.fn))(model, m.name)
    return Dict("var_units" => var_units)
end

function wflow_bmi(m::GetVarNbytes, model::Wflow.Model)
    var_nbytes = getfield(Wflow.BMI, Symbol(m.fn))(model, m.name)
    return Dict("var_nbytes" => var_nbytes)
end

function wflow_bmi(m::GetVarLocation, model::Wflow.Model)
    var_location = getfield(Wflow.BMI, Symbol(m.fn))(model, m.name)
    return Dict("var_location" => var_location)
end

function wflow_bmi(m::GetValue, model::Wflow.Model)
    value = getfield(Wflow.BMI, Symbol(m.fn))(model, m.name)
    return Dict("value" => value)
end

function wflow_bmi(m::GetValuePtr, model::Wflow.Model)
    value_ptr = getfield(Wflow.BMI, Symbol(m.fn))(model, m.name)
    return Dict("value_ptr" => value_ptr)
end

function wflow_bmi(m::GetValueAtIndices, model::Wflow.Model)
    value_at_indices = getfield(Wflow.BMI, Symbol(m.fn))(model, m.name, m.inds)
    return Dict("value_at_indices" => value_at_indices)
end

function wflow_bmi(m::SetValue, model::Wflow.Model)
    getfield(Wflow.BMI, Symbol(m.fn))(model, m.name, m.src)
    return nothing
end

function wflow_bmi(m::SetValueAtIndices, model::Wflow.Model)
    getfield(Wflow.BMI, Symbol(m.fn))(model, m.name, m.inds, m.src)
    return nothing
end

function wflow_bmi(m::GetGridType, model::Wflow.Model)
    grid_type = getfield(Wflow.BMI, Symbol(m.fn))(model, m.grid)
    return Dict("grid_type" => grid_type)
end

function wflow_bmi(m::GetGridShape, model::Wflow.Model)
    grid_shape = getfield(Wflow.BMI, Symbol(m.fn))(model, m.grid)
    return Dict("grid_shape" => grid_shape)
end

function wflow_bmi(m::GetGridRank, model::Wflow.Model)
    grid_rank = getfield(Wflow.BMI, Symbol(m.fn))(model, m.grid)
    return Dict("grid_rank" => grid_rank)
end

function wflow_bmi(m::GetGridSize, model::Wflow.Model)
    grid_size = getfield(Wflow.BMI, Symbol(m.fn))(model, m.grid)
    return Dict("grid_size" => grid_size)
end

function wflow_bmi(m::GetGridX, model::Wflow.Model)
    grid_x = getfield(Wflow.BMI, Symbol(m.fn))(model, m.grid)
    return Dict("grid_x" => grid_x)
end

function wflow_bmi(m::GetGridY, model::Wflow.Model)
    grid_y = getfield(Wflow.BMI, Symbol(m.fn))(model, m.grid)
    return Dict("grid_y" => grid_y)
end

function wflow_bmi(m::GetGridNodeCount, model::Wflow.Model)
    grid_node_count = getfield(Wflow.BMI, Symbol(m.fn))(model, m.grid)
    return Dict("grid_node_count" => grid_node_count)
end

function wflow_bmi(m::GetVarGrid, model::Wflow.Model)
    var_grid = getfield(Wflow.BMI, Symbol(m.fn))(model, m.name)
    return Dict("var_grid" => var_grid)
end

function wflow_bmi(m::Finalize, model::Wflow.Model)
    getfield(Wflow.BMI, Symbol(m.fn))(model)
    return Dict("status" => "OK")
end