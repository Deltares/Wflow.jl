module SharedHydrology

using StaticArrays: SVector

const MISSING_VALUE = Float64(NaN)

include("vegetation.jl")
include("soil/soil_structs.jl")
include("soil/soil_process.jl")

end # SharedHydrology