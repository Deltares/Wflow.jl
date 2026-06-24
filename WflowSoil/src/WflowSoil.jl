module WflowSoil

using StaticArrays: SVector, setindex
using FillArrays: Zeros

const MISSING_VALUE = Float64(NaN)

include("vegetation.jl")
include("soil_structs.jl")
include("soil_process.jl")
include("utils.jl")

end # WflowSoil
