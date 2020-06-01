using NCDatasets
using Dates
using CFTime
using Random

"Extract a NetCDF variable at a given time"
function get_at(var::NCDatasets.CFVariable, times::AbstractVector{<:TimeType}, t::TimeType)
    dim = findfirst(==("time"), NCDatasets.dimnames(var))
    i = findfirst(>=(t), times)
    i === nothing && throw(DomainError("time $t after dataset end $(last(times))"))
    return selectdim(var, dim, i)
end

"Get dynamic NetCDF input for the given time"
function forcing(dataset::NCDataset, time::DateTime)
    p_ = collect(get_at(dataset["P"], nctimes, starttime))
    temp_ = collect(get_at(dataset["TEMP"], nctimes, starttime))
    pet_ = collect(get_at(dataset["PET"], nctimes, starttime))

    # return only the active cells
    # TODO, the active cells need to be passed as input, not be derived independently here
    p = Float64.(filter(!ismissing, p_))
    temp = Float64.(filter(!ismissing, temp_))
    pet = Float64.(filter(!ismissing, pet_))
    @assert length(p) == length(temp) == length(pet)

    return (p = p, temp = temp, pet = pet)
end

"prepare an output dataset"
function setup_netcdf(output_path, nclon, nclat)
    ds = NCDataset(output_path, "c")
    defDim(ds, "time", Inf)  # unlimited
    defVar(
        ds,
        "lon",
        nclon,
        ("lon",),
        attrib = [
            "_FillValue" => NaN,
            "long_name" => "longitude",
            "units" => "degrees_east",
        ],
    )
    defVar(
        ds,
        "lat",
        nclat,
        ("lat",),
        attrib = [
            "_FillValue" => NaN,
            "long_name" => "latitude",
            "units" => "degrees_north",
        ],
    )
    defVar(
        ds,
        "time",
        Float64,
        ("time",),
        attrib = [
            "units" => CFTime.DEFAULT_TIME_UNITS,
            "calendar" => "proleptic_gregorian",
        ],
    )
    defVar(
        ds,
        "q",
        Float32,
        ("lon", "lat", "time"),
        attrib = ["_FillValue" => Float32(NaN)],
    )
    return ds
end

function grow_netcdf!(ds, var::AbstractString, time, A::AbstractArray)
    # index in the time dimension we want to add
    i = length(ds["time"]) + 1
    ds["time"][i] = time
    ds[var][:, :, i] = A
    return ds
end

@testset "SBM" begin
    model = Wflow.initialize_sbm_model(staticmaps_moselle_path, leafarea_moselle_path)
    param = model.vertical[1]
    @test param isa NamedTuple
    @test isbits(param)
    @test param.tt ≈ 1.2999999523162842
end

tomlpath = joinpath(@__DIR__, "config.toml")
tomldir = dirname(tomlpath)
config = Wflow.Config(TOML.parsefile(tomlpath))

# initialize a vector of SBM structs
sbm = Wflow.initialize_sbm_model(
    joinpath(tomldir, config.input.staticmaps),
    joinpath(tomldir, config.input.leafarea),
)




## reader

dataset = NCDataset(forcing_moselle_path)
nctimes = nomissing(dataset["time"][:])
nclon = Float64.(nomissing(dataset["lon"][:]))
nclat = Float64.(nomissing(dataset["lat"][:]))
dataset["P"]
dataset["PET"]
dataset["TEMP"]

# TODO check dimension ordering
dimnames(dataset["P"])

# times as defined by the simulation
starttime = DateTime(config.input.starttime)
Δt = Second(config.input.timestepsecs)
endtime = DateTime(config.input.endtime)
simtimes = starttime:Δt:endtime

p = get_at(dataset["P"], nctimes, starttime)
temp = get_at(dataset["TEMP"], nctimes, starttime)
pet = get_at(dataset["PET"], nctimes, starttime)

meteo = forcing(dataset, starttime)
@test keys(meteo) === (:p, :temp, :pet)
typeof(keys(meteo))

meteo.p
sbm

# TODO provide a mapping from other variable names in a Reader struct, which is used
# by the forcing function as well. Since it is fixed now we may as well use a NCDataset
# directly as a Reader in our model
# (p=:p, temp=:temp, pet=:pet)  # or dict
# struct Reader1
#     dataset::NCDataset
#     variables
# end

# TODO remove random string from the filename
# this makes it easier to develop for now, since we don't run into issues with open files
base, ext = splitext(config.output.path)
randomized_path = string(base, '_', randstring('a':'z', 4), ext)
output_path = joinpath(tomldir, randomized_path)
ds_out = setup_netcdf(output_path, nclon, nclat)

# q = rand(291, 313)
# grow_netcdf!(ds_out, "q", 1, q)

# close both input and output datasets
close(dataset)
close(ds_out)
rm(output_path)
