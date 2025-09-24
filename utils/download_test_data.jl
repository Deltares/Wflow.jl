using Downloads

const source_url = "https://github.com/visr/wflow-artifacts/releases/download"
const testdir = normpath(@__DIR__, "../Wflow/test")

# ensure test data is present
# this code is copied from runtests.jl, and is a temporary solution to get the data in place
const datadir = joinpath(testdir, "data")
const inputdir = joinpath(datadir, "input")
isdir(inputdir) || mkpath(inputdir)

"Download a test data file if it does not already exist"
function testdata(version, source_filename, target_filename)
    target_path = joinpath(inputdir, target_filename)
    url = string(source_url, '/', string('v', version), '/', source_filename)
    if isfile(target_path)
        println("- $target_filename already exists")
    else
        Downloads.download(url, target_path)
        println("- $target_filename downloaded from $source_filename")
    end
    return target_path
end

const data_paths = Dict{String, String}()

println("Downloading testdata from $source_url to $inputdir...")

for input in [
    ("staticmaps_rhine_path", v"0.1", "staticmaps.nc", "staticmaps-rhine.nc"),
    ("staticmaps_moselle_path", v"0.3.1", "staticmaps-moselle.nc", "staticmaps-moselle.nc"),
    ("forcing_moselle_path", v"0.2.6", "forcing-moselle.nc", "forcing-moselle.nc"),
    (
        "forcing_moselle_sed_path",
        v"0.2.3",
        "forcing-moselle-sed.nc",
        "forcing-moselle-sed.nc",
    ),
    (
        "staticmaps_moselle_sed_path",
        v"0.3.0",
        "staticmaps-moselle-sed.nc",
        "staticmaps-moselle-sed.nc",
    ),
    (
        "instates_moselle_sed_path",
        v"0.3.0",
        "instates-moselle-sed.nc",
        "instates-moselle-sed.nc",
    ),
    ("instates_moselle_path", v"0.3.1", "instates-moselle.nc", "instates-moselle.nc"),
    (
        "forcing_sbm_gw_path_part1",
        v"0.2.1",
        "forcing-sbm-groundwater-part1.nc",
        "forcing-sbm-groundwater-part1.nc",
    ),
    (
        "forcing_sbm_gw_path_part2",
        v"0.2.1",
        "forcing-sbm-groundwater-part2.nc",
        "forcing-sbm-groundwater-part2.nc",
    ),
    (
        "staticmaps_sbm_gw_path",
        v"0.2.3",
        "staticmaps-sbm-groundwater.nc",
        "staticmaps-sbm-groundwater.nc",
    ),
    (
        "instates_sbm_gw_path",
        v"0.2.3",
        "instates-example-sbm-gwf.nc",
        "instates-example-sbm-gwf.nc",
    ),
    ("lake_sh_1_path", v"0.2.1", "lake_sh_1.csv", "reservoir_sh_1.csv"),
    ("lake_sh_2_path", v"0.2.1", "lake_sh_2.csv", "reservoir_sh_2.csv"),
    ("lake_hq_2_path", v"0.2.1", "lake_hq_2.csv", "reservoir_hq_2.csv"),
    (
        "forcing_calendar_noleap_path",
        v"0.2.8",
        "forcing-calendar-noleap.nc",
        "forcing-calendar-noleap.nc",
    ),
    ("forcing_piave_path", v"0.2.9", "inmaps-era5-2010-piave.nc", "forcing-piave.nc"),
    ("staticmaps_piave_path", v"0.3.1", "staticmaps-piave.nc", "staticmaps-piave.nc"),
    ("instates_piave_path", v"0.3.1", "instates-piave.nc", "instates-piave.nc"),
    ("instates_piave_gwf_path", v"0.3.1", "instates-piave-gwf.nc", "instates-piave-gwf.nc"),
]
    data_paths[input[1]] = testdata(input[2:end]...)
end