using Wflow
using Downloads

# this assumes that the Wflow tests have already been run, so the data has been downloaded
testdir = abspath(dirname(pathof(Wflow)), "..", "test")

# ensure test data is present
# this code is copied from runtests.jl, and is a temporary solution to get the data in place
datadir = joinpath(testdir, "data")
inputdir = joinpath(datadir, "input")
isdir(inputdir) || mkpath(inputdir)

"Download a test data file if it does not already exist"
function testdata(version, source_filename, target_filename)
    target_path = joinpath(inputdir, target_filename)
    base_url = "https://github.com/visr/wflow-artifacts/releases/download"
    url = string(base_url, '/', string('v', version), '/', source_filename)
    isfile(target_path) || Downloads.download(url, target_path)
    return target_path
end

staticmaps_rhine_path = testdata(v"0.1", "staticmaps.nc", "staticmaps-rhine.nc")
staticmaps_moselle_path =
    testdata(v"0.2.8", "staticmaps-moselle.nc", "staticmaps-moselle.nc")
staticmaps_lahn_path = testdata(v"0.2.1", "staticmaps-lahn.nc", "staticmaps-lahn.nc")
staticmaps_meuse_path =
    testdata(v"0.2.8", "staticmaps_flex_meuse.nc", "staticmaps_flex_meuse.nc")
forcing_moselle_path = testdata(v"0.2.6", "forcing-moselle.nc", "forcing-moselle.nc")
forcing_lahn_path = testdata(v"0.2", "forcing-lahn.nc", "forcing-lahn.nc")
forcing_moselle_sed_path =
    testdata(v"0.2.3", "forcing-moselle-sed.nc", "forcing-moselle-sed.nc")
staticmaps_moselle_sed_path =
    testdata(v"0.2.3", "staticmaps-moselle-sed.nc", "staticmaps-moselle-sed.nc")
instates_moselle_sed_path =
    testdata(v"0.2", "instates-moselle-sed.nc", "instates-moselle-sed.nc")
instates_moselle_path = testdata(v"0.2.6", "instates-moselle.nc", "instates-moselle.nc")
forcing_sbm_gw_path = testdata(
    v"0.2.1",
    "forcing-sbm-groundwater-part1.nc",
    "forcing-sbm-groundwater-part1.nc",
)
forcing_sbm_gw_path = testdata(
    v"0.2.1",
    "forcing-sbm-groundwater-part2.nc",
    "forcing-sbm-groundwater-part2.nc",
)
forcing_meuse_path = testdata(v"0.2.8", "forcing_meuse.nc", "forcing_meuse.nc")
staticmaps_sbm_gw_path =
    testdata(v"0.2.2", "staticmaps-sbm-groundwater.nc", "staticmaps-sbm-groundwater.nc")
lake_sh_1_path = testdata(v"0.2.1", "lake_sh_1.csv", "lake_sh_1.csv")
lake_sh_2_path = testdata(v"0.2.1", "lake_sh_2.csv", "lake_sh_2.csv")
lake_hq_2_path = testdata(v"0.2.1", "lake_hq_2.csv", "lake_hq_2.csv")
