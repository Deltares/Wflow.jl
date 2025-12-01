using Downloads

const source_url = "https://github.com/visr/wflow-artifacts/releases/download"
const version = v"1.0.0"
const testdir = normpath(@__DIR__, "../Wflow/test")

# ensure test data is present
const datadir = joinpath(testdir, "data")
const inputdir = joinpath(datadir, "input")
isdir(inputdir) || mkpath(inputdir)

"Download a test data file if it does not already exist"
function testdata(filename)
    target_path = joinpath(inputdir, filename)
    url = string(source_url, '/', string('v', version), '/', filename)
    if isfile(target_path)
        println("- $filename already exists")
    else
        Downloads.download(url, target_path)
        println("- $filename downloaded")
    end
    return target_path
end

println("Downloading testdata from $source_url to $inputdir...")

for input in [
    "staticmaps-rhine.nc",
    "staticmaps-moselle.nc",
    "forcing-moselle.nc",
    "forcing-moselle-sed.nc",
    "staticmaps-moselle-sed.nc",
    "instates-moselle-sed.nc",
    "instates-moselle.nc",
    "forcing-sbm-groundwater-part1.nc",
    "forcing-sbm-groundwater-part2.nc",
    "staticmaps-sbm-groundwater.nc",
    "instates-example-sbm-gwf.nc",
    "reservoir_sh_1.csv",
    "reservoir_sh_2.csv",
    "reservoir_hq_2.csv",
    "forcing-calendar-noleap.nc",
    "forcing-piave.nc",
    "staticmaps-piave.nc",
    "instates-piave.nc",
    "instates-piave-gwf.nc",
]
    testdata(input)
end
