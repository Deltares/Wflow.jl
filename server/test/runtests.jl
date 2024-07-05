using ZMQ: ZMQ
using JSON3: JSON3
using StructTypes: StructTypes
using Wflow: Wflow
using WflowServer: WflowServer
import Statistics: mean
import Logging: with_logger, NullLogger
import Test: @testset, @test
using Downloads: Downloads

# ensure test data is present
testdir = @__DIR__
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

staticmaps_moselle_path = testdata(
    v"0.2.9", "staticmaps-moselle.nc", "staticmaps-moselle.nc"
)
forcing_moselle_path = testdata(v"0.2.6", "forcing-moselle.nc", "forcing-moselle.nc")
instates_moselle_path = testdata(v"0.2.6", "instates-moselle.nc", "instates-moselle.nc")

with_logger(NullLogger()) do
    @testset "Test client server Wflow ZMQ Server" begin
        include("client.jl")
    end
end
