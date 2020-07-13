@testset "configuration file" begin
    tomlpath = joinpath(@__DIR__, "config.toml")
    parsed_toml = Wflow.parsefile(tomlpath)
    config = Wflow.Config(parsed_toml)
    @test parsed_toml isa Dict{String,Any}
    @test config isa Wflow.Config
    @test getfield(config, :dict) === parsed_toml

    # test if the values are parsed as expected
    @test config.casename == "testcase"
    @test config.λ == 1.2
    @test config.input.starttime === DateTime(2000)
    @test config.input.endtime === DateTime(2000, 2)
    @test config.output.path == "data/specified_output.nc"
    @test config.output.parameters isa Vector
    @test config.output.parameters == [
        "satwaterdepth",
        "snow",
        "tsoil",
        "ustorelayerdepth",
        "snowwater",
        "canopystorage",
    ]
end
