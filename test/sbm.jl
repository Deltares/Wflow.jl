@testset "SBM" begin
    model = Wflow.initialize_sbm_model(staticmaps_moselle_path, leafarea_moselle_path)
    param = model.vertical[1]
    @test param isa Wflow.SBM{Float64,4,5}
    @test isbits(param)
    @test param.tt ≈ 1.2999999523162842
end
