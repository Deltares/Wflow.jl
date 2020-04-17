@testset "SBM" begin
    sbm = Wflow.initialize(staticmaps_moselle_path, leafarea_moselle_path)
    param = sbm[1]
    @test param isa Wflow.SBM{4,5}
    @test isbits(param)
    @test param.tt ≈ 1.2999999523162842
end
