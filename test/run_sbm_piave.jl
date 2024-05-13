
function run_piave(model, steps)
    q = zeros(n)
    ssf_vol = zeros(n)
    riv_vol = zeros(n)
    for i = 1:steps
        model = Wflow.run_timestep(model)
        ssf_vol[i] = mean(model.lateral.subsurface.volume)
        riv_vol[i] = mean(model.lateral.river.volume)
        q[i] = model.lateral.river.q_av[1]
    end
    return q, riv_vol, ssf_vol
end

tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.initialize_sbm_model(config)
q_demand, riv_vol_demand, ssf_vol_demand = run_piave(model, 30)
Wflow.close_files(model, delete_output = false)

tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.initialize_sbm_model(config)
q_, riv_vol, ssf_vol = run_piave(model, 30)
Wflow.close_files(model, delete_output = false)

@testset "piave with and without water demand" begin
    idx = 1:3:28
    @test q_demand[idx] ≈ [
        211.99857232972093f0,
        187.36607290941348f0,
        256.93915861230613f0,
        158.05552787697064f0,
        158.7873979875523f0,
        140.32256192256094f0,
        116.6701746161749f0,
        100.15204336618213f0,
        162.61886168901148f0,
        97.19096285459914f0,
    ]
    @test q_[idx] ≈ [
        214.11044869595924f0,
        194.62341296307127f0,
        262.8784653320827f0,
        165.0866387234042f0,
        168.35443098223672f0,
        149.6184845898368f0,
        125.9357923707458f0,
        108.93173996691478f0,
        168.72826201444647f0,
        104.30259811567151f0,
    ]
    @test riv_vol_demand[idx] ≈ [
        59533.69348978292f0,
        54620.38570310184f0,
        61322.43260348127f0,
        47995.855664860675f0,
        50515.64469041329f0,
        42348.99287033218f0,
        40331.780886294684f0,
        38081.953202775796f0,
        42249.57246495851f0,
        36981.497741891144f0,
    ]
    @test riv_vol[idx] ≈ [
        60170.899828551424f0,
        55441.96744406564f0,
        61876.88671129906f0,
        48934.40419781716f0,
        51529.963663951305f0,
        43536.64514794115f0,
        41396.32544232497f0,
        39318.78688860317f0,
        43039.982293072826f0,
        37921.92492566406f0,
    ]
    @test ssf_vol_demand[idx] ≈ [
        242578.74578032276f0,
        237097.96619152007f0,
        235395.91298919305f0,
        230407.79827618264f0,
        226702.3018279473f0,
        221608.53938398988f0,
        216901.31104175322f0,
        212403.77029684238f0,
        211237.6730674815f0,
        207698.64010900975f0,
    ]
    @test ssf_vol[idx] ≈ [
        242540.8197370262f0,
        237053.14368167287f0,
        235345.29805380543f0,
        230352.50258261946f0,
        226654.8064348857f0,
        221565.5289960294f0,
        216859.27766727592f0,
        212365.27882439757f0,
        211055.71212857193f0,
        207513.09867755277f0,
    ]
end
