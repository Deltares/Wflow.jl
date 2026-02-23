@testitem "Run SBM" begin
    using Wflow: to_SI, ABSOLUTE_DEGREES, MM_PER_DT, M3_PER_DAY, MM, Unit
    CM = Unit(; cm = 1)
    M_PER_DT = Unit(; m = 1, dt = -1)
    M_PER_DAY = Unit(; m = 1, d = -1)
    M2_PER_DAY = Unit(; m = 2, d = -1)
    using Dates
    using Statistics: mean
    include("testing_utils.jl")

    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    dt = Float64(config.time.timestepsecs)

    model = Wflow.Model(config)
    (; domain) = model

    Wflow.run_timestep!(model)

    # test if the first timestep was written to the CSV file
    flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
    @testset "CSV output" begin
        row = csv_first_row(model.writer.csv_path)

        @test row.time == DateTime("2000-01-02T00:00:00")
        @test row.Q ≈ 6.850474762647178
        @test row.storage ≈ 2.753500258633802e7
        @test row.temp_bycoord ≈ 2.390000104904175
        @test row.vwc_layer2_bycoord ≈ 0.25938707367146907
        @test row.temp_byindex ≈ 2.390000104904175
        @test row.Q_6336050 ≈ 0.003175570694335172
        @test row.Q_6336510 ≈ 0.012785802015566708
        @test row.Q_6836100 ≈ 0.008740434542207666
        @test row.Q_6336500 ≈ 0.00240582054420714
        @test row.Q_6836190 ≈ 0.002534862969007946
        @test row.Q_6336800 ≈ 0.004367717681790628
        @test row.Q_6336900 ≈ 0.003934027688070696
        @test row.Q_6336930 ≈ 0.013810820648092513
        @test row.Q_6336910 ≈ 0.003073485680460027
        @test row.Q_6136500 ≈ 0.00043538067672607716
        @test row.Q_6136520 ≈ 0.0004119095811616616
        @test row.Q_6136150 ≈ 0.0029657574491729643
        @test row.Q_6136151 ≈ 0.0027844732219454186
        @test row.Q_6136160 ≈ 3.472843959583444
        @test row.Q_6136202 ≈ 1.306258882043624
        @test row.recharge_1 ≈ -0.002257181032501202
    end

    @testset "NetCDF scalar output" begin
        ds = model.writer.dataset_scalar
        @test ds["time"][1] == DateTime("2000-01-02T00:00:00")
        @test ds["Q"][:][1:20] ≈ [
            0.7411495,
            1.7804804,
            1.7792363,
            1.5534992,
            13.25166,
            4.7680607,
            0.09143825,
            4.354929,
            0.010488547,
            6.0434036,
            6.0332,
            0.0059943455,
            0.00660558,
            0.006399016,
            0.0041626506,
            1.2043761,
            0.0022155696,
            1.7879109,
            8.354377,
            2.735113,
        ]
        @test ds["river_gauge__count"].attrib["cf_role"] == "timeseries_id"
        @test ds["temp_index"][:] ≈ [2.39]
        @test ds["temp_coord"][:] ≈ [2.39]
        @test keys(ds.dim) ==
              ["time", "layer", "river_gauge__count", "temp_bycoord", "temp_byindex"]
    end

    @testset "First timestep: interception" begin
        (; interception) = model.land

        @test model.clock.iteration == 1

        @test test_means(interception.parameters, Dict(:e_r => 0.2381398056521972))
        @test test_means(
            interception.parameters.vegetation_parameter_set,
            Dict(
                :kc => 1.0,
                :storage_wood => to_SI(0.19583714217432435, MM),
                :canopygapfraction => 0.5076891913529586,
                :storage_specific_leaf => to_SI(0.08942700997221185, MM),
                :leaf_area_index => 1.0624380387861794,
                :cmax => to_SI(0.2888534265722341, MM),
                :kext => 0.6741606033862981,
                :rootingdepth => to_SI(370.1426593311398, MM),
            ),
        )
        @test test_means(
            interception.variables,
            Dict(
                :stemflow => to_SI(0.0962877848464674, MM_PER_DT; dt_val = dt),
                :canopy_storage => 0.0,
                :canopy_potevap => to_SI(0.3142309029479176, MM_PER_DT; dt_val = dt),
                :interception_rate => to_SI(0.3138827922906597, MM_PER_DT; dt_val = dt),
                :throughfall => to_SI(1.6122206031467872, MM_PER_DT; dt_val = dt),
            ),
        )
    end

    @testset "First timestep: snow" begin
        (; snow) = model.land

        @test snow.variables.snow_storage[5] ≈ to_SI(3.768513390588815, MM)
        @test snow.parameters.tt[50063] ≈ to_SI(0.0, ABSOLUTE_DEGREES)

        @test test_means(
            snow.parameters,
            Dict(
                :whc => 0.1,
                :cfmax => to_SI(3.75653, Unit(; mm = 1, degC = -1, dt = -1); dt_val = dt),
                :tt => to_SI(0.0, ABSOLUTE_DEGREES),
                :ttm => to_SI(0.0, ABSOLUTE_DEGREES),
                :tti => 2.0,
            ),
        )
        @test test_means(
            snow.variables,
            Dict(
                :swe => to_SI(0.04059539652492171, MM),
                :runoff => to_SI(1.6679129914683328, MM_PER_DT; dt_val = dt),
                :snow_in => to_SI(3.688075238243213e-7, MM_PER_DT; dt_val = dt),
                :snow_melt => 0.0,
                :snow_out => to_SI(3.688075238243213e-7, MM_PER_DT; dt_val = dt),
                :snow_storage => to_SI(0.03801972367609432, MM),
                :snow_water => to_SI(0.0025756728488273866, MM),
            ),
        )
    end

    @testset "First timestep: runoff" begin
        (; runoff) = model.land

        @test test_means(
            runoff.boundary_conditions,
            Dict(
                :waterdepth_river => 0.0,
                :water_flux_surface => to_SI(1.6679129914683328, MM_PER_DT; dt_val = dt),
                :waterdepth_land => 0.0,
            ),
        )
        @test test_means(
            runoff.variables,
            Dict(
                :ae_openw_r => 0.0,
                :runoff_river => to_SI(0.02192394495624879, MM_PER_DT; dt_val = dt),
                :net_runoff_river => to_SI(0.02192394495624879, MM_PER_DT; dt_val = dt),
                :ae_openw_l => 0.0,
                :runoff_land => to_SI(0.001511577667599544, MM_PER_DT; dt_val = dt),
            ),
        )
    end

    @testset "First timestep: soil" begin
        (; soil) = model.land

        @test soil.variables.total_storage[50063] ≈ to_SI(559.886196102425, MM)
        @test soil.variables.total_storage[429] ≈ to_SI(599.6375723706255, MM) # river cell
        @test soil.parameters.theta_s[50063] ≈ 0.48755401372909546
        @test soil.parameters.theta_r[50063] ≈ 0.15943120419979095

        @test test_means(
            soil.boundary_conditions,
            Dict(
                :potential_soilevaporation =>
                    to_SI(0.31617207529997243, MM_PER_DT; dt_val = dt),
                :water_flux_surface => to_SI(1.6444774688444848, MM_PER_DT; dt_val = dt),
                :potential_transpiration =>
                    to_SI(0.0003481106572579187, MM_PER_DT; dt_val = dt),
            ),
        )
        @test test_means(
            soil.parameters,
            Dict(
                :theta_s => 0.4409211971535584,
                :soilthickness => to_SI(1837.833713668117, MM),
                :h1 => 0.0,
                :cap_n => 2.0,
                :sumlayers =>
                    to_SI.(
                        SVector((
                            0.0,
                            100.0,
                            399.97241953958815,
                            1156.6640303676966,
                            1984.2496012545644,
                        )),
                        Ref(MM),
                    ),
                :c => SVector((
                    9.428788533549843,
                    9.821687673542042,
                    10.240060684266773,
                    10.24248827673959,
                )),
                :h4 => to_SI(-16000.0, CM),
                :w_soil => 0.11249999999999985,
                :hb => to_SI(-10.0, CM),
                :kvfrac => SVector((1.0, 1.0, 1.0, 1.0)),
                :h3_low => to_SI(-1000.0, CM),
                :soil_fraction => 0.49508064804427593,
                :cap_hmax => to_SI(2000.0, MM),
                :act_thickl =>
                    to_SI.(
                        SVector((
                            100.0,
                            299.97241953958815,
                            756.6640303676965,
                            784.2496012545644,
                        )),
                        Ref(MM),
                    ),
                :pathfrac => 0.0129736199199603,
                :h3_high => to_SI(-400.0, CM),
                :infiltcappath => to_SI(5.0, MM_PER_DT; dt_val = dt),
                :infiltcapsoil => to_SI(418.6617132363326, MM_PER_DT; dt_val = dt),
                :cf_soil => 0.038,
                :theta_fc => 0.30670581485821735,
                :maxleakage => 0.0,
                :h2 => to_SI(-100.0, CM),
                :rootdistpar => to_SI(-500.0, Unit(; mm = -1)),
                :alpha_h1 => 1.0,
                :soilwatercapacity => to_SI(505.7182462428663, MM),
                :rootfraction => SVector((
                    0.2769545250424278,
                    0.7118645592933596,
                    0.011101016337365556,
                    0.0,
                )),
                :theta_r => 0.16574578136033039,
            ),
        )
        @test test_means(
            soil.variables,
            Dict(
                :transfer => to_SI(1.1890024449871777e-8, MM_PER_DT; dt_val = dt),
                :drainable_waterdepth => to_SI(210.2789323014635, MM),
                :infiltexcess => 0.0,
                :actinfilt => to_SI(1.6444774688444848, MM_PER_DT; dt_val = dt),
                :ustorelayerdepth => to_SI(
                    [1.6280832346593193, 0.6460862670135344, 0.23103976228101625, 0.0],
                    MM,
                ),
                :h3 => to_SI(-1000.0, CM),
                :recharge => to_SI(-0.0022571810325012027, MM_PER_DT; dt_val = dt),
                :ustorelayerthickness => to_SI(
                    [99.3509423527741, 190.22782953886028, 89.8015331913427, 0.0],
                    MM,
                ),
                :soilevap => to_SI(0.02122698830889416, MM_PER_DT; dt_val = dt),
                :net_runoff => to_SI(0.04033701711287029, MM_PER_DT; dt_val = dt),
                :actinfiltpath => to_SI(0.02144198548416092, MM_PER_DT; dt_val = dt),
                :ustorecapacity => to_SI(74.42171691684314, MM),
                :ae_ustore => 0.0,
                :vwc_root => 0.24099999255961077,
                :runoff => to_SI(0.04033701711287029, MM_PER_DT; dt_val = dt),
                :exfiltsatwater => to_SI(0.03882543944527075, MM_PER_DT; dt_val = dt),
                :total_soilwater_storage => to_SI(431.2965293260232, MM),
                :soilevapsat => to_SI(0.0019107351725580284, MM_PER_DT; dt_val = dt),
                :ustoredepth => to_SI(2.50520926395387, MM),
                :rootstore => to_SI(29.03737411323857, MM),
                :infiltsoilpath => to_SI(1.6444774688444848, MM_PER_DT; dt_val = dt),
                :zi => to_SI(278.73782910028507, MM),
                :vwc_percroot => 54.55735722947444,
                :total_storage => to_SI(431.55733671295263, MM),
                :actleakage => 0.0,
                :actevap => to_SI(0.33545623834952143, MM_PER_DT; dt_val = dt),
                :satwaterdepth => to_SI(428.79132006206936, MM),
                :transpiration => to_SI(0.0003464577499676236, MM_PER_DT; dt_val = dt),
                :actinfiltsoil => to_SI(1.6230354833603235, MM_PER_DT; dt_val = dt),
                :actcapflux => 0.0,
                :vwc_perc =>
                    [41.86055714460927, 63.23920727815791, 99.92432179078311, 100.0],
                :tsoil => to_SI(9.243668273365222, ABSOLUTE_DEGREES),
                :f_infiltration_reduction => 1.0,
                :actevapsat => to_SI(0.0003464577499676236, MM_PER_DT; dt_val = dt),
                :vwc => [
                    0.1842331927167544,
                    0.27939237676857903,
                    0.4405326102482261,
                    0.4383915015384879,
                ],
            ),
        )
    end

    @testset "First timestep: subsurface routing" begin
        (; subsurface_flow) = model.routing

        @test test_means(
            subsurface_flow.parameters,
            Dict(
                :soilthickness => 1.8378337136681173,
                :specific_yield => 0.134215382295341,
                :khfrac => 100.0,
                :specific_yield_dyn => 0.18809013342973288,
            ),
        )
        @test test_means(
            subsurface_flow.parameters.kh_profile,
            Dict(:f => 3.303715296489842, :kh_0 => to_SI(41.86617132363326, M_PER_DAY)),
        )
        @test test_means(
            subsurface_flow.variables,
            Dict(
                :zi => 0.2787378291002851,
                :to_river => to_SI(83.98143291277904, M3_PER_DAY),
                :ssfin => to_SI(1205.6076607372684, M3_PER_DAY),
                :storage => 118195.04455118616,
                :ssfmax => to_SI(3.149179231329745, M2_PER_DAY),
                :ssf => to_SI(1289.6589215389067, M2_PER_DAY),
                :exfiltwater => to_SI(3.8825439445270756e-5, M_PER_DT; dt_val = dt),
            ),
        )
    end

    @testset "First timestep: overland routing" begin
        (; overland_flow) = model.routing

        @test test_means(
            overland_flow.boundary_conditions,
            Dict(:inwater => 0.00026861922353565927),
        )
        @test test_means(
            overland_flow.parameters,
            Dict(
                :alpha => 19.61031861512466,
                :slope => 0.09473195535994473,
                :alpha_term => 1.4127165363578293,
                :mannings_n => 0.48407901249960045,
            ),
        )
        @test test_means(overland_flow.variables, Dict(:to_river => 9.421773254130464e-5))
        @test test_means(
            overland_flow.variables.flow,
            Dict(
                :qin => 0.0006926507774505128,
                :storage => 15.068288821912116,
                :h => 2.649435687032212e-5,
                :qin_av => 0.0006878506670979908,
                :qlat => 2.637385225086365e-7,
                :q => 0.0007875271901973019,
                :q_av => 0.0007820683996392956,
            ),
        )
    end

    @testset "First timestep: river routing" begin
        (; river_flow) = model.routing

        @test test_means(
            river_flow.boundary_conditions,
            Dict(
                :actual_external_abstraction_av => 0.0,
                :external_inflow => 0.0,
                :abstraction => 0.0,
                :inwater => 0.010415800149478108,
            ),
        )
        @test test_means(river_flow.parameters, Dict(:bankfull_depth => 1.1024306883272645))
        @test test_means(
            river_flow.parameters.flow,
            Dict(
                :alpha => 7.270795158143269,
                :slope => 0.0031017516045770657,
                :alpha_term => 1.2676399832938348,
                :mannings_n => 0.030500946218238616,
            ),
        )
        @test test_means(
            river_flow.variables,
            Dict(
                :qin => 0.30373899023076173,
                :storage => 947.2681057776836,
                :h => 0.019348596909199465,
                :qin_av => 0.15655823559102536,
                :qlat => 1.808757248203136e-5,
                :q => 0.3031964891388029,
                :q_av => 0.15601028451622473,
            ),
        )
    end

    # run the second timestep
    Wflow.run_timestep!(model)

    @testset "Second timestep: snow" begin
        (; snow) = model.land

        @test snow.variables.snow_storage[5] ≈ to_SI(3.843412524052313, MM)

        @test test_means(
            snow.variables,
            Dict(
                :swe => to_SI(0.03583038054459421, MM; dt_val = dt),
                :runoff => to_SI(0.09007336794521474, MM_PER_DT; dt_val = dt),
                :snow_in => to_SI(3.7770667468749113e-7, MM_PER_DT; dt_val = dt),
                :snow_melt => to_SI(0.005351364198880822, MM_PER_DT; dt_val = dt),
                :snow_out => to_SI(3.777066746874912e-7, MM_PER_DT; dt_val = dt),
                :snow_storage => to_SI(0.03461317061870947, MM),
                :snow_water => to_SI(0.0012172099258847414, MM),
            ),
        )
    end

    @testset "Second timestep: interception" begin
        (; interception) = model.land

        @test test_means(
            interception.variables,
            Dict(
                :canopy_storage => 0.0,
                :stemflow => to_SI(0.0075251278256163005, MM_PER_DT; dt_val = dt),
                :canopy_potevap => to_SI(0.26344947521068923, MM_PER_DT; dt_val = dt),
                :interception_rate => to_SI(0.0697346938492787, MM_PER_DT; dt_val = dt),
                :throughfall => to_SI(0.07778322413927095, MM_PER_DT; dt_val = dt),
            ),
        )
    end

    @testset "Second timestep: runoff" begin
        (; runoff) = model.land

        @test test_means(
            runoff.variables,
            Dict(
                :ae_openw_r => to_SI(0.006319657683443494, MM_PER_DT; dt_val = dt),
                :runoff_river => to_SI(0.0011240257971819377, MM_PER_DT; dt_val = dt),
                :net_runoff_river => to_SI(-0.00519563188626155, MM_PER_DT; dt_val = dt),
                :ae_openw_l => to_SI(0.0003684619419575865, MM_PER_DT; dt_val = dt),
                :runoff_land => to_SI(0.00014222791218190106, MM_PER_DT; dt_val = dt),
            ),
        )
    end

    @testset "Second timestep: soil" begin
        (; soil) = model.land

        @test soil.parameters.theta_s[50063] ≈ 0.48755401372909546
        @test soil.parameters.theta_r[50063] ≈ 0.15943120419979095
        # Total storage is affected by all modules, not just soil
        @test soil.variables.total_storage[50063] ≈ to_SI(559.9633189802773, MM)
        @test soil.variables.total_storage[429] ≈ to_SI(620.6174448020174, MM)  # river cell

        @test test_means(
            soil.variables,
            Dict(
                :transfer => to_SI(0.013450017131698456, MM_PER_DT; dt_val = dt),
                :drainable_waterdepth => to_SI(209.50390590615132, MM),
                :infiltexcess => 0.0,
                :actinfilt => to_SI(0.0886192861442016, MM_PER_DT; dt_val = dt),
                :ustorelayerdepth => to_SI(
                    [
                        1.7064168532989077,
                        1.0896449436965638,
                        0.6820793025801899,
                        0.00010165257547301471,
                    ],
                    MM,
                ),
                :h3 => to_SI(-1000.0, CM),
                :recharge => to_SI(-0.16424404920993543, MM_PER_DT; dt_val = dt),
                :ustorelayerthickness => to_SI(
                    [
                        99.40922626886153,
                        192.70470155860653,
                        171.82135284713667,
                        30.507142365925574,
                    ],
                    MM,
                ),
                :soilevap => to_SI(0.01873822620551053, MM_PER_DT; dt_val = dt),
                :net_runoff => to_SI(0.2209225219909784, MM_PER_DT; dt_val = dt),
                :actinfiltpath => to_SI(0.001450556592956387, MM_PER_DT; dt_val = dt),
                :ustorecapacity => to_SI(74.89923253956516, MM),
                :ae_ustore => to_SI(0.0050466712209524854, MM_PER_DT; dt_val = dt),
                :vwc_root => 0.24143642302182922,
                :runoff => to_SI(0.221290983932936, MM_PER_DT; dt_val = dt),
                :exfiltsatwater => to_SI(0.22096092792910482, MM_PER_DT; dt_val = dt),
                :total_soilwater_storage => to_SI(430.8190137033012, MM),
                :soilevapsat => to_SI(0.001778482147882288, MM_PER_DT; dt_val = dt),
                :excesswatersoil => to_SI(0.00018782809164930448, MM_PER_DT; dt_val = dt),
                :ustoredepth => to_SI(3.4782427521511345, MM),
                :rootstore => to_SI(29.196404187600344, MM),
                :infiltsoilpath => to_SI(0.0886192861442016, MM_PER_DT; dt_val = dt),
                :zi => to_SI(283.1713993535127, MM),
                :vwc_percroot => 54.646709985548235,
                # Total storage is affected by all modules, not just soil
                :total_storage => to_SI(431.44809617518115, MM),
                :actleakage => 0.0,
                :actevap => to_SI(0.2744312454440505, MM_PER_DT; dt_val = dt),
                :satwaterdepth => to_SI(427.34077095115, MM),
                :transpiration => to_SI(0.17927020576386024, MM_PER_DT; dt_val = dt),
                :actinfiltsoil => to_SI(0.0871687295512452, MM_PER_DT; dt_val = dt),
                :actcapflux => to_SI(0.0016920496508438876, MM_PER_DT; dt_val = dt),
                :vwc_perc => [
                    42.213917397074624,
                    63.201354128983084,
                    99.7692631420567,
                    99.99996110706128,
                ],
                :tsoil => to_SI(8.55200654167359, ABSOLUTE_DEGREES),
                :excesswater => to_SI(0.00018782809164930394, MM_PER_DT; dt_val = dt),
                :f_infiltration_reduction => 1.0,
                :actevapsat => to_SI(0.17422353454290776, MM_PER_DT; dt_val = dt),
                :vwc => [
                    0.18585985815835235,
                    0.27925583454823527,
                    0.4397834998887272,
                    0.4383913077433366,
                ],
            ),
        )
    end

    @testset "Second timestep: subsurface routing" begin
        (; subsurface_flow) = model.routing
        (; ssf) = subsurface_flow.variables

        @test ssf[domain.land.network.order[1]] ≈ to_SI(699.3636285243076, M3_PER_DAY)
        @test ssf[domain.land.network.order[end - 100]] ≈
              to_SI(2395.6159482448143, M3_PER_DAY)
        @test ssf[domain.land.network.order[end]] ≈ to_SI(287.61501877867994, M3_PER_DAY)

        @test test_means(
            subsurface_flow.variables,
            Dict(
                :zi => 0.2831713993535127,
                :to_river => to_SI(83.4857164905592, M3_PER_DAY),
                :ssfin => to_SI(1164.889898167789, M3_PER_DAY),
                :storage => 117753.99735856777,
                :ssfmax => to_SI(3.149179231329745, M2_PER_DAY),
                :ssf => to_SI(1248.442951721258, M2_PER_DAY),
                :exfiltwater => to_SI(0.00022096092792910484, M_PER_DT; dt_val = dt),
            ),
        )
    end

    @testset "Second timestep: overland routing" begin
        (; overland_flow) = model.routing

        q = Wflow.get_average(overland_flow.variables.q_av)
        @test q[26625] ≈ 0.0
        @test q[39308] ≈ 0.0
        @test q[domain.land.network.order[end]] ≈ Wflow.KIN_WAVE_MIN_FLOW

        @test test_means(
            overland_flow.boundary_conditions,
            Dict(:inwater => 0.001470763758351302),
        )
        @test test_means(overland_flow.variables, Dict(:to_river => 0.0006258080185533339))
        @test test_means(
            overland_flow.variables.flow,
            Dict(
                :qin => 0.006912379626717226,
                :storage => 88.07246474045597,
                :h => 0.00015534152871148443,
                :qin_av => 0.004661418794902121,
                :qlat => 1.6550140762075134e-6,
                :q => 0.00797347889594406,
                :q_av => 0.005287226813455454,
            ),
        )
    end

    @testset "river flow" begin
        q = Wflow.get_average(model.routing.river_flow.variables.q_av)
        @test sum(q) ≈ 3696.5789619579173
        @test q[1622] ≈ 0.0007502850311515928
        @test q[43] ≈ 11.458528971675033
        @test q[domain.river.network.order[end]] ≈ 0.04397648668235761
    end

    @testset "reservoir simple" begin
        res = model.routing.river_flow.boundary_conditions.reservoir
        @test res.variables.outflow[1] ≈ 0.2174998614438593
        @test Wflow.get_average(res.variables.outflow_av)[1] ≈ 0.21749986282401396
        @test Wflow.get_average(res.boundary_conditions.inflow)[1] ≈ 0.0005130607130643568
        @test res.variables.storage[1] ≈ 2.751299001489657f7
        @test res.variables.actevap.cumulative_material[1] ≈
              to_SI(0.5400000810623169, MM; dt_val = dt)
        @test res.boundary_conditions.precipitation[1] ≈
              to_SI(0.17999997735023499, MM_PER_DT; dt_val = dt)
        @test res.boundary_conditions.evaporation[1] ≈
              to_SI(0.5400000810623169, MM_PER_DT; dt_val = dt)
    end

    Wflow.close_files(model; delete_output = false)

    # test without lateral snow transport
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.snow_gravitational_transport__flag = false

    model = Wflow.Model(config)
    Wflow.run_timestep!(model)

    @testset "lateral snow transport off" begin
        snow = model.land.snow
        @test snow.variables.snow_storage[5] ≈ to_SI(3.7686103651001375, MM)
        @test mean(snow.variables.snow_storage) ≈ to_SI(0.03801972367609432, MM)
        @test mean(snow.variables.snow_water) ≈ to_SI(0.0025756728488273866, MM)
        @test mean(snow.variables.swe) ≈ to_SI(0.0405953965249217, MM)
    end

    # test without snow model
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.snow__flag = false
    pop!(config.output.netcdf_grid.variables, "snowpack_dry_snow__leq_depth")
    pop!(config.output.netcdf_grid.variables, "snowpack_liquid_water__depth")
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)

    @testset "snow model not included" begin
        snow = model.land.snow
        @test typeof(model.land.snow) == Wflow.NoSnowModel
    end
    Wflow.close_files(model; delete_output = false)
end

@testitem "Set pit and multithreading multiple basins" begin
    using Dates
    # test for setting a pit and multithreading multiple basins (by setting 2 extra pits
    # resulting in 3 basins)
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.pit__flag = true
    config.input.basin_pit_location__mask = "wflow_pits"
    config.time.endtime = DateTime(2000, 1, 9)
    config.logging.loglevel = "info"

    model = Wflow.run(config)

    @testset "timing" begin
        # clock has been reset
        @test model.clock.time == Wflow.cftime(config.time.starttime, config.time.calendar)
        @test model.clock.iteration == 0
    end

    @testset "river flow at basin outlets and downstream of one pit" begin
        q = model.routing.river_flow.variables.q_av.average
        @test q[4009] ≈ 8.426694842173548 # pit/ outlet, CartesianIndex(141, 228)
        @test q[4020] ≈ 0.006370691658310787 # downstream of pit 4009, CartesianIndex(141, 229)
        @test q[2508] ≈ 131.40631419288573 # pit/ outlet
        @test q[5808] ≈ 0.11941498848157915  # pit/ outlet
    end
end

@testitem "Changing forcing and cyclic LAI parameter" begin
    using Wflow: to_SI, MM_PER_DT
    # Run unchanged
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    dt = Wflow.tosecond(model.clock.dt)

    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    precip = copy(model.land.atmospheric_forcing.precipitation)
    evap = copy(model.land.atmospheric_forcing.potential_evaporation)
    lai = copy(model.land.vegetation_parameters.leaf_area_index)
    res_evap = copy(
        model.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.evaporation,
    )

    # Run changed
    config.input.forcing["atmosphere_water__precipitation_volume_flux"] =
        Dict("scale" => 2.0, "netcdf_variable_name" => "precip")
    config.input.forcing["land_surface_water__potential_evaporation_volume_flux"] =
        Dict("scale" => 3.0, "offset" => 1.50, "netcdf_variable_name" => "pet")
    config.input.cyclic["vegetation__leaf_area_index"] =
        Dict("scale" => 1.6, "netcdf_variable_name" => "LAI")
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "changed dynamic parameters" begin
        res = model.routing.river_flow.boundary_conditions.reservoir
        land = model.land
        @test land.atmospheric_forcing.precipitation[2] / precip[2] ≈ 2.0f0
        @test (
            land.atmospheric_forcing.potential_evaporation[100] -
            to_SI(1.50, MM_PER_DT; dt_val = dt)
        ) / evap[100] ≈ 3.0f0
        @test land.vegetation_parameters.leaf_area_index[100] / lai[100] ≈ 1.6f0
        @test (
            res.boundary_conditions.evaporation[2] - to_SI(1.50, MM_PER_DT; dt_val = dt)
        ) / res_evap[2] ≈ 3.0f0
    end
end

@testitem "Cyclic river and reservoir external inflow (kinematic wave routing)" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()

    config.input.cyclic["river_water__external_inflow_volume_flow_rate"] = "inflow"
    config.input.cyclic["reservoir_water__external_inflow_volume_flow_rate"] = "reservoir_inflow"

    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "kinematic wave routing: river and reservoir external inflow (cyclic)" begin
        (; reservoir) = model.routing.river_flow.boundary_conditions
        @test model.routing.river_flow.boundary_conditions.external_inflow[44] ≈ 0.75
        @test model.routing.river_flow.boundary_conditions.actual_external_abstraction_av.average[44] ==
              0.0
        @test Wflow.get_average(model.routing.river_flow.variables.q_av)[44] ≈
              10.156053077341845
        @test reservoir.boundary_conditions.external_inflow[2] == -1.0
        @test reservoir.boundary_conditions.actual_external_abstraction_av.average[2] == 1.0
        @test reservoir.boundary_conditions.inflow.average[2] ≈ -0.9054049318713108
        @test reservoir.variables.outflow_av.average[2] ≈ 3.000999922024245
    end
end

@testitem "Cyclic river and reservoir external inflow (local inertial routing)" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()

    config.input.cyclic["river_water__external_inflow_volume_flow_rate"] = "inflow"
    config.input.cyclic["reservoir_water__external_inflow_volume_flow_rate"] = "reservoir_inflow"

    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "local inertial routing: river and reservoir external inflow (cyclic)" begin
        (; reservoir) = model.routing.river_flow.boundary_conditions
        @test model.routing.river_flow.boundary_conditions.external_inflow[44] ≈ 0.75
        @test model.routing.river_flow.boundary_conditions.actual_external_abstraction_av.average[44] ==
              0.0
        @test model.routing.river_flow.variables.q_av.average[44] ≈ 10.119602100591411
        @test reservoir.boundary_conditions.external_inflow[2] == -1.0
        @test reservoir.boundary_conditions.actual_external_abstraction_av.average[2] == 1.0
        @test reservoir.boundary_conditions.inflow.average[2] ≈ -0.9090882985601229
        @test reservoir.variables.outflow_av.average[2] ≈ 3.000999922022744
    end
end

@testitem "External negative inflow" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    config.dir_output = mktempdir()
    model.routing.river_flow.boundary_conditions.external_inflow[44] = -10.0
    (; actual_external_abstraction_av, external_inflow) =
        model.routing.river_flow.boundary_conditions
    (; q_av) = model.routing.river_flow.variables
    @testset "river external negative inflow" begin
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_av.average[44] ≈ 1.5965227273142157
        @test q_av.average[44] ≈ 1.4328200140906533
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_av.average[44] ≈ 5.416747895349456
        @test q_av.average[44] ≈ 4.000700150630984
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_av.average[44] ≈ 9.756847289612736
        @test external_inflow[44] == -10.0
        @test q_av.average[44] ≈ 7.275452569433593
    end
end

@testitem "Fixed forcing (precipitation = 2.5)" begin
    using Wflow: to_SI, MM_PER_DT
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.input.forcing["atmosphere_water__precipitation_volume_flux"] = 2.5
    model = Wflow.Model(config)
    dt = Wflow.tosecond(model.clock.dt)
    Wflow.load_fixed_forcing!(model)

    @testset "fixed precipitation forcing (initialize)" begin
        @test maximum(model.land.atmospheric_forcing.precipitation) ≈
              to_SI(2.5, MM_PER_DT; dt_val = dt)
        @test minimum(model.land.atmospheric_forcing.precipitation) ≈ 0.0
        @test all(
            isapprox.(
                model.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.precipitation,
                to_SI(2.5, MM_PER_DT; dt_val = dt),
            ),
        )
    end

    Wflow.run_timestep!(model)

    @testset "fixed precipitation forcing (first timestep)" begin
        @test maximum(model.land.atmospheric_forcing.precipitation) ≈
              to_SI(2.5, MM_PER_DT; dt_val = dt)
        @test minimum(model.land.atmospheric_forcing.precipitation) ≈ 0.0
        @test all(
            isapprox.(
                model.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.precipitation,
                to_SI(2.5, MM_PER_DT; dt_val = dt),
            ),
        )
    end

    Wflow.close_files(model; delete_output = false)
end

@testitem "Local-inertial option for river flow river_routing" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()

    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "river flow and depth (local inertial)" begin
        q = model.routing.river_flow.variables.q_av.average
        @test sum(q) ≈ 3692.1984174051863
        @test q[1622] ≈ 7.260572774917902e-5
        @test q[43] ≈ 11.248731903973153
        @test q[501] ≈ 3.163313530437266
        h = model.routing.river_flow.variables.h
        @test h[1622] ≈ 0.00191277920611667
        @test h[43] ≈ 0.447403557220214
        @test h[501] ≈ 0.37810388857945015
        q_channel = model.routing.river_flow.variables.q_channel_av.average
        @test q ≈ q_channel
    end
end

@testitem "External negative inflow local-inertial river flow" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    model.routing.river_flow.boundary_conditions.external_inflow[44] = -10.0
    (; actual_external_abstraction_av, external_inflow, reservoir) =
        model.routing.river_flow.boundary_conditions
    (; q_av) = model.routing.river_flow.variables
    @testset "river external negative inflow (local inertial)" begin
        Wflow.run_timestep!(model)
        @test Wflow.get_average(actual_external_abstraction_av)[44] ≈ 3.0038267868094692
        @test Wflow.get_average(q_av)[44] ≈ 0.0007694313441619584
        Wflow.run_timestep!(model)
        @test Wflow.get_average(actual_external_abstraction_av)[44] ≈ 9.420226369186521
        @test Wflow.get_average(q_av)[44] ≈ 0.009169625606576883
        Wflow.run_timestep!(model)
        @test Wflow.get_average(actual_external_abstraction_av)[44] ≈ 9.999999999999991
        @test Wflow.get_average(external_inflow)[44] == -10.0
        @test Wflow.get_average(q_av)[44] ≈ 6.917023458482667
    end
    Wflow.close_files(model; delete_output = false)
end

@testitem "Local-inertial option for river and overland flow" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-land-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()

    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "river and overland flow and depth (local inertial)" begin
        q = model.routing.river_flow.variables.q_av.average
        @test sum(q) ≈ 2415.8565080484427
        @test q[1622] ≈ 7.289953869913158e-5
        @test q[43] ≈ 5.235679976913814
        @test q[501] ≈ 1.395242045016401
        h = model.routing.river_flow.variables.h
        @test h[1622] ≈ 0.001915431360203216
        @test h[43] ≈ 0.30691170109894245
        @test h[501] ≈ 0.3039452016829932
        qx = model.routing.overland_flow.variables.qx
        qy = model.routing.overland_flow.variables.qy
        @test qx[[26, 35, 631]] ≈ [0.17466039055941732, 0.002077783779618663, 0.0]
        @test qy[[26, 35, 631]] ≈ [0.12287975269080074, 0.019011432839373375, 0.0]
        h = model.routing.overland_flow.variables.h
        @test h[[26, 35, 631]] ≈ [0.07163764409112827, 0.009015112665828464, 0.0]
    end

    Wflow.close_files(model; delete_output = false)
end

@testitem "Local-inertial option for river flow including 1D floodplain schematization" begin
    using Wflow: to_SI, M3_PER_DAY, MM, Unit
    M_PER_DAY = Unit(; m = 1, d = -1)
    M2_PER_DAY = Unit(; m = 2, d = -1)
    tomlpath = joinpath(@__DIR__, "sbm_river-floodplain-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)

    (; flow_length, flow_length) = model.domain.river.parameters
    fp = model.routing.river_flow.floodplain.parameters.profile
    river = model.routing.river_flow
    dh = diff(fp.depth)
    Δv = diff(fp.storage[:, 3])
    Δa = diff(fp.a[:, 3])

    @testset "river flow (local inertial) floodplain schematization" begin
        # floodplain geometry checks (index 3)
        @test fp.storage[:, 3] ≈ [0.0, 8641.0, 19011.0, 31685.0, 51848.0, 80653.0]
        @test fp.width[:, 3] ≈ [
            30.0,
            99.28617594254938,
            119.15260323159785,
            145.6258527827648,
            231.6754039497307,
            330.9730700179533,
        ]
        @test fp.p[:, 3] ≈ [
            69.28617594254938,
            70.28617594254938,
            91.15260323159785,
            118.62585278276481,
            205.6754039497307,
            305.9730700179533,
        ]
        @test fp.a[:, 3] ≈ [
            0.0,
            49.64308797127469,
            109.21938958707361,
            182.032315978456,
            297.8700179533214,
            463.35655296229805,
        ]
        @test dh .* fp.width[2:end, 3] * flow_length[3] ≈ Δv
        @test fp.a[:, 3] * flow_length[3] ≈ fp.storage[:, 3]
        # flood depth from flood storage (8000.0)
        flood_vol = 8000.0
        river.variables.storage[3] = flood_vol + river.parameters.bankfull_storage[3]
        i1, i2 = Wflow.interpolation_indices(flood_vol, fp.storage[:, 3])
        @test (i1, i2) == (1, 2)
        flood_depth = Wflow.flood_depth(fp, flood_vol, flow_length[3], 3)
        @test flood_depth ≈ 0.46290938548779076
        @test (flood_depth - fp.depth[i1]) * fp.width[i2, 3] * flow_length[3] +
              fp.storage[i1, 3] ≈ flood_vol
        # flood depth from flood storage (12000.0)
        flood_vol = 12000.0
        river.variables.storage[3] = flood_vol + river.parameters.bankfull_storage[3]
        i1, i2 = Wflow.interpolation_indices(flood_vol, fp.storage[:, 3])
        @test (i1, i2) == (2, 3)
        flood_depth = Wflow.flood_depth(fp, flood_vol, flow_length[3], 3)
        @test flood_depth ≈ 0.6619575699132112
        @test (flood_depth - fp.depth[i1]) * fp.width[i2, 3] * flow_length[3] +
              fp.storage[i1, 3] ≈ flood_vol
        # test extrapolation of segment
        flood_vol = 95000.0
        river.variables.storage[3] = flood_vol + river.parameters.bankfull_storage[3]
        i1, i2 = Wflow.interpolation_indices(flood_vol, fp.storage[:, 3])
        @test (i1, i2) == (6, 6)
        flood_depth = Wflow.flood_depth(fp, flood_vol, flow_length[3], 3)
        @test flood_depth ≈ 2.749036625585836
        @test (flood_depth - fp.depth[i1]) * fp.width[i2, 3] * flow_length[3] +
              fp.storage[i1, 3] ≈ flood_vol
        river.variables.storage[3] = 0.0 # reset storage
        # flow area and wetted perimeter based on hf
        h = 0.5
        i1, i2 = Wflow.interpolation_indices(h, fp.depth)
        @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
              49.64308797127469
        @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 70.28617594254938
        h = 1.5
        i1, i2 = Wflow.interpolation_indices(h, fp.depth)
        @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
              182.032315978456
        @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 118.62585278276481
        h = 1.7
        i1, i2 = Wflow.interpolation_indices(h, fp.depth)
        @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
              228.36739676840216
        @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 119.02585278276482
        h = 3.2
        i1, i2 = Wflow.interpolation_indices(h, fp.depth)
        @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
              695.0377019748654
        @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 307.3730700179533
        h = 4.0
        i1, i2 = Wflow.interpolation_indices(h, fp.depth)
        @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
              959.816157989228
        @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 308.9730700179533
        @test Wflow.flow_area(fp.width[i2, 4], fp.a[i1, 4], fp.depth[i1], h) ≈
              407.6395313908081
        @test Wflow.wetted_perimeter(fp.p[i1, 4], fp.depth[i1], h) ≈ 90.11775307900271
    end

    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "river flow (local inertial) with floodplain schematization simulation" begin
        q = Wflow.get_average(model.routing.river_flow.variables.q_av)
        @test sum(q) ≈ 3682.763689938466
        @test q[1622] ≈ 7.260572782240141e-5
        @test q[43] ≈ 11.248731903973168
        @test q[501] ≈ 3.121950941199987
        @test q[5808] ≈ 0.0022169051490145645
        h = model.routing.river_flow.variables.h
        @test h[1622] ≈ 0.0019127792060060795
        @test h[43] ≈ 0.44740355722021324
        @test h[501] ≈ 0.37317174589196656
        @test h[5808] ≈ 0.0073332218703272825
    end

    # set boundary condition local inertial routing from netCDF file
    config.input.static["model_boundary_condition_river__length"] = "riverlength_bc"
    config.input.static["model_boundary_condition_river_bank_water__depth"] = "riverdepth_bc"
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "change boundary condition for local inertial routing (including floodplain)" begin
        q = Wflow.get_average(model.routing.river_flow.variables.q_av)
        @test sum(q) ≈ 3682.950549280907
        @test q[1622] ≈ 7.260572782240141e-5
        @test q[43] ≈ 11.248731903973168
        @test q[501] ≈ 3.121950941199987
        @test q[5808] ≈ 0.05466707293198311
        h = model.routing.river_flow.variables.h
        @test h[1622] ≈ 0.0019127792060060795
        @test h[43] ≈ 0.44740355722021324
        @test h[501] ≈ 0.37317174589196656
        @test h[5808] ≈ 2.0000266345482336
    end
    Wflow.close_files(model; delete_output = false)

    # test different ksat profiles
    @testset "ksat profiles (SBM)" begin
        tomlpath = joinpath(@__DIR__, "sbm_config.toml")
        function get_config(profile)
            config = Wflow.Config(tomlpath)
            config.dir_output = mktempdir()
            config.model.saturated_hydraulic_conductivity_profile = profile
            config.input.static["soil_layer_water__vertical_saturated_hydraulic_conductivity"] = "kv"
            config.input.static["soil_exponential_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth"] =
                400.0
            config.input.static["soil_layered_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth"] =
                400.0
            config
        end

        i = 100

        @testset "exponential profile" begin
            config = get_config(Wflow.VerticalConductivityProfile.exponential)
            config.dir_output = mktempdir()
            model = Wflow.Model(config)
            (; soil) = model.land
            (; kv_profile) = soil.parameters
            (; subsurface_flow) = model.routing
            z = soil.variables.zi[i]
            kvfrac = soil.parameters.kvfrac
            kv_z = Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, z, i, 2)
            @test kv_z ≈ kvfrac[i][2] * kv_profile.kv_0[i] * exp(-kv_profile.f[i] * z)
            @test subsurface_flow.variables.ssfmax[i] ≈ to_SI(28.32720603576582, M2_PER_DAY)
            @test subsurface_flow.variables.ssf[i] ≈ to_SI(11683.330684556406, M3_PER_DAY)
        end

        @testset "exponential constant profile" begin
            config = get_config(Wflow.VerticalConductivityProfile.exponential_constant)
            config.dir_output = mktempdir()
            model = Wflow.Model(config)
            (; soil) = model.land
            (; kv_profile) = soil.parameters
            (; subsurface_flow) = model.routing
            z = soil.variables.zi[i]
            kvfrac = soil.parameters.kvfrac
            kv_z = Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, z, i, 2)
            @test kv_z ≈
                  kvfrac[i][2] *
                  kv_profile.exponential.kv_0[i] *
                  exp(-kv_profile.exponential.f[i] * z)
            kv_400 = Wflow.hydraulic_conductivity_at_depth(
                kv_profile,
                kvfrac,
                to_SI(400.0, MM),
                i,
                2,
            )
            kv_1000 = Wflow.hydraulic_conductivity_at_depth(
                kv_profile,
                kvfrac,
                to_SI(1000.0, MM),
                i,
                3,
            )
            @test kv_400 ≈ kv_1000
            @test all(kv_profile.z_exp .== to_SI(400.0, MM))
            @test subsurface_flow.variables.ssfmax[i] ≈ to_SI(49.38558575188426, M2_PER_DAY)
            @test subsurface_flow.variables.ssf[i] ≈ to_SI(24810.460986497365, M3_PER_DAY)
        end

        @testset "layered profile" begin
            config = get_config(Wflow.VerticalConductivityProfile.layered)
            config.dir_output = mktempdir()
            model = Wflow.Model(config)
            (; soil) = model.land
            (; kv_profile) = soil.parameters
            (; subsurface_flow) = model.routing
            z = soil.variables.zi[i]
            kvfrac = soil.parameters.kvfrac
            @test Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, z, i, 2) ≈
                  kv_profile.kv[i][2]
            Wflow.kh_layered_profile!(soil, subsurface_flow, kv_profile)
            @test subsurface_flow.parameters.kh_profile.kh[i] ≈
                  to_SI(47.508932674632355, M_PER_DAY)
            @test subsurface_flow.variables.ssfmax[i] ≈
                  to_SI(30.237094380100316, M2_PER_DAY)
            @test subsurface_flow.variables.ssf[i] ≈ to_SI(14546.518932613191, M3_PER_DAY)
        end

        config = get_config(Wflow.VerticalConductivityProfile.layered_exponential)

        @testset "layered exponential profile" begin
            config.dir_output = mktempdir()
            model = Wflow.Model(config)
            (; soil) = model.land
            (; kv_profile) = soil.parameters
            (; subsurface_flow) = model.routing
            z = soil.variables.zi[i]
            kvfrac = soil.parameters.kvfrac
            @test Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, z, i, 2) ≈
                  kv_profile.kv[i][2]
            @test kv_profile.nlayers_kv[i] == 2
            Wflow.kh_layered_profile!(soil, subsurface_flow, kv_profile)
            @test subsurface_flow.parameters.kh_profile.kh[i] ≈
                  to_SI(33.76026208801769, M_PER_DAY)
            @test all(kv_profile.z_layered[1:10] .== to_SI(400.0, MM))
            @test subsurface_flow.variables.ssfmax[i] ≈ to_SI(23.4840490395906, M2_PER_DAY)
            @test subsurface_flow.variables.ssf[i] ≈ to_SI(10336.88327617503, M3_PER_DAY)
        end

        @testset "river flow layered exponential profile" begin
            config.dir_output = mktempdir()
            model = Wflow.Model(config)
            Wflow.run_timestep!(model)
            Wflow.run_timestep!(model)
            q = Wflow.get_average(model.routing.river_flow.variables.q_av)
            @test sum(q) ≈ 3032.617045063436
            @test q[1622] ≈ 0.0006987386860043929
            @test q[43] ≈ 8.710529495056752
        end

        Wflow.close_files(model; delete_output = false)
    end
end

@testitem "run wflow sbm" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.time.endtime = "2000-01-05"
    Wflow.run(config)
end

@testitem "water balance sbm (kinematic wave routing)" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.water_mass_balance__flag = true
    model = Wflow.Model(config)
    (; land_water_balance, routing) = model.mass_balance
    (; overland_water_balance, river_water_balance, subsurface_water_balance) = routing
    Wflow.run_timestep!(model)
    @testset "water balance first timestep" begin
        @test all(e -> abs(e) < 1e-9, land_water_balance.error)
        @test all(re -> abs(re) < 1e-9, land_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, overland_water_balance.error)
        @test all(re -> abs(re) < 6.6e11, overland_water_balance.relative_error)
        inds = findall(
            x -> x > 1e-3,
            Wflow.get_average(model.routing.overland_flow.variables.q_av),
        )
        @test all(re -> abs(re) < 1e-9, routing.overland_water_balance.relative_error[inds])
        @test all(e -> abs(e) < 1.e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1e-9, subsurface_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 1e-9, land_water_balance.error)
        @test all(re -> abs(re) < 1e-9, land_water_balance.relative_error)
        @test all(e -> abs(e) < 1.e-9, routing.overland_water_balance.error)
        @test all(re -> abs(re) < 5.4e11, routing.overland_water_balance.relative_error)
        inds = findall(
            x -> x > 1e-3,
            Wflow.get_average(model.routing.overland_flow.variables.q_av),
        )
        @test all(re -> abs(re) < 1e-9, routing.overland_water_balance.relative_error[inds])
        @test all(e -> abs(e) < 3e-5, river_water_balance.error)
        @test all(re -> abs(re) < 12.2, river_water_balance.relative_error)
        inds = findall(
            x -> x > 1e-3,
            Wflow.get_average(model.routing.river_flow.variables.q_av),
        )
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error[inds])
        @test all(e -> abs(e) < 1e-9, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1e-9, subsurface_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end

@testitem "water balance river local inertial routing" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.water_mass_balance__flag = true
    model = Wflow.Model(config)
    (; river_water_balance) = model.mass_balance.routing
    Wflow.run_timestep!(model)
    @testset "water balance first timestep" begin
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end

@testitem "water balance river local inertial routing with floodplain" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-floodplain-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.water_mass_balance__flag = true
    model = Wflow.Model(config)
    (; river_water_balance) = model.mass_balance.routing
    Wflow.run_timestep!(model)
    @testset "water balance first timestep" begin
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end

@testitem "water balance river and land local inertial routing" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-land-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.water_mass_balance__flag = true
    model = Wflow.Model(config)
    (; overland_water_balance) = model.mass_balance.routing
    Wflow.run_timestep!(model)
    @testset "water balance first timestep" begin
        @test all(e -> abs(e) < 1e-9, overland_water_balance.error)
        @test all(re -> abs(re) < 1e-9, overland_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 1e-9, overland_water_balance.error)
        @test all(re -> abs(re) < 1e-9, overland_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end
