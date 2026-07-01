@testitem "Run SBM" begin
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
    flush(model.writer.csv_writer.output_io)  # ensure the buffer is written fully to disk
    @testset "CSV output" begin
        row = csv_first_row(model.writer.csv_writer.output_path)

        @test row.time == DateTime("2000-01-02T00:00:00")
        @test row.Q ≈ 6.504946776799171
        @test row.storage ≈ 2.753500258633802e7
        @test row.temp_bycoord ≈ 2.390000104904175
        @test row.vwc_layer2_bycoord ≈ 0.25938707367146907
        @test row.temp_byindex ≈ 2.390000104904175
        @test row.Q_6336050 ≈ 0.0031759508929959544
        @test row.Q_6336510 ≈ 0.01278682587610146
        @test row.Q_6836100 ≈ 0.008741053540441504
        @test row.Q_6336500 ≈ 0.002406203339582039
        @test row.Q_6836190 ≈ 0.0025349488458642046
        @test row.Q_6336800 ≈ 0.004368075580441462
        @test row.Q_6336900 ≈ 0.003934294374817718
        @test row.Q_6336930 ≈ 0.01381216902016
        @test row.Q_6336910 ≈ 0.0030740079361179225
        @test row.Q_6136500 ≈ 0.00043545328783658126
        @test row.Q_6136520 ≈ 0.000412024150923924
        @test row.Q_6136150 ≈ 0.0029661013246156565
        @test row.Q_6136151 ≈ 0.0027847974999332813
        @test row.Q_6136160 ≈ 3.154622606318124
        @test row.Q_6136202 ≈ 1.1754540303915215
        @test row.recharge_1 ≈ -0.002257181032501202
    end

    @testset "NetCDF scalar output" begin
        ds = model.writer.scalar_writer.output_dataset
        @test ds["time"][1] == DateTime("2000-01-02T00:00:00")
        @test ds["Q"][:][1:20] ≈ [
            0.6442139,
            1.6144257,
            1.6131837,
            1.5442443,
            12.546244,
            4.767698,
            0.09143856,
            4.166856,
            0.010489026,
            5.6183286,
            5.608439,
            0.0059947344,
            0.0066055492,
            0.0063990145,
            0.004162496,
            1.2043678,
            0.0022157491,
            1.7636175,
            7.7863965,
            2.709964,
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

        @test test_means(
            interception.parameters,
            Dict(:evaporation_to_precipitation_ratio => 0.2381398056521972),
        )
        @test test_means(
            interception.parameters.vegetation_parameter_set,
            Dict(
                :crop_coefficient => 1.0,
                :storage_wood => 0.00019583714217432435,
                :canopy_gap_fraction => 0.5076891913529586,
                :storage_specific_leaf => 8.942700997221185e-5,
                :leaf_area_index => 1.0624380387861794,
                :maximum_canopy_storage => 0.0002888534265722341,
                :light_extinction_coefficient => 0.6741606033862981,
                :rooting_depth => 0.3701426593311398,
            ),
        )
        @test test_means(
            interception.variables,
            Dict(
                :stemflow => 1.1144419542415207e-9,
                :canopy_storage => 0.0,
                :canopy_potevap => 3.636931747082379e-9,
                :interception_rate => 3.6329026885493017e-9,
                :throughfall => 1.8659960684569295e-8,
            ),
        )
    end

    @testset "First timestep: snow" begin
        (; snow) = model.land

        @test snow.variables.snow_storage[5] ≈ 0.003768513390588815
        @test snow.parameters.temperature_threshold_snowfall[50063] ≈ 273.14999999999998

        @test test_means(
            snow.parameters,
            Dict(
                :water_holding_capacity => 0.1,
                :degree_day_factor => 4.3478356481481481e-8,
                :temperature_threshold_snowfall => 273.14999999999998,
                :temperature_threshold_melt => 273.14999999999998,
                :temperature_interval_snowfall => 2.0,
            ),
        )
        @test test_means(
            snow.variables,
            Dict(
                :snow_water_equivalent => 4.059539652492171e-5,
                :runoff => 1.9304548512364964e-8,
                :snow_in => 4.268605599818533e-15,
                :snow_melt => 0.0,
                :snow_out => 4.268605599818533e-15,
                :snow_storage => 3.801972367609432e-5,
                :snow_water => 2.5756728488273866e-6,
            ),
        )
    end

    @testset "First timestep: runoff" begin
        (; runoff) = model.land

        @test test_means(
            runoff.boundary_conditions,
            Dict(
                :waterdepth_river => 0.0,
                :water_flux_surface => 1.9304548512364964e-8,
                :waterdepth_land => 0.0,
            ),
        )
        @test test_means(
            runoff.variables,
            Dict(
                :actual_open_water_evaporation_river => 0.0,
                :runoff_river => 2.5374936291954614e-10,
                :net_runoff_river => 2.5374936291954614e-10,
                :actual_open_water_evaporation_land => 0.0,
                :runoff_land => 1.749511189351324e-11,
            ),
        )
    end

    @testset "First timestep: soil" begin
        (; soil) = model.land

        @test soil.variables.total_storage[50063] ≈ 0.559886196102425
        @test soil.variables.total_storage[429] ≈ 0.5989948259600009 # river cell
        @test soil.parameters.theta_s[50063] ≈ 0.48755401372909546
        @test soil.parameters.theta_r[50063] ≈ 0.15943120419979095

        @test test_means(
            soil.boundary_conditions,
            Dict(
                :potential_soilevaporation => 3.6593990196756077e-9,
                :water_flux_surface => 1.9033304037551905e-8,
                :potential_transpiration => 4.029058533077763e-12,
            ),
        )
        @test test_means(
            soil.parameters,
            Dict(
                :theta_s => 0.4409211971535584,
                :soil_thickness => 1.837833713668117,
                :h1 => 0.0,
                :cap_n => 2.0,
                :cumulative_layer_depth => SVector((
                    0.0,
                    0.1,
                    0.39997241953958815,
                    1.1566640303676966,
                    1.9842496012545644,
                )),
                :brooks_corey_exponent => SVector((
                    9.428788533549843,
                    9.821687673542042,
                    10.240060684266773,
                    10.24248827673959,
                )),
                :h4 => -160,
                :w_soil => 0.11249999999999985,
                :air_entry_pressure => -0.10000000000000001,
                :vertical_hydraulic_conductivity_factor => SVector((1.0, 1.0, 1.0, 1.0)),
                :h3_low => -10,
                :soil_fraction => 0.49508064804427593,
                :cap_hmax => 2.0,
                :actual_layer_thickness => SVector((
                    0.1,
                    0.29997241953958815,
                    0.7566640303676965,
                    0.7842496012545644,
                )),
                :compacted_soil_area_fraction => 0.0129736199199603,
                :h3_high => -4,
                :infiltration_capacity_compacted_soil => 5.7870370370370404e-8,
                :infiltration_capacity_soil => 4.8456216809760705e-6,
                :cf_soil => 0.038,
                :theta_fc => 0.30670581485821735,
                :maximum_leakage => 0.0,
                :h2 => -1,
                :wet_root_distribution_parameter => -500000.0,
                :alpha_h1 => 1.0,
                :soil_water_capacity => 0.5057182462428663,
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
                :transfer => 1.3761602372536785e-16,
                :drainable_water_depth => 0.21027893230146349,
                :infiltration_excess => 0.0,
                :actual_infiltration => 1.9033304037551905e-8,
                :unsaturated_layer_depth => SVector((
                    0.0016280832346593193,
                    0.0006460862670135344,
                    0.00023103976228101625,
                    0.0,
                )),
                :h3 => -10,
                :recharge => -2.6124780468763912e-11,
                :unsaturated_layer_thickness => SVector((
                    0.0993509423527741,
                    0.19022782953886028,
                    0.0898015331913427,
                    0.0,
                )),
                :soil_evaporation => 2.4568273505664539e-10,
                :net_runoff => 4.6686362399155382e-10,
                :actual_infiltration_compacted_soil => 2.4817112828889956e-10,
                :unsaturated_store_capacity => 0.074421716916843156,
                :actual_evaporation_unsaturated_store => 0.0,
                :volumetric_water_content_root_zone => 0.24099999255961077,
                :runoff => 4.6686362399155382e-10,
                :exfiltration_saturated_water => 4.4936851209804049e-10,
                :total_soil_water_storage => 0.43129652932602325,
                :soil_evaporation_saturated_zone => 2.211499042312533e-11,
                :unsaturated_store_depth => 0.0025052092639538705,
                :root_zone_storage => 0.029037374113238562,
                :infiltration => 1.9033304037551905e-8,
                :water_table_depth => 0.27873782910028516,
                :relative_volumetric_water_content_root_zone => 54.55735722947444,
                :total_storage => 0.43155732894128607,
                :actual_leakage => 0.0,
                :actual_evapotranspiration => 3.8825953512676096e-9,
                :saturated_water_depth => 0.42879132006206927,
                :transpiration => 4.0099276616623107e-12,
                :actual_infiltration_soil => 1.8785132909263007e-8,
                :actual_capillary_flux => 0.0,
                :relative_volumetric_water_content =>
                    [41.86055714460927, 63.23920727815791, 99.92432179078311, 100.0],
                :soil_surface_temperature => 282.39366827336522,
                :f_infiltration_reduction => 1.0,
                :actual_evaporation_saturated_zone => 4.0099276616623107e-12,
                :volumetric_water_content => [
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
                :soil_thickness => 1.8378337136681173,
                :specific_yield => 0.134215382295341,
                :horizontal_to_vertical_hydraulic_conductivity_ratio => 100.0,
            ),
        )
        @test test_means(
            subsurface_flow.parameters.kh_profile,
            Dict(
                :hydraulic_conductivity_scale_parameter => 3.303715296489842,
                :kh_0 => 0.00048456216809760719,
            ),
        )
        @test test_means(
            subsurface_flow.variables,
            Dict(
                :water_table_depth => 0.2787378291002851,
                :to_river_average => 0.00097200732537938705,
                :q_in => 0.01395379236964431,
                :storage => 118195.04455118616,
                :q_max => 3.6448833695946121e-5,
                :q => 0.01492660788818179,
                :exfiltwater_average => 4.493685120980405e-10,
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
        @test test_means(
            overland_flow.variables,
            Dict(:to_river_average => 7.2848055222914435e-5),
        )
        @test test_means(
            overland_flow.variables.flow,
            Dict(
                :qin => 0.00089086351171593122,
                :storage => 16.914628942221,
                :h => 2.969784057393903e-5,
                :qin_average => 0.00059166067952359249,
                :qlat => 2.6373852250863608e-7,
                :q => 0.0010160565843809637,
                :q_average => 0.00066450873474650687,
            ),
        )
    end

    @testset "First timestep: river routing" begin
        (; river_flow) = model.routing

        @test test_means(
            river_flow.boundary_conditions,
            Dict(
                :actual_external_abstraction_average => 0.0,
                :external_inflow => 0.0,
                :abstraction => 0.0,
                :inwater => 0.010231632451838139,
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
                :qin => 0.29084047475454816,
                :storage => 931.35597624692161,
                :h => 0.019029209285808049,
                :qin_average => 0.14990440354072712,
                :qlat => 1.7628453567888727e-5,
                :q => 0.29029797394167772,
                :q_average => 0.14935645293415184,
            ),
        )
    end

    # run the second timestep
    Wflow.run_timestep!(model)

    @testset "Second timestep: snow" begin
        (; snow) = model.land

        @test snow.variables.snow_storage[5] ≈ 0.003843412524052313

        @test test_means(
            snow.variables,
            Dict(
                :snow_water_equivalent => 3.583038054459421e-5,
                :runoff => 1.0425158326992447e-9,
                :snow_in => 4.371605031101146e-15,
                :snow_melt => 6.19370856351947e-11,
                :snow_out => 4.371605031101141e-15,
                :snow_storage => 3.461317061870947e-5,
                :snow_water => 1.2172099258847414e-6,
            ),
        )
    end

    @testset "Second timestep: interception" begin
        (; interception) = model.land

        @test test_means(
            interception.variables,
            Dict(
                :canopy_storage => 0.0,
                :stemflow => 8.709638687055903e-11,
                :canopy_potevap => 3.049183740864458e-9,
                :interception_rate => 8.071145121444295e-10,
                :throughfall => 9.002687979082284e-10,
            ),
        )
    end

    @testset "Second timestep: runoff" begin
        (; runoff) = model.land

        @test test_means(
            runoff.variables,
            Dict(
                :actual_open_water_evaporation_river => 7.314418628894106e-11,
                :runoff_river => 1.3009557837753908e-11,
                :net_runoff_river => -6.013462845118714e-11,
                :actual_open_water_evaporation_land => 4.264605809694288e-12,
                :runoff_land => 1.6461563909942252e-12,
            ),
        )
    end

    @testset "Second timestep: soil" begin
        (; soil) = model.land

        @test soil.parameters.theta_s[50063] ≈ 0.48755401372909546
        @test soil.parameters.theta_r[50063] ≈ 0.15943120419979095
        # Total storage is affected by all modules, not just soil
        @test soil.variables.total_storage[50063] ≈ 0.5599633189802773
        @test soil.variables.total_storage[429] ≈ 0.6209840187839495  # river cell

        @test test_means(
            soil.variables,
            Dict(
                :transfer => 1.5567149457984323e-10,
                :drainable_water_depth => 0.20950390590615132,
                :infiltration_excess => 0.0,
                :actual_infiltration => 1.0256861822245556e-9,
                :unsaturated_layer_depth => SVector((
                    0.0017064168532989077,
                    0.0010896449436965638,
                    0.0006820793025801899,
                    1.0165257547301471e-7,
                )),
                :h3 => -10,
                :recharge => -1.9009727917816603e-9,
                :unsaturated_layer_thickness => SVector((
                    0.099409226268861529,
                    0.19270470155860653,
                    0.17182135284713667,
                    0.030507142365925574,
                )),
                :soil_evaporation => 2.1687761811933483e-10,
                :net_runoff => 2.5569736341548426e-9,
                :actual_infiltration_compacted_soil => 1.6788849455513738e-11,
                :unsaturated_store_capacity => 0.074899232539565163,
                :actual_evaporation_unsaturated_store => 5.8410546538801915e-11,
                :volumetric_water_content_root_zone => 0.24143642302182922,
                :runoff => 2.5612382399645373e-9,
                :exfiltration_saturated_water => 2.5574181473276022e-9,
                :total_soil_water_storage => 0.4308190137033012,
                :soil_evaporation_saturated_zone => 2.0584284119007964e-11,
                :excess_water_soil => 2.1739362459410241e-12,
                :unsaturated_store_depth => 0.0034782427521511345,
                :root_zone_storage => 0.029196404187600344,
                :infiltration => 1.0256861822245556e-9,
                :water_table_depth => 0.2831713993535127,
                :relative_volumetric_water_content_root_zone => 54.646709985548235,
                # Total storage is affected by all modules, not just soil
                :total_storage => 0.43144808509560123,
                :actual_leakage => 0.0,
                :actual_evapotranspiration => 3.1762875630098438e-9,
                :saturated_water_depth => 0.42734077095115,
                :transpiration => 2.0748866407854196e-9,
                :actual_infiltration_soil => 1.0088973327690417e-9,
                :actual_capillary_flux => 1.9583907995878329e-11,
                :relative_volumetric_water_content => [
                    42.213917397074624,
                    63.201354128983084,
                    99.7692631420567,
                    99.99996110706128,
                ],
                :soil_surface_temperature => 281.70200654167354,
                :saturation_excess_water => 2.1739362459410189e-12,
                :f_infiltration_reduction => 1.0,
                :actual_evaporation_saturated_zone => 2.016476094246616e-9,
                :volumetric_water_content => [
                    0.18585985815835235,
                    0.27925583454823527,
                    0.4397834998887272,
                    0.4383913077433366,
                ],
            ),
        )
    end

    @testset "subsurface flow" begin
        q = model.routing.subsurface_flow.variables.q_cumulative
        @test sum(q) ≈ 6.250079949202134e7
        @test q[domain.land.network.order[1]] ≈ 699.3636285243076
        @test q[domain.land.network.order[end - 100]] ≈ 2395.6159482448143
        @test q[domain.land.network.order[end]] ≈ 287.61501877867994
    end

    @testset "Second timestep: overland routing" begin
        (; overland_flow) = model.routing

        q = overland_flow.variables.q_average
        @test q[26625] ≈ 0.0
        @test q[39308] ≈ 0.0
        @test q[domain.land.network.order[end]] ≈ Wflow.KIN_WAVE_MIN_FLOW

        @test test_means(
            overland_flow.boundary_conditions,
            Dict(:inwater => 0.001470763758351302),
        )
        @test test_means(
            overland_flow.variables,
            Dict(:to_river_average => 0.00064316539719498872),
        )
        @test test_means(
            overland_flow.variables.flow,
            Dict(
                :qin => 0.0069527815375203259,
                :storage => 88.41912734612589,
                :h => 0.00015596174801490019,
                :qin_average => 0.0047461547800561193,
                :qlat => 1.6550140762075119e-6,
                :q => 0.0080249579262317433,
                :q_average => 0.0053893201772511092,
            ),
        )
    end

    @testset "river flow" begin
        q = model.routing.river_flow.variables.q_average
        @test sum(q) ≈ 3684.852140874334
        @test q[1622] ≈ 0.00075026328229499
        @test q[43] ≈ 11.640174323692957
        @test q[domain.river.network.order[end]] ≈ 0.043976809930571785
    end

    @testset "reservoir simple" begin
        res = model.routing.river_flow.boundary_conditions.reservoir
        @test res.variables.outflow[1] ≈ 0.2174998614438593
        @test res.variables.outflow_average[1] ≈ 0.21749986282401396
        @test res.boundary_conditions.inflow_average[1] ≈ 0.0005130587494645402
        @test res.variables.storage[1] ≈ 2.751299001489657f7
        @test res.variables.actevap_cumulative[1] ≈ 0.00054000002145767148
        @test res.boundary_conditions.precipitation[1] ≈ 2.0833334161175625e-9
        @test res.boundary_conditions.evaporation[1] ≈ 6.250000287445232e-9
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
        @test snow.variables.snow_storage[5] ≈ 0.0037686103651001375
        @test mean(snow.variables.snow_storage) ≈ 3.801972367609432e-5
        @test mean(snow.variables.snow_water) ≈ 2.5756728488273866e-6
        @test mean(snow.variables.snow_water_equivalent) ≈ 4.05953965249217e-5
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
        (; snow) = model.land
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
        q = model.routing.river_flow.variables.q_average
        @test q[4009] ≈ 8.426693946008998 # pit/ outlet, CartesianIndex(141, 228)
        @test q[4020] ≈ 0.006370691658310787 # downstream of pit 4009, CartesianIndex(141, 229)
        @test q[2508] ≈ 131.41740124407315 # pit/ outlet
        @test q[5808] ≈ 0.11941506586989471  # pit/ outlet
    end
end

@testitem "Changing forcing and cyclic LAI parameter" begin
    # Run unchanged
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    dt = Wflow.tosecond(model.clock.dt)

    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    precipitation = copy(model.land.atmospheric_forcing.precipitation)
    potential_evaporation = copy(model.land.atmospheric_forcing.potential_evaporation)
    leaf_area_index = copy(model.land.vegetation_parameters.leaf_area_index)
    evaporation = copy(
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
        (; reservoir) = model.routing.river_flow.boundary_conditions
        (; land) = model
        @test land.atmospheric_forcing.precipitation[2] / precipitation[2] ≈ 2.0f0
        @test (land.atmospheric_forcing.potential_evaporation[100] - 1.736111111111111e-8) /
              potential_evaporation[100] ≈ 3.0f0
        @test land.vegetation_parameters.leaf_area_index[100] / leaf_area_index[100] ≈ 1.6f0
        @test (reservoir.boundary_conditions.evaporation[2] - 1.736111111111111e-8) /
              evaporation[2] ≈ 3.0f0
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
        @test model.routing.river_flow.boundary_conditions.actual_external_abstraction_average[44] ==
              0.0
        @test model.routing.river_flow.variables.q_average[44] ≈ 10.365364897089025
        @test reservoir.boundary_conditions.external_inflow[2] == -1.0
        @test reservoir.boundary_conditions.actual_external_abstraction_average[2] == 1.0
        @test reservoir.boundary_conditions.inflow_average[2] ≈ -0.9054042878240134
        @test reservoir.variables.outflow_average[2] ≈ 3.000999922024245
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
        @test model.routing.river_flow.boundary_conditions.actual_external_abstraction_average[44] ==
              0.0
        @test model.routing.river_flow.variables.q_average[44] ≈ 10.326210164434384
        @test reservoir.boundary_conditions.external_inflow[2] == -1.0
        @test reservoir.boundary_conditions.actual_external_abstraction_average[2] == 1.0
        @test reservoir.boundary_conditions.inflow_average[2] ≈ -0.9090852656640709
        @test reservoir.variables.outflow_average[2] ≈ 3.000999922022744
    end
end

@testitem "External negative inflow" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    config.dir_output = mktempdir()
    model.routing.river_flow.boundary_conditions.external_inflow[44] = -10.0
    (; actual_external_abstraction_average, external_inflow) =
        model.routing.river_flow.boundary_conditions
    (; q_average) = model.routing.river_flow.variables
    @testset "river external negative inflow" begin
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_average[44] ≈ 1.4800592763999205
        @test q_average[44] ≈ 1.2991010772732907
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_average[44] ≈ 5.509341540685494
        @test q_average[44] ≈ 4.122878209438459
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_average[44] ≈ 9.768399128401446
        @test external_inflow[44] == -10.0
        @test q_average[44] ≈ 7.295244400365228
    end
end

@testitem "Fixed forcing (precipitation = 2.5)" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.input.forcing["atmosphere_water__precipitation_volume_flux"] = 2.5
    model = Wflow.Model(config)
    dt = Wflow.tosecond(model.clock.dt)
    Wflow.load_fixed_forcing!(model)

    @testset "fixed precipitation forcing (initialize)" begin
        @test maximum(model.land.atmospheric_forcing.precipitation) ≈ 2.8935185185185185e-8
        @test minimum(model.land.atmospheric_forcing.precipitation) ≈ 0.0
        @test all(
            isapprox.(
                model.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.precipitation,
                2.8935185185185185e-8,
            ),
        )
    end

    Wflow.run_timestep!(model)

    @testset "fixed precipitation forcing (first timestep)" begin
        @test maximum(model.land.atmospheric_forcing.precipitation) ≈ 2.8935185185185185e-8
        @test minimum(model.land.atmospheric_forcing.precipitation) ≈ 0.0
        @test all(
            isapprox.(
                model.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.precipitation,
                2.8935185185185185e-8,
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
        q = model.routing.river_flow.variables.q_average
        @test sum(q) ≈ 3674.9826484892055
        @test q[1622] ≈ 7.266165539770211e-5
        @test q[43] ≈ 11.412913790413869
        @test q[501] ≈ 2.7237177111088933
        h = model.routing.river_flow.variables.h
        @test h[1622] ≈ 0.0019132839771304254
        @test h[43] ≈ 0.4517927241864443
        @test h[501] ≈ 0.3736919888456138
        q_channel = model.routing.river_flow.variables.q_channel_average
        @test q ≈ q_channel
    end
end

@testitem "External negative inflow local-inertial river flow" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    model.routing.river_flow.boundary_conditions.external_inflow[44] = -10.0
    (; actual_external_abstraction_average, external_inflow, reservoir) =
        model.routing.river_flow.boundary_conditions
    (; q_average) = model.routing.river_flow.variables
    @testset "river external negative inflow (local inertial)" begin
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_average[44] ≈ 2.7605062760705117
        @test q_average[44] ≈ 0.0007404199646763214
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_average[44] ≈ 9.632823305808886
        @test q_average[44] ≈ 0.008032294694728919
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_average[44] ≈ 9.999999999999991
        @test external_inflow[44] == -10.0
        @test q_average[44] ≈ 6.944679106931151
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
        q = model.routing.river_flow.variables.q_average
        @test sum(q) ≈ 2415.8565080484427
        @test q[1622] ≈ 7.289953869913158e-5
        @test q[43] ≈ 5.235679976913814
        @test q[501] ≈ 1.395242045016401
        h = model.routing.river_flow.variables.h
        @test h[1622] ≈ 0.001915431360203216
        @test h[43] ≈ 0.30691170109894245
        @test h[501] ≈ 0.3039452016829932
        (; qx, qy) = model.routing.overland_flow.variables
        @test qx[[26, 35, 631]] ≈ [0.17466039055941732, 0.002077783779618663, 0.0]
        @test qy[[26, 35, 631]] ≈ [0.12287975269080074, 0.019011432839373375, 0.0]
        (; h) = model.routing.overland_flow.variables
        @test h[[26, 35, 631]] ≈ [0.07163764409112827, 0.009015112665828464, 0.0]
    end

    Wflow.close_files(model; delete_output = false)
end

@testitem "Local-inertial option for river flow including 1D floodplain schematization" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-floodplain-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)

    (; flow_length, flow_length) = model.domain.river.parameters
    (; profile) = model.routing.river_flow.floodplain.parameters
    (; river_flow) = model.routing
    dh = diff(profile.depth)
    Δv = diff(profile.storage[:, 3])
    Δa = diff(profile.flow_area[:, 3])

    @testset "river flow (local inertial) floodplain schematization" begin
        # floodplain geometry checks (index 3)
        @test profile.storage[:, 3] ≈ [0.0, 8641.0, 19011.0, 31685.0, 51848.0, 80653.0]
        @test profile.width[:, 3] ≈ [
            30.0,
            99.28617594254938,
            119.15260323159785,
            145.6258527827648,
            231.6754039497307,
            330.9730700179533,
        ]
        @test profile.wetted_perimeter[:, 3] ≈ [
            69.28617594254938,
            70.28617594254938,
            91.15260323159785,
            118.62585278276481,
            205.6754039497307,
            305.9730700179533,
        ]
        @test profile.flow_area[:, 3] ≈ [
            0.0,
            49.64308797127469,
            109.21938958707361,
            182.032315978456,
            297.8700179533214,
            463.35655296229805,
        ]
        @test dh .* profile.width[2:end, 3] * flow_length[3] ≈ Δv
        @test profile.flow_area[:, 3] * flow_length[3] ≈ profile.storage[:, 3]
        # flood depth from flood storage (8000.0)
        flood_vol = 8000.0
        river_flow.variables.storage[3] =
            flood_vol + river_flow.parameters.bankfull_storage[3]
        i1, i2 = Wflow.interpolation_indices(flood_vol, profile.storage[:, 3])
        @test (i1, i2) == (1, 2)
        flood_depth = Wflow.flood_depth(profile, flood_vol, flow_length[3], 3)
        @test flood_depth ≈ 0.46290938548779076
        @test (flood_depth - profile.depth[i1]) * profile.width[i2, 3] * flow_length[3] +
              profile.storage[i1, 3] ≈ flood_vol
        # flood depth from flood storage (12000.0)
        flood_vol = 12000.0
        river_flow.variables.storage[3] =
            flood_vol + river_flow.parameters.bankfull_storage[3]
        i1, i2 = Wflow.interpolation_indices(flood_vol, profile.storage[:, 3])
        @test (i1, i2) == (2, 3)
        flood_depth = Wflow.flood_depth(profile, flood_vol, flow_length[3], 3)
        @test flood_depth ≈ 0.6619575699132112
        @test (flood_depth - profile.depth[i1]) * profile.width[i2, 3] * flow_length[3] +
              profile.storage[i1, 3] ≈ flood_vol
        # test extrapolation of segment
        flood_vol = 95000.0
        river_flow.variables.storage[3] =
            flood_vol + river_flow.parameters.bankfull_storage[3]
        i1, i2 = Wflow.interpolation_indices(flood_vol, profile.storage[:, 3])
        @test (i1, i2) == (6, 6)
        flood_depth = Wflow.flood_depth(profile, flood_vol, flow_length[3], 3)
        @test flood_depth ≈ 2.749036625585836
        @test (flood_depth - profile.depth[i1]) * profile.width[i2, 3] * flow_length[3] +
              profile.storage[i1, 3] ≈ flood_vol
        river_flow.variables.storage[3] = 0.0 # reset storage
        # flow area and wetted perimeter based on hf
        h = 0.5
        i1, i2 = Wflow.interpolation_indices(h, profile.depth)
        @test Wflow.flow_area(
            profile.width[i2, 3],
            profile.flow_area[i1, 3],
            profile.depth[i1],
            h,
        ) ≈ 49.64308797127469
        @test Wflow.wetted_perimeter(
            profile.wetted_perimeter[i1, 3],
            profile.depth[i1],
            h,
        ) ≈ 70.28617594254938
        h = 1.5
        i1, i2 = Wflow.interpolation_indices(h, profile.depth)
        @test Wflow.flow_area(
            profile.width[i2, 3],
            profile.flow_area[i1, 3],
            profile.depth[i1],
            h,
        ) ≈ 182.032315978456
        @test Wflow.wetted_perimeter(
            profile.wetted_perimeter[i1, 3],
            profile.depth[i1],
            h,
        ) ≈ 118.62585278276481
        h = 1.7
        i1, i2 = Wflow.interpolation_indices(h, profile.depth)
        @test Wflow.flow_area(
            profile.width[i2, 3],
            profile.flow_area[i1, 3],
            profile.depth[i1],
            h,
        ) ≈ 228.36739676840216
        @test Wflow.wetted_perimeter(
            profile.wetted_perimeter[i1, 3],
            profile.depth[i1],
            h,
        ) ≈ 119.02585278276482
        h = 3.2
        i1, i2 = Wflow.interpolation_indices(h, profile.depth)
        @test Wflow.flow_area(
            profile.width[i2, 3],
            profile.flow_area[i1, 3],
            profile.depth[i1],
            h,
        ) ≈ 695.0377019748654
        @test Wflow.wetted_perimeter(
            profile.wetted_perimeter[i1, 3],
            profile.depth[i1],
            h,
        ) ≈ 307.3730700179533
        h = 4.0
        i1, i2 = Wflow.interpolation_indices(h, profile.depth)
        @test Wflow.flow_area(
            profile.width[i2, 3],
            profile.flow_area[i1, 3],
            profile.depth[i1],
            h,
        ) ≈ 959.816157989228
        @test Wflow.wetted_perimeter(
            profile.wetted_perimeter[i1, 3],
            profile.depth[i1],
            h,
        ) ≈ 308.9730700179533
        @test Wflow.flow_area(
            profile.width[i2, 4],
            profile.flow_area[i1, 4],
            profile.depth[i1],
            h,
        ) ≈ 407.6395313908081
        @test Wflow.wetted_perimeter(
            profile.wetted_perimeter[i1, 4],
            profile.depth[i1],
            h,
        ) ≈ 90.11775307900271
    end

    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "river flow (local inertial) with floodplain schematization simulation" begin
        q = model.routing.river_flow.variables.q_average
        @test sum(q) ≈ 3665.2261366163852
        @test q[1622] ≈ 7.26616553794276e-5
        @test q[43] ≈ 11.412913790413914
        @test q[501] ≈ 2.684242648974844
        @test q[5808] ≈ 0.00221158569545013
        h = model.routing.river_flow.variables.h
        @test h[1622] ≈ 0.001913283977158018
        @test h[43] ≈ 0.45179272418645067
        @test h[501] ≈ 0.3685482876724819
        @test h[5808] ≈ 0.007318613650752014
    end

    # set boundary condition local inertial routing from netCDF file
    config.input.static["model_boundary_condition_river__length"] = "riverlength_bc"
    config.input.static["model_boundary_condition_river_bank_water__depth"] = "riverdepth_bc"
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "change boundary condition for local inertial routing (including floodplain)" begin
        q = model.routing.river_flow.variables.q_average
        @test sum(q) ≈ 3665.4127709930654
        @test q[1622] ≈ 7.26616553794276e-5
        @test q[43] ≈ 11.412913790413914
        @test q[501] ≈ 2.684242648974844
        @test q[5808] ≈ 0.05460501423141849
        h = model.routing.river_flow.variables.h
        @test h[1622] ≈ 0.001913283977158018
        @test h[43] ≈ 0.45179272418645067
        @test h[501] ≈ 0.3685482876724819
        @test h[5808] ≈ 2.0000269327006985
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
            z = soil.variables.water_table_depth[i]
            vertical_hydraulic_conductivity_factor =
                soil.parameters.vertical_hydraulic_conductivity_factor
            kv_z = Wflow.hydraulic_conductivity_at_depth(
                kv_profile,
                vertical_hydraulic_conductivity_factor,
                z,
                i,
                2,
            )
            @test kv_z ≈
                  vertical_hydraulic_conductivity_factor[i][2] *
                  kv_profile.kv_0[i] *
                  exp(-kv_profile.hydraulic_conductivity_scale_parameter[i] * z)
            @test subsurface_flow.variables.q_max[i] ≈ 0.00032786118096951182
            @test subsurface_flow.variables.q[i] ≈ 0.13522373477495839
        end

        @testset "exponential constant profile" begin
            config = get_config(Wflow.VerticalConductivityProfile.exponential_constant)
            config.dir_output = mktempdir()
            model = Wflow.Model(config)
            (; soil) = model.land
            (; kv_profile) = soil.parameters
            (; subsurface_flow) = model.routing
            z = soil.variables.water_table_depth[i]
            vertical_hydraulic_conductivity_factor =
                soil.parameters.vertical_hydraulic_conductivity_factor
            kv_z = Wflow.hydraulic_conductivity_at_depth(
                kv_profile,
                vertical_hydraulic_conductivity_factor,
                z,
                i,
                2,
            )
            @test kv_z ≈
                  vertical_hydraulic_conductivity_factor[i][2] *
                  kv_profile.exponential.kv_0[i] *
                  exp(-kv_profile.exponential.hydraulic_conductivity_scale_parameter[i] * z)
            kv_400 = Wflow.hydraulic_conductivity_at_depth(
                kv_profile,
                vertical_hydraulic_conductivity_factor,
                0.4,
                i,
                2,
            )
            kv_1000 = Wflow.hydraulic_conductivity_at_depth(
                kv_profile,
                vertical_hydraulic_conductivity_factor,
                1.0,
                i,
                3,
            )
            @test kv_400 ≈ kv_1000
            @test all(kv_profile.z_exp .== 0.4)
            @test subsurface_flow.variables.q_max[i] ≈ 0.00057159242768384559
            @test subsurface_flow.variables.q[i] ≈ 0.28715811326964541
        end

        @testset "layered profile" begin
            config = get_config(Wflow.VerticalConductivityProfile.layered)
            config.dir_output = mktempdir()
            model = Wflow.Model(config)
            (; soil) = model.land
            (; kv_profile) = soil.parameters
            (; subsurface_flow) = model.routing
            z = soil.variables.water_table_depth[i]
            vertical_hydraulic_conductivity_factor =
                soil.parameters.vertical_hydraulic_conductivity_factor
            @test Wflow.hydraulic_conductivity_at_depth(
                kv_profile,
                vertical_hydraulic_conductivity_factor,
                z,
                i,
                2,
            ) ≈ kv_profile.kv[i][2]
            Wflow.kh_layered_profile!(soil, subsurface_flow, kv_profile)
            @test subsurface_flow.parameters.kh_profile.kh[i] ≈ 0.00054987190595639297
            @test subsurface_flow.variables.q_max[i] ≈ 0.00034996637014004996
            @test subsurface_flow.variables.q[i] ≈ 0.16836248764598602
        end

        config = get_config(Wflow.VerticalConductivityProfile.layered_exponential)

        @testset "layered exponential profile" begin
            config.dir_output = mktempdir()
            model = Wflow.Model(config)
            (; soil) = model.land
            (; kv_profile) = soil.parameters
            (; subsurface_flow) = model.routing
            z = soil.variables.water_table_depth[i]
            vertical_hydraulic_conductivity_factor =
                soil.parameters.vertical_hydraulic_conductivity_factor
            @test Wflow.hydraulic_conductivity_at_depth(
                kv_profile,
                vertical_hydraulic_conductivity_factor,
                z,
                i,
                2,
            ) ≈ kv_profile.kv[i][2]
            @test kv_profile.nlayers_kv[i] == 2
            Wflow.kh_layered_profile!(soil, subsurface_flow, kv_profile)
            @test subsurface_flow.parameters.kh_profile.kh[i] ≈ 0.00039074377416687143
            @test all(kv_profile.z_layered[1:10] .== 0.4)
            @test subsurface_flow.variables.q_max[i] ≈ 0.00027180612314340973
            @test subsurface_flow.variables.q[i] ≈ 0.11963985273350729
        end

        @testset "river flow layered exponential profile" begin
            config.dir_output = mktempdir()
            model = Wflow.Model(config)
            Wflow.run_timestep!(model)
            Wflow.run_timestep!(model)
            q = model.routing.river_flow.variables.q_average
            @test sum(q) ≈ 3024.5079459603826
            @test q[1622] ≈ 0.0006987116204789493
            @test q[43] ≈ 8.767295043872874
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
        inds = findall(x -> x > 1e-3, model.routing.overland_flow.variables.q_average)
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
        inds = findall(x -> x > 1e-3, model.routing.overland_flow.variables.q_average)
        @test all(re -> abs(re) < 1e-9, routing.overland_water_balance.relative_error[inds])
        @test all(e -> abs(e) < 2e-4, river_water_balance.error)
        @test all(re -> abs(re) < 24.0, river_water_balance.relative_error)
        inds = findall(x -> x > 1e-3, model.routing.river_flow.variables.q_average)
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
        @test all(
            rating_curve_exponent -> abs(rating_curve_exponent) < 1e-9,
            river_water_balance.error,
        )
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(
            rating_curve_exponent -> abs(rating_curve_exponent) < 1e-9,
            river_water_balance.error,
        )
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
        @test all(
            rating_curve_exponent -> abs(rating_curve_exponent) < 1e-9,
            river_water_balance.error,
        )
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(
            rating_curve_exponent -> abs(rating_curve_exponent) < 1e-9,
            river_water_balance.error,
        )
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
        @test all(
            rating_curve_exponent -> abs(rating_curve_exponent) < 1e-9,
            overland_water_balance.error,
        )
        @test all(re -> abs(re) < 1e-9, overland_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(
            rating_curve_exponent -> abs(rating_curve_exponent) < 1e-9,
            overland_water_balance.error,
        )
        @test all(re -> abs(re) < 1e-9, overland_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end
