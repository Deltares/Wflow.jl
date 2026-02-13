const MISSING_VALUE = NaN

using Parameters

# Minimal type definitions
struct SbmSoilModel
    variables::Any
end
struct AtmosphericForcing
    temperature::Any
    net_shortwave_radiation::Any
    net_longwave_radiation::Any
    net_radiation::Any
    wind_speed::Any
end
struct VegetationParameters
    canopy_height::Any
end
struct Config
    input::Any
    time::Any
end

include("../src/land_surface_temperature.jl")

# Test data
n = 1
soil_model = SbmSoilModel((; actevap = [2.0]))
atmospheric_forcing = AtmosphericForcing([20.0], [80.0], [20.0], [100.0], [2.0])  # T, SW, LW, wind
vegetation_parameters = VegetationParameters([0.12])
config = Config(Dict("wind_altitude" => 2.0), (; timestepsecs = 3600.0))

lst_model = LandSurfaceTemperatureModel(n)
update!(lst_model, soil_model, atmospheric_forcing, vegetation_parameters, config)

using Test

@testset "LST output" begin
    @test lst_model.variables.aerodynamic_resistance[1] ≈ 14.479951355887868 atol = 1e-8
    @test lst_model.variables.latent_heat_flux[1] ≈ 1.3630555555555555 atol = 1e-8
    @test lst_model.variables.sensible_heat_flux[1] ≈ 98.63694444444444 atol = 1e-8
    @test lst_model.variables.latent_heat_of_vaporization[1] ≈ 2453.5 atol = 1e-8
    @test lst_model.variables.land_surface_temperature[1] ≈ 21.16012440446662 atol = 1e-8
    @test isnan(lst_model.variables.net_radiation[1])
    @test isnan(lst_model.variables.net_shortwave_radiation[1])
    @test isnan(lst_model.variables.net_longwave_radiation[1])
end
