@testitem "unit: land surface temperature" begin
    temperature = 293.15
    wind_measurement_height = 2.0
    wind_speed = 2.0
    skin_layer_height = 0.12
    actevap = 3.333333333333333e-8
    net_radiation = 300.0
    dt = 3600.0

    d0 = 2.0 / 3.0 * skin_layer_height
    z0m = 0.123 * skin_layer_height
    z0h = 0.1 * z0m

    latent_heat_of_vaporization = Wflow.compute_latent_heat_of_vaporization(temperature)
    latent_heat_flux = Wflow.compute_latent_heat_flux(temperature, actevap, dt)
    sensible_heat_flux = Wflow.compute_sensible_heat_flux(net_radiation, latent_heat_flux)
    aerodynamic_resistance = Wflow.compute_aerodynamic_resistance(
        wind_speed,
        wind_measurement_height,
        skin_layer_height,
        d0,
        z0m,
        z0h,
    )

    lst = Wflow.compute_land_surface_temperature(
        sensible_heat_flux,
        aerodynamic_resistance,
        temperature,
    )

    @test latent_heat_of_vaporization ≈ 2.4535e6
    @test latent_heat_flux ≈ 81.78333333333332
    @test sensible_heat_flux ≈ 188.2166666666667
    @test aerodynamic_resistance ≈ 103.83203500394342
    @test lst ≈ 309.0240335235324
end
