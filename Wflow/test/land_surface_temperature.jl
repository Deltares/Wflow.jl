@testitem "unit: land surface temperature" begin
    temperature = 293.15
    wind_measurement_height = 2.0
    wind_speed = 2.0
    canopy_height = 0.12
    actevap = 3.333333333333333e-8
    net_radiation = 300.0
    dt = 3600.0

    latent_heat_of_vaporization = Wflow.compute_latent_heat_of_vaporization(temperature)
    latent_heat_flux = Wflow.compute_latent_heat_flux(temperature, actevap, dt)
    sensible_heat_flux = Wflow.compute_sensible_heat_flux(net_radiation, latent_heat_flux)
    aerodynamic_resistance =
        Wflow.wind_and_aero_resistance(wind_speed, wind_measurement_height, canopy_height)

    lst = Wflow.compute_land_surface_temperature(
        sensible_heat_flux,
        aerodynamic_resistance,
        temperature,
    )

    @test latent_heat_of_vaporization ≈ 2.4535e6
    @test latent_heat_flux ≈ 81.78333333333332
    @test sensible_heat_flux ≈ 188.2166666666667
    @test aerodynamic_resistance ≈ 14.479951355887868
    @test lst ≈ 295.3637217404412
end
