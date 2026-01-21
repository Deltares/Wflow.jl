# From Eq. 7:2.2.24 in https://swat.tamu.edu/media/99192/swat2009-theory.pdf
# [m s⁻¹ mm⁻²]
const STOKES_FACTOR = 411 / 3600

"""
Compute the sediment fall velocity [m s⁻¹] given the mean sediment diameter d50 [m]
Based on Stokes' law and assuming 22 °C and sediment density 1.2 t m⁻³
Eq. 7:2.2.24 in https://swat.tamu.edu/media/99192/swat2009-theory.pdf
"""
fall_velocity(d50) = STOKES_FACTOR * from_SI(d50, MM)^2

"""
    reservoir_deposition_camp(
        input,
        q,
        waterlevel,
        wb_area,
        wb_trapping_efficiency,
        dm,
        slope,
    )

Deposition of sediment in reservoirs from Camp 1945.

# Arguments
- `input` (sediment input [t dt⁻¹ = kg s⁻¹])
- `q` (discharge [m³ dt⁻¹ => m³ s⁻¹])
- `waterlevel` (water level [m])
- `reservoir_area` (reservoir area [m²])
- `res_trapping_efficiency` (reservoir trapping efficiency [-])
- `dm` (mean diameter [m])
- `slope` (slope [-])

# Output
- `deposition` (deposition [t dt⁻¹ => kg s⁻¹])
"""
function reservoir_deposition_camp(
    input::Float64,
    q::Float64,
    waterlevel::Float64,
    reservoir_area::Float64,
    res_trapping_efficiency::Float64,
    dm::Float64,
    slope::Float64,
)
    # Compute critical velocity
    # [m s⁻¹] = [m³ s⁻¹] / [m²]
    reservoir_critical_velocity = q / reservoir_area
    # Natural deposition
    # [kg s⁻¹] = [kg s⁻¹] * min([-], [m s⁻¹] / [m s⁻¹])
    deposition = input * min(1.0, fall_velocity(dm) / reservoir_critical_velocity)

    # Check if particles are traveling in suspension or bed load using Rouse number
    # [m]
    dsuspf = sqrt(
        1.2 * 0.41 * sqrt(GRAVITATIONAL_ACCELERATION * waterlevel * slope) / STOKES_FACTOR,
    )
    # If bed load, we have extra deposition depending on the reservoir type
    if dm > dsuspf
        deposition = max(deposition, res_trapping_efficiency * input)
    end

    return deposition
end
