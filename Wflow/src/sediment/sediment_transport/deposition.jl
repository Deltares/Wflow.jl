# From Eq. 7:2.2.24 in https://swat.tamu.edu/media/99192/swat2009-theory.pdf
# [m s⁻¹ mm⁻²]
const STOKES_FACTOR = 411 / 3600

"""
    reservoir_deposition_camp(
        input,
        q,
        waterlevel,
        res_area,
        res_trapping_efficiency,
        dm,
        slope,
    )

Deposition of sediment in reservoirs from Camp 1945.

# Arguments
- `input` (sediment input [t Δt⁻¹])
- `q` (discharge [m³ Δt⁻¹])
- `waterlevel` (water level [m])
- `res_area` (reservoir area [m²])
- `res_trapping_efficiency` (reservoir trapping efficiency [-])
- `dm` (mean diameter [m])
- `slope` (slope [-])

# Output
- `deposition` (deposition [t Δt⁻¹])
"""
function reservoir_deposition_camp(
    input::Float64,
    q::Float64,
    waterlevel::Float64,
    res_area::Float64,
    res_trapping_efficiency::Float64,
    dm::Float64,
    slope::Float64,
)
    # Compute critical velocity
    vcres = q / res_area
    DCres = STOKES_FACTOR / vcres
    # Natural deposition
    deposition = input * min(1.0, (DCres * (dm / 1000)^2))

    # Check if particles are travelling in suspension or bed load using Rouse number
    dsuspf =
        1e3 * sqrt(
            1.2 * 0.41 * sqrt(GRAVITATIONAL_ACCELERATION * waterlevel * slope) /
            STOKES_FACTOR,
        )
    # If bed load, we have extra deposition depending on the reservoir type
    if dm > dsuspf
        deposition = max(deposition, res_trapping_efficiency * input)
    end

    return deposition
end
