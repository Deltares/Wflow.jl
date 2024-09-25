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

Deposition of sediment in waterbodies from Camp 1945.

# Arguments
- `input` (sediment input [t Δt⁻¹])
- `q` (discharge [m³ Δt⁻¹])
- `waterlevel` (water level [m])
- `wb_area` (waterbody area [m²])
- `wb_trapping_efficiency` (waterbody trapping efficiency [-])
- `dm` (mean diameter [m])
- `slope` (slope [-])

# Output
- `deposition` (deposition [t Δt⁻¹])
"""
function reservoir_deposition_camp(
    input,
    q,
    waterlevel,
    wb_area,
    wb_trapping_efficiency,
    dm,
    slope,
)
    # Compute critical velocity
    vcres = q / wb_area
    DCres = 411 / 3600 / vcres
    # Natural deposition
    deposition = input * min(1.0, (DCres * (dm / 1000)^2))

    # Check if the particles was travelling in suspension or bed load using Rouse number
    dsuspf = 1e3 * (1.2 * 3600 * 0.41 / 411 * (9.81 * waterlevel * slope)^0.5)^0.5
    # If bed load, we have extra deposition depending on the reservoir type 
    if dm > dsuspf
        deposition = max(deposition, wb_trapping_efficiency * input)
    end

    return deposition
end