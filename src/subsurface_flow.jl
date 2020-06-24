Base.@kwdef struct LateralSSF{T}
    kh₀::Vector{T}                          # Horizontal hydraulic conductivity at soil surface [mm Δt⁻¹]
    f::Vector{T}                            # A scaling parameter [mm⁻¹] (controls exponential decline of kh₀)
    soilthickness::Vector{T}                # Soil thickness [mm]
    θₑ::Vector{T}                           # Effective porosity [-]
    Δt::T=1.0                               # model time step
    βₗ::Vector{T}                           # Slope [m m⁻¹]
    dl::Vector{T}                           # Drain length [mm]
    dw::Vector{T}                           # Flow width [mm]
    zi::Vector{T} = fill(mv, length(f))     # Pseudo-water table depth [mm] (top of the saturated zone)
    exfiltwater::Vector{T} = fill(mv, length(f))  # Exfiltration [mm]  (groundwater above surface level, saturated excess conditions)
    recharge::Vector{T} = fill(mv, length(f))     # Net recharge to saturated store [mm]
    ssf::Vector{T} =
        ((kh₀ .* βₗ) ./ f) .* (exp.(-f .* zi) - exp.(-f .* soilthickness)) .* dw    # Subsurface flow [mm³ Δt⁻¹]
    ssfmax::Vector{T} = ((kh₀ .* βₗ) ./ f) .* (1.0 .- exp.(-f .* soilthickness))     # Maximum subsurface flow [mm² Δt⁻¹]
    to_river::Vector{T} = zeros(length(f))  # Part of subsurface flow [mm³ Δt⁻¹] that flows to the river
    pits::Vector{Int64} = zeros(Int64,length(f))
end

"""
    statenames(::Type{LateralSSF})

Returns Array{Symbol,1} for extracting model state fields from SBM struct.
"""
function statenames(::Type{LateralSSF})

    # depends on ini file settings (optional: glaciers, snow, irrigation)
    states = [:ssf]
    # TODO: (warm) states read from netcdf file or cold state (reinit=1, setting in ini file)

end

function update(ssf::LateralSSF, dag, toposort, frac_toriver, river)
    for v in toposort
        upstream_nodes = inneighbors(dag, v)
        if Bool(river[v]) & (ssf.pits[v] == 0)
            ssfin = isempty(upstream_nodes) ? 0.0 :
                sum(ssf.ssf[i] * (1.0 - frac_toriver[i]) for i in upstream_nodes if ssf.pits[i] == 0)
            ssf.to_river[v] = isempty(upstream_nodes) ? 0.0 :
                sum(ssf.ssf[i] * frac_toriver[i] for i in upstream_nodes if ssf.pits[i] == 0)
        elseif Bool(river[v]) & (ssf.pits[v] == 1)
            ssf.to_river[v] = isempty(upstream_nodes) ? 0.0 : sum(ssf.ssf[i] for i in upstream_nodes)
            ssfin = 0.0
        else
            ssfin = isempty(upstream_nodes) ? 0.0 : sum(ssf.ssf[i] for i in upstream_nodes)
        end
        ssf.ssf[v], ssf.zi[v], ssf.exfiltwater[v] = kinematic_wave_ssf(
            ssfin,
            ssf.ssf[v],
            ssf.zi[v],
            ssf.recharge[v],
            ssf.kh₀[v],
            ssf.βₗ[v],
            ssf.θₑ[v],
            ssf.f[v],
            ssf.soilthickness[v],
            ssf.Δt,
            ssf.dl[v],
            ssf.dw[v],
            ssf.ssfmax[v],
        )
    end
end
