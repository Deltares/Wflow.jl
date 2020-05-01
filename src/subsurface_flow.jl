const mv = NaN

Base.@kwdef struct SubSurfaceFlow{T,N}
    k₀::Vector{T}               # Horizontal hydraulic conductivity at soil surface [mm Δt⁻¹]
    f::Vector{T}                # A scaling parameter [mm⁻¹] (controls exponential decline of k₀)
    soilthickness::Vector{T}    # Soil thickness [mm]
    θₑ::Vector{T}               # Effective porosity [-]
    Δt::Dates.Second
    βₗ::Vector{T}               # Slope [m m⁻¹]
    dl::Vector{T}               # Drain length [mm]
    dw::Vector{T}                # Flow width [mm]
    zi::Vector{T} = fill(mv, N)              # Pseudo-water table depth [mm] (top of the saturated zone)
    exfiltwater::Vector{T} = fill(mv, N)
    ssfin::Vector{T} = fill(mv, N)
    subsurfaceflow::Vector{T} = ((k₀ .* βₗ) ./ f) .* (exp.(-f .* zi) - exp.(-f .* soilthickness)) .* dw
    ssfmax::Vector{T} = ((k₀ .* βₗ) ./ f) .* (exp.(-f .* zi) - exp.(-f .* soilthickness)) .* dw
end

"""
    statenames(::Type{SubSurfaceFlow})

Returns Array{Symbol,1} for extracting model state fields from SBM struct.
"""
function statenames(::Type{SubSurfaceFlow})

    # depends on ini file settings (optional: glaciers, snow, irrigation)
    states = [:subsurfaceflow]
    #TODO: (warm) states read from netcdf file or cold state (reinit=1, setting in ini file)

end
