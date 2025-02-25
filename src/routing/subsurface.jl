"Exponential depth profile of horizontal hydraulic conductivity at the soil surface"
struct KhExponential{T}
    # Horizontal hydraulic conductivity at soil surface [m d⁻¹]
    kh_0::Vector{T}
    # A scaling parameter [m⁻¹] (controls exponential decline of kh_0)
    f::Vector{T}
end

"Exponential constant depth profile of horizontal hydraulic conductivity"
struct KhExponentialConstant{T}
    # Exponential horizontal hydraulic conductivity profile type
    exponential::KhExponential
    # Depth [m] from soil surface for which exponential decline of kv_0 is valid
    z_exp::Vector{T}
end

"Layered depth profile of horizontal hydraulic conductivity"
struct KhLayered{T}
    # Horizontal hydraulic conductivity [m d⁻¹]
    kh::Vector{T}
end

"Struct for storing lateral subsurface flow model parameters"
@with_kw struct LateralSsfParameters{T, Kh}
    kh_profile::Kh              # Horizontal hydraulic conductivity profile type [-]  
    khfrac::Vector{T}           # A muliplication factor applied to vertical hydraulic conductivity `kv` [-]
    soilthickness::Vector{T}    # Soil thickness [m]
    theta_s::Vector{T}          # Saturated water content (porosity) [-]
    theta_r::Vector{T}          # Residual water content [-]   
    slope::Vector{T}            # Slope [m m⁻¹]
    flow_length::Vector{T}      # Flow length [m]
    flow_width::Vector{T}       # Flow width [m]
end

"Initialize lateral subsurface flow model parameters"
function LateralSsfParameters(
    dataset,
    config,
    indices;
    soil,
    slope,
    flow_length,
    flow_width,
)
    lens = lens_input_parameter(
        config,
        "subsurface_water__horizontal-to-vertical_saturated_hydraulic_conductivity_ratio",
    )
    khfrac = ncread(dataset, config, lens; sel = indices, defaults = 1.0, type = Float64)
    n_cells = length(khfrac)

    (; theta_s, theta_r, soilthickness) = soil
    soilthickness = soilthickness .* 0.001

    kh_profile_type =
        get(config.model, "saturated_hydraulic_conductivity_profile", "exponential")::String
    dt = Second(config.time.timestepsecs) / BASETIMESTEP
    if kh_profile_type == "exponential"
        (; kv_0, f) = soil.kv_profile
        kh_0 = khfrac .* kv_0 .* 0.001 .* dt
        kh_profile = KhExponential(kh_0, f .* 1000.0)
    elseif kh_profile_type == "exponential_constant"
        (; z_exp) = soil.kv_profile
        (; kv_0, f) = soil.kv_profile.exponential
        kh_0 = khfrac .* kv_0 .* 0.001 .* dt
        exp_profile = KhExponential(kh_0, f .* 1000.0)
        kh_profile = KhExponentialConstant(exp_profile, z_exp .* 0.001)
    elseif kh_profile_type == "layered" || kh_profile_type == "layered_exponential"
        kh_profile = KhLayered(fill(MISSING_VALUE, n_cells))
    end
    parameters = LateralSsfParameters(
        kh_profile,
        khfrac,
        soilthickness,
        theta_s,
        theta_r,
        slope,
        flow_length,
        flow_width,
    )
    return parameters
end

"Struct for storing lateral subsurface flow model variables"
@with_kw struct LateralSsfVariables{T}
    zi::Vector{T}           # Pseudo-water table depth [m] (top of the saturated zone)
    exfiltwater::Vector{T}  # Exfiltration [m Δt⁻¹] (groundwater above surface level, saturated excess conditions)
    recharge::Vector{T}     # Net recharge to saturated store [m² Δt⁻¹]
    ssf::Vector{T}          # Subsurface flow [m³ d⁻¹]
    ssfin::Vector{T}        # Inflow from upstream cells [m³ d⁻¹]
    ssfmax::Vector{T}       # Maximum subsurface flow [m² d⁻¹]
    to_river::Vector{T}     # Part of subsurface flow [m³ d⁻¹] that flows to the river
    storage::Vector{T}      # Subsurface storage [m³]
end

"Initialize lateral subsurface flow model variables"
function LateralSsfVariables(ssf, zi, xl, yl)
    n = length(zi)
    storage = @. (ssf.theta_s - ssf.theta_r) * (ssf.soilthickness - zi) * (xl * yl)
    variables = LateralSsfVariables(;
        zi,
        exfiltwater = fill(MISSING_VALUE, n),
        recharge = fill(MISSING_VALUE, n),
        ssf = fill(MISSING_VALUE, n),
        ssfin = fill(MISSING_VALUE, n),
        ssfmax = fill(MISSING_VALUE, n),
        to_river = zeros(n),
        storage,
    )
    return variables
end

"Struct for storing lateral subsurface flow model boundary conditions"
@with_kw struct LateralSsfBC{T}
    recharge::Vector{T} # Net recharge to saturated store [m² Δt⁻¹]
end

"Lateral subsurface flow model"
@with_kw struct LateralSSF{T, Kh} <: AbstractSubsurfaceFlowModel
    boundary_conditions::LateralSsfBC{T}
    parameters::LateralSsfParameters{T, Kh}
    variables::LateralSsfVariables{T}
end

"Initialize lateral subsurface flow model"
function LateralSSF(
    dataset,
    config,
    indices;
    soil,
    slope,
    flow_length,
    flow_width,
    x_length,
    y_length,
)
    parameters = LateralSsfParameters(
        dataset,
        config,
        indices;
        soil = soil.parameters,
        slope,
        flow_length,
        flow_width,
    )
    zi = 0.001 * soil.variables.zi
    variables = LateralSsfVariables(parameters, zi, x_length, y_length)
    boundary_conditions = LateralSsfBC(; recharge = fill(MISSING_VALUE, length(zi)))
    ssf = LateralSSF(; boundary_conditions, parameters, variables)
    return ssf
end

"Update lateral subsurface model for a single timestep"
function update!(model::LateralSSF, network, dt)
    (;
        order_of_subdomains,
        order_subdomain,
        subdomain_indices,
        upstream_nodes,
        area,
        frac_to_river,
    ) = network

    (; recharge) = model.boundary_conditions
    (; ssfin, ssf, ssfmax, to_river, zi, exfiltwater, storage) = model.variables
    (; slope, theta_s, theta_r, soilthickness, flow_length, flow_width, kh_profile) =
        model.parameters

    ns = length(order_of_subdomains)
    for k in 1:ns
        threaded_foreach(eachindex(order_of_subdomains[k]); basesize = 1) do i
            m = order_of_subdomains[k][i]
            for (n, v) in zip(subdomain_indices[m], order_subdomain[m])
                # for a river cell without a reservoir or lake part of the upstream
                # subsurface flow goes to the river (frac_to_river) and part goes to the
                # subsurface flow reservoir (1.0 - frac_to_river) upstream nodes with a
                # reservoir or lake are excluded
                ssfin[v] = sum_at(
                    i -> ssf[i] * (1.0 - frac_to_river[i]),
                    upstream_nodes[n],
                    eltype(ssfin),
                )
                to_river[v] = sum_at(
                    i -> ssf[i] * frac_to_river[i],
                    upstream_nodes[n],
                    eltype(to_river),
                )
                ssf[v], zi[v], exfiltwater[v] = kinematic_wave_ssf(
                    ssfin[v],
                    ssf[v],
                    zi[v],
                    recharge[v],
                    slope[v],
                    theta_s[v] - theta_r[v],
                    soilthickness[v],
                    dt,
                    flow_length[v],
                    flow_width[v],
                    ssfmax[v],
                    kh_profile,
                    v,
                )
                storage[v] =
                    (theta_s[v] - theta_r[v]) * (soilthickness[v] - zi[v]) * area[v]
            end
        end
    end
    return nothing
end

"""
Struct for storing groundwater exchange variables for coupling with an external groundwater
model.
"""
@with_kw struct GroundwaterExchangeVariables{T}
    exfiltwater::Vector{T}  # Exfiltration [m Δt⁻¹] (groundwater above surface level, saturated excess conditions)
    zi::Vector{T}           # Pseudo-water table depth [m] (top of the saturated zone)
    to_river::Vector{T}     # Part of subsurface flow [m³ d⁻¹] that flows to the river
    ssf::Vector{T}          # Subsurface flow [m³ d⁻¹]
end

"Initialize groundwater exchange variables"
function GroundwaterExchangeVariables(n)
    variables = GroundwaterExchangeVariables{Float64}(;
        exfiltwater = fill(MISSING_VALUE, n),
        zi = fill(MISSING_VALUE, n),
        to_river = fill(MISSING_VALUE, n),
        ssf = zeros(n),
    )
    return variables
end

"Groundwater exchange"
@with_kw struct GroundwaterExchange{T} <: AbstractSubsurfaceFlowModel
    variables::GroundwaterExchangeVariables{T}
end

"Initialize groundwater exchange"
function GroundwaterExchange(n)
    variables = GroundwaterExchangeVariables(n)
    ssf = GroundwaterExchange{Float64}(; variables)
    return ssf
end

# wrapper methods
get_water_depth(model::Union{LateralSSF, GroundwaterExchange}) = model.variables.zi
get_exfiltwater(model::Union{LateralSSF, GroundwaterExchange}) = model.variables.exfiltwater

get_flux_to_river(model::Union{LateralSSF, GroundwaterExchange}) =
    model.variables.to_river ./ tosecond(BASETIMESTEP) # [m³ s⁻¹]