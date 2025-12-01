"Struct for storing lateral subsurface flow model variables"
@with_kw struct LateralSsfVariables
    zi::Vector{Float64}           # Pseudo-water table depth [m] (top of the saturated zone)
    exfiltwater::Vector{Float64}  # Exfiltration [m Δt⁻¹] (groundwater above surface level, saturated excess conditions)
    ssf::Vector{Float64}          # Subsurface flow [m³ d⁻¹]
    ssfin::Vector{Float64}        # Inflow from upstream cells [m³ d⁻¹]
    ssfmax::Vector{Float64}       # Maximum subsurface flow [m² d⁻¹]
    to_river::Vector{Float64}     # Part of subsurface flow [m³ d⁻¹] that flows to the river
    storage::Vector{Float64}      # Subsurface storage [m³]
end

"Struct for storing lateral subsurface flow model parameters"
@with_kw struct LateralSsfParameters{Kh}
    kh_profile::Kh                      # Horizontal hydraulic conductivity profile type [-]
    khfrac::Vector{Float64}             # A muliplication factor applied to vertical hydraulic conductivity `kv` [-]
    soilthickness::Vector{Float64}      # Soil thickness [m]
    theta_s::Vector{Float64}            # Saturated water content (porosity) [-]
    theta_r::Vector{Float64}            # Residual water content [-]
end

"Struct for storing lateral subsurface flow model boundary conditions"
@with_kw struct LateralSsfBC
    recharge::Vector{Float64} # Net recharge to saturated store [m² d⁻¹]
end

"Lateral subsurface flow model"
@with_kw struct LateralSSF{Kh} <: AbstractSubsurfaceFlowModel
    boundary_conditions::LateralSsfBC
    parameters::LateralSsfParameters{Kh}
    variables::LateralSsfVariables
end

"Exponential depth profile of horizontal hydraulic conductivity at the soil surface"
struct KhExponential
    # Horizontal hydraulic conductivity at soil surface [m d⁻¹]
    kh_0::Vector{Float64}
    # A scaling parameter [m⁻¹] (controls exponential decline of kh_0)
    f::Vector{Float64}
end

"Exponential constant depth profile of horizontal hydraulic conductivity"
struct KhExponentialConstant
    # Exponential horizontal hydraulic conductivity profile type
    exponential::KhExponential
    # Depth [m] from soil surface for which exponential decline of kv_0 is valid
    z_exp::Vector{Float64}
end

"Layered depth profile of horizontal hydraulic conductivity"
struct KhLayered
    # Horizontal hydraulic conductivity [m d⁻¹]
    kh::Vector{Float64}
end

"Initialize lateral subsurface flow model parameters"
function LateralSsfParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    soil::SbmSoilParameters,
)
    khfrac = ncread(
        dataset,
        config,
        "subsurface_water__horizontal_to_vertical_saturated_hydraulic_conductivity_ratio";
        optional = false,
        sel = indices,
        type = Float64,
    )

    (; theta_s, theta_r, soilthickness) = soil
    soilthickness = soilthickness .* 0.001

    kh_profile_type = config.model.saturated_hydraulic_conductivity_profile
    factor_dt = BASETIMESTEP / Second(config.time.timestepsecs)
    if kh_profile_type == VerticalConductivityProfile.exponential
        (; kv_0, f) = soil.kv_profile
        kh_0 = khfrac .* kv_0 .* 0.001 .* factor_dt
        kh_profile = KhExponential(kh_0, f .* 1000.0)
    elseif kh_profile_type == VerticalConductivityProfile.exponential_constant
        (; z_exp) = soil.kv_profile
        (; kv_0, f) = soil.kv_profile.exponential
        kh_0 = khfrac .* kv_0 .* 0.001 .* factor_dt
        exp_profile = KhExponential(kh_0, f .* 1000.0)
        kh_profile = KhExponentialConstant(exp_profile, z_exp .* 0.001)
    elseif kh_profile_type == VerticalConductivityProfile.layered ||
           kh_profile_type == VerticalConductivityProfile.layered_exponential
        n_cells = length(khfrac)
        kh_profile = KhLayered(fill(MISSING_VALUE, n_cells))
    end
    ssf_parameters =
        LateralSsfParameters(; kh_profile, khfrac, soilthickness, theta_s, theta_r)
    return ssf_parameters
end

"Initialize lateral subsurface flow model variables"
function LateralSsfVariables(
    ssf::LateralSsfParameters,
    zi::Vector{Float64},
    area::Vector{Float64},
)
    n = length(zi)
    storage = @. (ssf.theta_s - ssf.theta_r) * (ssf.soilthickness - zi) * area
    variables = LateralSsfVariables(;
        zi,
        exfiltwater = fill(MISSING_VALUE, n),
        ssf = fill(MISSING_VALUE, n),
        ssfin = fill(MISSING_VALUE, n),
        ssfmax = fill(MISSING_VALUE, n),
        to_river = zeros(n),
        storage,
    )
    return variables
end

"Initialize lateral subsurface flow model"
function LateralSSF(
    dataset::NCDataset,
    config::Config,
    domain::DomainLand,
    soil::SbmSoilModel,
)
    (; indices) = domain.network
    (; area) = domain.parameters
    parameters = LateralSsfParameters(dataset, config, indices, soil.parameters)
    zi = 0.001 * soil.variables.zi
    variables = LateralSsfVariables(parameters, zi, area)
    boundary_conditions = LateralSsfBC(; recharge = fill(MISSING_VALUE, length(zi)))
    ssf = LateralSSF(; boundary_conditions, parameters, variables)
    return ssf
end

"Update lateral subsurface model for a single timestep"
function update!(model::LateralSSF, domain::DomainLand, dt::Float64)
    (; order_of_subdomains, order_subdomain, subdomain_indices, upstream_nodes) =
        domain.network
    (; flow_length, flow_width, area, flow_fraction_to_river, slope) = domain.parameters

    (; recharge) = model.boundary_conditions
    (; ssfin, ssf, ssfmax, to_river, zi, exfiltwater, storage) = model.variables
    (; theta_s, theta_r, soilthickness, kh_profile) = model.parameters

    ns = length(order_of_subdomains)
    for k in 1:ns
        threaded_foreach(eachindex(order_of_subdomains[k]); basesize = 1) do i
            m = order_of_subdomains[k][i]
            for (n, v) in zip(subdomain_indices[m], order_subdomain[m])
                # for a river cell without a reservoir part of the upstream subsurface flow
                # goes to the river (flow_fraction_to_river) and part goes to the subsurface
                # flow reservoir (1.0 - flow_fraction_to_river) upstream nodes with a
                # reservoir are excluded
                ssfin[v] = sum_at(
                    i -> ssf[i] * (1.0 - flow_fraction_to_river[i]),
                    upstream_nodes[n],
                    eltype(ssfin),
                )
                to_river[v] = sum_at(
                    i -> ssf[i] * flow_fraction_to_river[i],
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

# wrapper methods
get_water_depth(model::LateralSSF) = model.variables.zi
get_exfiltwater(model::LateralSSF) = model.variables.exfiltwater

get_flux_to_river(model::LateralSSF, inds::Vector{Int}) =
    model.variables.to_river[inds] ./ tosecond(BASETIMESTEP) # [m³ s⁻¹]

get_inflow(model::LateralSSF) = model.variables.ssfin
get_outflow(model::LateralSSF) = model.variables.ssf
get_storage(model::LateralSSF) = model.variables.storage
