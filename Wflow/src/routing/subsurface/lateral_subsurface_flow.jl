"Struct for storing lateral subsurface flow model variables"
@kwdef struct LateralSsfVariables
    n::Int
    # Pseudo-water table depth [m] (top of the saturated zone)
    water_table_depth::Vector{Float64}
    # Hydraulic head [m]
    head::Vector{Float64}
    # Cumulative exfiltration [m] (groundwater above surface level, saturated excess conditions)
    exfiltwater_cumulative::Vector{Float64} = zeros(n)
    # Average exfiltration [m s⁻¹] (groundwater above surface level, saturated excess conditions)
    exfiltwater_average::Vector{Float64} = zeros(n)
    # Subsurface flow [m³ s⁻¹]
    q::Vector{Float64} = fill(MISSING_VALUE, n)
    # Cumulative subsurface flow [m³] for model timestep Δt
    q_cumulative::Vector{Float64} = zeros(n)
    # Average subsurface flow [m³ s⁻¹] for model timestep Δt
    q_average::Vector{Float64} = zeros(n)
    # Inflow from upstream cells [m³ s⁻¹]
    q_in::Vector{Float64} = fill(MISSING_VALUE, n)
    # cumulative inflow from upstream cells [m³] for model timestep dt
    q_in_cumulative::Vector{Float64} = zeros(n)
    # Average inflow from upstream cells [m³ s⁻¹] for model timestep dt
    q_in_average::Vector{Float64} = zeros(n)
    # Maximum subsurface flow [m s⁻¹]
    q_max::Vector{Float64} = fill(MISSING_VALUE, n)
    # Cumulative of the part of subsurface flow [m³ s⁻¹] that flows to the river
    to_river_cumulative::Vector{Float64} = zeros(n)
    # Average of the part of subsurface flow [m³ s⁻¹] that flows to the river
    to_river_average::Vector{Float64} = zeros(n)
    # Net flow for boundaries subsurface flow [m³ s⁻¹]
    q_net_bnds::Vector{Float64} = fill(MISSING_VALUE, n)
    # Cumulative net flow (total) [m³]
    q_net_cumulative::Vector{Float64} = zeros(n)
    # Average net flow (total) [m³ s⁻¹]
    q_net_average::Vector{Float64} = zeros(n)
    # Subsurface storage that can be released [m³]
    storage::Vector{Float64}
end

"Struct for storing lateral subsurface flow model parameters"
@kwdef struct LateralSsfParameters{Kh}
    # Horizontal hydraulic conductivity profile type [-]
    kh_profile::Kh
    # A multiplication factor applied to vertical hydraulic conductivity `kv` [-]
    horizontal_to_vertical_hydraulic_conductivity_ratio::Vector{Float64}
    # Soil thickness [m]
    soil_thickness::Vector{Float64}
    # Specific yield (theta_s - theta_fc) [-]
    specific_yield::Vector{Float64}
    # Area of cell [m²]
    area::Vector{Float64}
    # Top of subsurface flow layer [m]
    top::Vector{Float64}
end

"Lateral subsurface flow model"
@kwdef struct LateralSSFModel{Kh, B <: SubsurfaceFlowBC} <: AbstractSubsurfaceFlowModel
    timestepping::TimeStepping
    boundary_conditions::B
    parameters::LateralSsfParameters{Kh}
    variables::LateralSsfVariables
end

"Exponential depth profile of horizontal hydraulic conductivity at the soil surface"
struct KhExponential
    # Horizontal hydraulic conductivity at soil surface [m s⁻¹]
    kh_0::Vector{Float64}
    # A scaling parameter [m⁻¹] (controls exponential decline of kh_0)
    hydraulic_conductivity_scale_parameter::Vector{Float64}
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
    area::Vector{Float64},
)
    elevation = ncread(dataset, config, "land_surface__elevation", Routing; sel = indices)
    horizontal_to_vertical_hydraulic_conductivity_ratio = ncread(
        dataset,
        config,
        "subsurface_water__horizontal_to_vertical_saturated_hydraulic_conductivity_ratio",
        Routing;
        sel = indices,
    )

    (; theta_s, theta_fc, soil_thickness) = soil

    kh_profile_type = config.model.saturated_hydraulic_conductivity_profile

    if kh_profile_type == VerticalConductivityProfile.exponential
        (; kv_0, hydraulic_conductivity_scale_parameter) = soil.kv_profile
        kh_0 = horizontal_to_vertical_hydraulic_conductivity_ratio .* kv_0
        kh_profile = KhExponential(kh_0, hydraulic_conductivity_scale_parameter)
    elseif kh_profile_type == VerticalConductivityProfile.exponential_constant
        (; z_exp) = soil.kv_profile
        (; kv_0, hydraulic_conductivity_scale_parameter) = soil.kv_profile.exponential
        kh_0 = horizontal_to_vertical_hydraulic_conductivity_ratio .* kv_0
        exp_profile = KhExponential(kh_0, hydraulic_conductivity_scale_parameter)
        kh_profile = KhExponentialConstant(exp_profile, z_exp)
    elseif kh_profile_type == VerticalConductivityProfile.layered ||
           kh_profile_type == VerticalConductivityProfile.layered_exponential
        n_cells = length(horizontal_to_vertical_hydraulic_conductivity_ratio)
        kh_profile = KhLayered(fill(MISSING_VALUE, n_cells))
    end
    specific_yield = @. lower_bound_drainable_porosity(theta_s, theta_fc)
    ssf_parameters = LateralSsfParameters(;
        kh_profile,
        horizontal_to_vertical_hydraulic_conductivity_ratio,
        soil_thickness = copy(soil_thickness),
        specific_yield,
        area,
        top = elevation,
    )
    return ssf_parameters
end

"Initialize lateral subsurface flow model variables"
function LateralSsfVariables(ssf::LateralSsfParameters, water_table_depth::Vector{Float64})
    n = length(water_table_depth)
    storage = @. ssf.specific_yield * (ssf.soil_thickness - water_table_depth) * ssf.area
    head = ssf.top - water_table_depth
    variables = LateralSsfVariables(; n, water_table_depth, storage, head)
    return variables
end

"Initialize lateral subsurface flow model"
function LateralSSFModel(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil::SbmSoilModel,
)
    (; land, river, drain) = domain
    (; indices) = land.network
    (; area) = domain.land.parameters
    n = length(indices)
    timestepping = init_kinematic_wave_timestepping(config, n; domain = "subsurface")
    parameters = LateralSsfParameters(dataset, config, indices, soil.parameters, area)
    water_table_depth = copy(soil.variables.diagnostic.water_table_depth)
    variables = LateralSsfVariables(parameters, water_table_depth)
    recharge = RechargeModel(; n)
    if config.model.river_subsurface_exchange_head_based__flag
        river = GwfRiverModel(dataset, config, river.network.indices)
    else
        river = nothing
    end
    if config.model.drain__flag
        drain = DrainageModel(dataset, config, drain.network.indices)
    else
        drain = nothing
    end
    boundary_conditions = SubsurfaceFlowBC(; recharge, river, drain)
    ssf_model = LateralSSFModel(; timestepping, boundary_conditions, parameters, variables)
    return ssf_model
end

function update_fluxes!(subsurface_flow_model::LateralSSFModel, domain::Domain, dt::Float64)
    for bc in get_boundaries(subsurface_flow_model.boundary_conditions)
        indices = get_boundary_index(bc, domain)
        flux!(bc, subsurface_flow_model, indices, dt)
    end
    return nothing
end

function flux_to_river!(
    subsurface_flow_model::LateralSSFModel,
    domain::NetworkRiver,
    dt::Float64,
)
    (; to_river_average, to_river_cumulative) = subsurface_flow_model.variables
    (; river) = subsurface_flow_model.boundary_conditions
    if isnothing(river)
        @. to_river_average = to_river_cumulative / dt
    else
        inds = domain.land_indices
        to_river_average[inds] = -river.variables.flux_average
    end
    return nothing
end

function kinwave_subsurface_update!(
    subsurface_flow_model::LateralSSFModel,
    soil_model::SbmSoilModel,
    domain::Domain,
    dt::Float64,
)
    (; order_of_subdomains, order_subdomain, subdomain_indices, upstream_nodes) =
        domain.land.network
    (; flow_length, flow_width, area, flow_fraction_to_river, slope) =
        domain.land.parameters

    (;
        q_in,
        q_in_cumulative,
        q,
        q_cumulative,
        to_river_cumulative,
        water_table_depth,
        head,
        exfiltwater_cumulative,
        q_max,
        storage,
        q_net_bnds,
        q_net_cumulative,
    ) = subsurface_flow_model.variables
    (; specific_yield, top, soil_thickness, kh_profile) = subsurface_flow_model.parameters
    (; river) = subsurface_flow_model.boundary_conditions

    ns = length(order_of_subdomains)
    for hydraulic_conductivity in 1:ns
        threaded_foreach(
            eachindex(order_of_subdomains[hydraulic_conductivity]);
            basesize = 1,
        ) do i
            m = order_of_subdomains[hydraulic_conductivity][i]
            for (n, v) in zip(subdomain_indices[m], order_subdomain[m])
                if isnothing(river)
                    # for a river cell without a reservoir part of the upstream subsurface flow
                    # goes to the river (flow_fraction_to_river) and part goes to the subsurface
                    # flow reservoir (1.0 - flow_fraction_to_river) upstream nodes with a
                    # reservoir are excluded
                    q_in[v] = sum_at(
                        i -> q[i] * (1.0 - flow_fraction_to_river[i]),
                        upstream_nodes[n],
                    )
                    to_river_cumulative[v] +=
                        sum_at(i -> q[i] * flow_fraction_to_river[i], upstream_nodes[n]) *
                        dt
                else
                    q_in[v] = sum_at(i -> q[i], upstream_nodes[n])
                end

                q[v], water_table_depth[v], _exfiltwater, netflux = kinematic_wave_ssf(
                    q_in[v],
                    q[v],
                    water_table_depth[v],
                    q_net_bnds[v],
                    slope[v],
                    specific_yield[v],
                    soil_thickness[v],
                    dt,
                    flow_length[v],
                    flow_width[v],
                    q_max[v],
                    kh_profile,
                    soil_model,
                    v,
                )
                q_in_cumulative[v] += q_in[v] * dt
                q_cumulative[v] += q[v] * dt
                exfiltwater_cumulative[v] += _exfiltwater * dt
                q_net_cumulative[v] += netflux * area[v] * dt
                head[v] = top[v] - water_table_depth[v]
                storage[v] =
                    specific_yield[v] * (soil_thickness[v] - water_table_depth[v]) * area[v]
            end
        end
    end
end

"""
Update lateral subsurface model for a single timestep `dt`. Timestepping within `dt` is
either with a fixed timestep `dt_fixed` or adaptive.
"""
function update_subsurface_flow_model!(
    subsurface_flow_model::LateralSSFModel,
    soil_model::SbmSoilModel,
    domain::Domain,
    dt::Float64,
)
    (; to_river_cumulative) = subsurface_flow_model.variables
    (; adaptive) = subsurface_flow_model.timestepping

    to_river_cumulative .= 0.0
    set_flux_vars!(subsurface_flow_model)
    t = 0.0
    while t < dt
        subsurface_flow_model.variables.q_net_bnds .= 0.0
        dt_s =
            adaptive ? stable_timestep(subsurface_flow_model, domain.land) :
            subsurface_flow_model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        update_fluxes!(subsurface_flow_model, domain, dt_s)
        kinwave_subsurface_update!(subsurface_flow_model, soil_model, domain, dt_s)
        t += dt_s
    end
    average_flux_vars!(subsurface_flow_model, dt)
    flux_to_river!(subsurface_flow_model, domain.river.network, dt)
    return nothing
end

"""
Compute a stable timestep size for the kinematice wave method for a lateral subsurface flow
model using a nonlinear scheme (Chow et al., 1988).

A stable time step is computed for each vector element based on the Courant timestep size
criterion. Li et al. (1975) found that the nonlinear scheme is unconditionally stable and
that a wide range of dt/dx values can be used without loss of accuracy.
"""
function stable_timestep(subsurface_flow_model::LateralSSFModel, domain::DomainLand)
    (; water_table_depth) = subsurface_flow_model.variables
    (; specific_yield, kh_profile) = subsurface_flow_model.parameters
    (; flow_length, slope) = domain.parameters
    (; stable_timesteps, alpha_coefficient) = subsurface_flow_model.timestepping

    n = length(water_table_depth)
    stable_timesteps .= Inf
    hydraulic_conductivity = 0
    for i in 1:n
        if water_table_depth[i] > 0.0
            hydraulic_conductivity += 1
            c = ssf_celerity(
                water_table_depth[i],
                slope[i],
                specific_yield[i],
                kh_profile,
                i,
            )
            stable_timesteps[hydraulic_conductivity] = (flow_length[i] / c)
        end
    end

    dt_min = if hydraulic_conductivity > 0
        minimum(@view(stable_timesteps[1:hydraulic_conductivity]))
    else
        0.5
    end

    return dt_min * alpha_coefficient
end

get_flux_to_river(subsurface_flow_model::LateralSSFModel, inds::Vector{Int}) =
    subsurface_flow_model.variables.to_river_average[inds]

# wrapper method
get_water_depth(subsurface_flow_model::LateralSSFModel) =
    subsurface_flow_model.variables.water_table_depth
