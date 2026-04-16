"Struct for storing lateral subsurface flow model variables"
@with_kw struct LateralSsfVariables
    n::Int
    # Pseudo-water table depth [m] (top of the saturated zone)
    zi::Vector{Float64}
    # Hydraulic head [m]
    head::Vector{Float64}
    # Exfiltration [m Δt⁻¹ => m s⁻¹] (groundwater above surface level, saturated excess conditions)
    exfiltwater::Vector{Float64} = fill(MISSING_VALUE, n)
    # Subsurface flow [m³ d⁻¹ => m³ s⁻¹]
    q::Vector{Float64} = fill(MISSING_VALUE, n)
    # Average subsurface flow [m³ d⁻¹ => m³ s⁻¹] for model timestep Δt
    q_av::AverageVector = AverageVector(; n)
    # Inflow from upstream cells [m³ d⁻¹ => m³ s⁻¹]
    q_in::Vector{Float64} = fill(MISSING_VALUE, n)
    # Average inflow from upstream cells [m³ d⁻¹ => m³ s⁻¹] for model timestep dt
    q_in_av::AverageVector = AverageVector(; n)
    # Maximum subsurface flow [m² d⁻¹ => m s⁻¹]
    q_max::Vector{Float64} = fill(MISSING_VALUE, n)
    # Part of subsurface flow [m³ d⁻¹ => m³ s⁻¹] that flows to the river
    to_river::AverageVector = AverageVector(; n)
    # Net flow for boundaries subsurface flow [m³ d⁻¹ => m³ s⁻¹]
    q_net_bnds::Vector{Float64} = fill(MISSING_VALUE, n)
    # Average net flow (total) [m³ d⁻¹ => m³ s⁻¹]
    q_net_av::AverageVector = AverageVector(; n)
    # Subsurface storage that can be released [m³]
    storage::Vector{Float64}
end

"Struct for storing lateral subsurface flow model parameters"
@with_kw struct LateralSsfParameters{Kh}
    # Horizontal hydraulic conductivity profile type [-]
    kh_profile::Kh
    # A multiplication factor applied to vertical hydraulic conductivity `kv` [-]
    khfrac::Vector{Float64}
    # Soil thickness [m]
    soilthickness::Vector{Float64}
    # Specific yield (theta_s - theta_fc) [-]
    specific_yield::Vector{Float64}
    # Area of cell [m²]
    area::Vector{Float64}
    # Top of subsurface flow layer [m]
    top::Vector{Float64}
end

"Lateral subsurface flow model"
@with_kw struct LateralSSFModel{Kh, B <: SubsurfaceFlowBC} <: AbstractSubsurfaceFlowModel
    timestepping::TimeStepping
    boundary_conditions::B
    parameters::LateralSsfParameters{Kh}
    variables::LateralSsfVariables
end

"Exponential depth profile of horizontal hydraulic conductivity at the soil surface"
struct KhExponential
    # Horizontal hydraulic conductivity at soil surface [m d⁻¹ => m s⁻¹]
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
    area::Vector{Float64},
)
    elevation = ncread(dataset, config, "land_surface__elevation", Routing; sel = indices)
    khfrac = ncread(
        dataset,
        config,
        "subsurface_water__horizontal_to_vertical_saturated_hydraulic_conductivity_ratio",
        Routing;
        sel = indices,
    )

    (; theta_s, theta_fc, soilthickness) = soil

    kh_profile_type = config.model.saturated_hydraulic_conductivity_profile

    if kh_profile_type == VerticalConductivityProfile.exponential
        (; kv_0, f) = soil.kv_profile
        # [m s⁻¹] = [-] * [m s⁻¹]
        kh_0 = khfrac .* kv_0
        kh_profile = KhExponential(kh_0, f)
    elseif kh_profile_type == VerticalConductivityProfile.exponential_constant
        (; z_exp) = soil.kv_profile
        (; kv_0, f) = soil.kv_profile.exponential
        # [m s⁻¹] = [-] * [m s⁻¹]
        kh_0 = khfrac .* kv_0
        exp_profile = KhExponential(kh_0, f)
        kh_profile = KhExponentialConstant(exp_profile, z_exp)
    elseif kh_profile_type == VerticalConductivityProfile.layered ||
           kh_profile_type == VerticalConductivityProfile.layered_exponential
        n_cells = length(khfrac)
        kh_profile = KhLayered(fill(MISSING_VALUE, n_cells))
    end
    specific_yield = @. lower_bound_drainable_porosity(theta_s, theta_fc)
    ssf_parameters = LateralSsfParameters(;
        kh_profile,
        khfrac,
        soilthickness,
        specific_yield,
        area,
        top = elevation,
    )
    return ssf_parameters
end

"Initialize lateral subsurface flow model variables"
function LateralSsfVariables(ssf::LateralSsfParameters, zi::Vector{Float64})
    n = length(zi)
    # [m³] = [-] * ([m] - [m]) * [m²]
    storage = @. ssf.specific_yield * (ssf.soilthickness - zi) * ssf.area
    # [m] = [m] - [m]
    head = ssf.top - zi
    variables = LateralSsfVariables(; n, zi, storage, head)
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
    zi = soil.variables.zi
    variables = LateralSsfVariables(parameters, zi)
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
    (; to_river) = subsurface_flow_model.variables
    (; river) = subsurface_flow_model.boundary_conditions
    if isnothing(river)
        average!(to_river, dt)
    else
        inds = domain.land_indices
        get_average(to_river)[inds] = -get_average(river.variables.flux_av)
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
        q_in_av,
        q,
        q_av,
        to_river,
        zi,
        head,
        exfiltwater,
        q_max,
        storage,
        q_net_bnds,
        q_net_av,
    ) = subsurface_flow_model.variables
    (; specific_yield, top, soilthickness, kh_profile) = subsurface_flow_model.parameters
    (; river) = subsurface_flow_model.boundary_conditions

    ns = length(order_of_subdomains)
    for k in 1:ns
        threaded_foreach(eachindex(order_of_subdomains[k]); basesize = 1) do i
            m = order_of_subdomains[k][i]
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
                    add_to_cumulative!(
                        to_river,
                        v,
                        sum_at(i -> q[i] * flow_fraction_to_river[i], upstream_nodes[n]),
                        dt,
                    )
                else
                    q_in[v] = sum_at(i -> q[i], upstream_nodes[n])
                end

                q[v], zi[v], _exfiltwater, netflux = kinematic_wave_ssf(
                    q_in[v],
                    q[v],
                    zi[v],
                    q_net_bnds[v],
                    slope[v],
                    specific_yield[v],
                    soilthickness[v],
                    dt,
                    flow_length[v],
                    flow_width[v],
                    q_max[v],
                    kh_profile,
                    soil_model,
                    v,
                )
                # [m³] += [m³ s⁻¹] * [s]
                add_to_cumulative!(q_in_av, v, q_in[v], dt)
                # [m³] += [m³ s⁻¹] * [s]
                add_to_cumulative!(q_av, v, q[v], dt)
                # [m s⁻¹] += [m s⁻¹]
                exfiltwater[v] += _exfiltwater
                # [m³] += [m s⁻¹] * [s]
                add_to_cumulative!(q_net_av, v, netflux * area[v], dt)
                # [m] = [m] - [m]
                head[v] = top[v] - zi[v]
                # [m³] = [-] * ([m] - [m]) * [m²]
                storage[v] = specific_yield[v] * (soilthickness[v] - zi[v]) * area[v]
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
    (; to_river) = subsurface_flow_model.variables
    (; adaptive) = subsurface_flow_model.timestepping

    zero!(to_river)
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
    (; zi) = subsurface_flow_model.variables
    (; specific_yield, kh_profile) = subsurface_flow_model.parameters
    (; flow_length, slope) = domain.parameters
    (; stable_timesteps, alpha_coefficient) = subsurface_flow_model.timestepping

    n = length(zi)
    stable_timesteps .= Inf
    k = 0
    for i in 1:n
        if zi[i] > 0.0
            k += 1
            # [m s⁻¹]
            c = ssf_celerity(zi[i], slope[i], specific_yield[i], kh_profile, i)
            # [s] = [m] / [m s⁻¹]
            stable_timesteps[k] = (flow_length[i] / c)
        end
    end

    dt_min = if k > 0
        minimum(@view(stable_timesteps[1:k]))
    else
        0.5
    end

    # [s] = [s] * [-]
    return dt_min * alpha_coefficient
end

get_flux_to_river(subsurface_flow_model::LateralSSFModel, inds::Vector{Int}) =
    get_average(subsurface_flow_model.variables.to_river)[inds]

# wrapper method
get_water_depth(subsurface_flow_model::LateralSSFModel) = subsurface_flow_model.variables.zi
