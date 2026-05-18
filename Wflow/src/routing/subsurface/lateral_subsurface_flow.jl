"Struct for storing lateral subsurface flow model variables"
@with_kw struct LateralSsfVariables
    n_cells::Int
    zi::Vector{Float64}                                     # Pseudo-water table depth [m] (top of the saturated zone)
    head::Vector{Float64}                                   # Hydraulic head [m]
    exfiltwater::Vector{Float64} = fill(MISSING_VALUE, n_cells)   # Exfiltration [m Δt⁻¹] (groundwater above surface level, saturated excess conditions)
    q::Vector{Float64} = fill(MISSING_VALUE, n_cells)             # Subsurface flow [m³ d⁻¹]
    q_av::Vector{Float64} = fill(MISSING_VALUE, n_cells)          # Average subsurface flow [m³ d⁻¹] for model timestep Δt
    q_in::Vector{Float64} = fill(MISSING_VALUE, n_cells)          # Inflow from upstream cells [m³ d⁻¹]
    q_in_av::Vector{Float64} = fill(MISSING_VALUE, n_cells)       # Average inflow from upstream cells [m³ d⁻¹] for model timestep Δt
    q_max::Vector{Float64} = fill(MISSING_VALUE, n_cells)         # Maximum subsurface flow [m² d⁻¹]
    to_river::Vector{Float64} = fill(MISSING_VALUE, n_cells)      # Part of subsurface flow [m³ d⁻¹] that flows to the river
    q_net_bnds::Vector{Float64} = fill(MISSING_VALUE, n_cells)    # Net flow for boundaries subsurface flow [m³ d⁻¹]
    q_net_av::Vector{Float64} = fill(MISSING_VALUE, n_cells)      # Average net flow (total) [m³ d⁻¹]
    storage::Vector{Float64}                                # Subsurface storage that can be released [m³]
end

"Struct for storing lateral subsurface flow model parameters"
@with_kw struct LateralSsfParameters{Kh}
    kh_profile::Kh                      # Horizontal hydraulic conductivity profile type [-]
    khfrac::Vector{Float64}             # A multiplication factor applied to vertical hydraulic conductivity `kv` [-]
    soilthickness::Vector{Float64}      # Soil thickness [m]
    specific_yield::Vector{Float64}     # Specific yield (theta_s - theta_fc) [-]
    area::Vector{Float64}               # Area of cell [m²]
    top::Vector{Float64}                # Top of subsurface flow layer [m]
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
    land_indices_2d::Vector{CartesianIndex{2}},
    soil::SbmSoilParameters,
    area::Vector{Float64},
)
    elevation =
        ncread(dataset, config, "land_surface__elevation", Routing; sel = land_indices_2d)
    khfrac = ncread(
        dataset,
        config,
        "subsurface_water__horizontal_to_vertical_saturated_hydraulic_conductivity_ratio",
        Routing;
        sel = land_indices_2d,
    )

    (; theta_s, theta_fc, soilthickness) = soil
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
    n_cells = length(zi)
    storage = @. ssf.specific_yield * (ssf.soilthickness - zi) * ssf.area
    head = @. ssf.top - zi
    variables = LateralSsfVariables(; n_cells, zi, storage, head)
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
    (; land_indices_2d) = land.network
    (; area) = domain.land.parameters
    n_cells = length(land_indices_2d)
    timestepping = init_kinematic_wave_timestepping(config, n_cells; domain = "subsurface")
    parameters =
        LateralSsfParameters(dataset, config, land_indices_2d, soil.parameters, area)
    zi = 0.001 * soil.variables.zi
    variables = LateralSsfVariables(parameters, zi)
    recharge = RechargeModel(; n_cells)
    if config.model.river_subsurface_exchange_head_based__flag
        river = GwfRiverModel(dataset, config, river.network.river_indices_2d)
    else
        river = nothing
    end
    if config.model.drain__flag
        drain = DrainageModel(dataset, config, drain.network.drain_indices_2d)
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
        to_river ./= dt
    else
        to_river[domain.cell_indices_containing_river] .= -river.variables.flux_av
    end
    return nothing
end

function kinwave_subsurface_update!(
    subsurface_flow_model::LateralSSFModel,
    soil_model::SbmSoilModel,
    domain::Domain,
    dt::Float64,
)
    (; order_of_subdomains, order_subdomain, subdomain_global_order, upstream_nodes) =
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

    n_subdomain_sets = length(order_of_subdomains)
    for subdomain_set_idx in 1:n_subdomain_sets
        threaded_foreach(
            eachindex(order_of_subdomains[subdomain_set_idx]);
            basesize = 1,
        ) do in_subdomain_set_idx
            subdomain_idx = order_of_subdomains[subdomain_set_idx][in_subdomain_set_idx]
            for (land_global_traversion_idx, cell_idx) in
                zip(subdomain_global_order[subdomain_idx], order_subdomain[subdomain_idx])
                if isnothing(river)
                    # for a river cell without a reservoir part of the upstream subsurface flow
                    # goes to the river (flow_fraction_to_river) and part goes to the subsurface
                    # flow reservoir (1.0 - flow_fraction_to_river) upstream nodes with a
                    # reservoir are excluded
                    q_in[cell_idx] = sum_at(
                        cell_idx_other ->
                            q[cell_idx_other] *
                            (1.0 - flow_fraction_to_river[cell_idx_other]),
                        upstream_nodes[land_global_traversion_idx],
                    )
                    to_river[cell_idx] +=
                        sum_at(
                            cell_idx_other ->
                                q[cell_idx_other] * flow_fraction_to_river[cell_idx_other],
                            upstream_nodes[land_global_traversion_idx],
                        ) * dt
                else
                    q_in[cell_idx] = sum_at(
                        cell_idx_other -> q[cell_idx_other],
                        upstream_nodes[land_global_traversion_idx],
                    )
                end

                q[cell_idx], zi[cell_idx], _exfiltwater, netflux = kinematic_wave_ssf(
                    q_in[cell_idx],
                    q[cell_idx],
                    zi[cell_idx],
                    q_net_bnds[cell_idx],
                    slope[cell_idx],
                    specific_yield[cell_idx],
                    soilthickness[cell_idx],
                    dt,
                    flow_length[cell_idx],
                    flow_width[cell_idx],
                    q_max[cell_idx],
                    kh_profile,
                    soil_model,
                    cell_idx,
                )
                q_in_av[cell_idx] += q_in[cell_idx] * dt
                q_av[cell_idx] += q[cell_idx] * dt
                exfiltwater[cell_idx] += _exfiltwater
                q_net_av[cell_idx] += netflux * area[cell_idx]
                head[cell_idx] = top[cell_idx] - zi[cell_idx]
                storage[cell_idx] =
                    specific_yield[cell_idx] *
                    (soilthickness[cell_idx] - zi[cell_idx]) *
                    area[cell_idx]
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

    to_river .= 0.0
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
    (; zi, n_cells) = subsurface_flow_model.variables
    (; specific_yield, kh_profile) = subsurface_flow_model.parameters
    (; flow_length, slope) = domain.parameters
    (; stable_timesteps, alpha_coefficient) = subsurface_flow_model.timestepping

    stable_timesteps .= Inf
    stable_timestep_idx = 0
    for cell_idx in 1:n_cells
        if zi[cell_idx] > 0.0
            stable_timestep_idx += 1
            c = ssf_celerity(
                zi[cell_idx],
                slope[cell_idx],
                specific_yield[cell_idx],
                kh_profile,
                cell_idx,
            )
            stable_timesteps[stable_timestep_idx] = (flow_length[cell_idx] / c)
        end
    end

    dt_min = if stable_timestep_idx > 0
        minimum(@view(stable_timesteps[1:stable_timestep_idx]))
    else
        0.5
    end

    dt_min = alpha_coefficient * dt_min

    return dt_min
end

function get_flux_to_river(subsurface_flow_model::LateralSSFModel, inds::Vector{Int})
    dt = tosecond(BASETIMESTEP) # conversion to [m³ s⁻¹]
    flux = subsurface_flow_model.variables.to_river[inds] ./ dt
    return flux
end

# wrapper method
get_water_depth(subsurface_flow_model::LateralSSFModel) = subsurface_flow_model.variables.zi
