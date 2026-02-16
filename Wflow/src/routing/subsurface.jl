"Struct for storing lateral subsurface flow model variables"
@with_kw struct LateralSsfVariables
    n::Int
    zi::Vector{Float64}                                    # Pseudo-water table depth [m] (top of the saturated zone)
    head::Vector{Float64}                                  # Hydraulic head [m]
    exfiltwater::Vector{Float64} = fill(MISSING_VALUE, n)  # Exfiltration [m Δt⁻¹] (groundwater above surface level, saturated excess conditions)
    ssf::Vector{Float64} = fill(MISSING_VALUE, n)          # Subsurface flow [m³ d⁻¹]
    ssf_av::Vector{Float64} = fill(MISSING_VALUE, n)       # Average subsurface flow [m³ d⁻¹] for model timestep Δt
    ssfin::Vector{Float64} = fill(MISSING_VALUE, n)        # Inflow from upstream cells [m³ d⁻¹]
    ssfin_av::Vector{Float64} = fill(MISSING_VALUE, n)     # Average inflow from upstream cells [m³ d⁻¹] for model timestep Δt
    ssfmax::Vector{Float64} = fill(MISSING_VALUE, n)       # Maximum subsurface flow [m² d⁻¹]
    to_river::Vector{Float64} = fill(MISSING_VALUE, n)     # Part of subsurface flow [m³ d⁻¹] that flows to the river
    q_net::Vector{Float64} = fill(MISSING_VALUE, n)        # Net flow (boundaries) [m³ d⁻¹]
    net_flux::Vector{Float64} = fill(MISSING_VALUE, n)     # Net flux [m Δt⁻¹]
    storage::Vector{Float64}                               # Subsurface storage that can be released [m³]
end

"Struct for storing lateral subsurface flow model parameters"
@with_kw struct LateralSsfParameters{Kh}
    kh_profile::Kh                      # Horizontal hydraulic conductivity profile type [-]
    khfrac::Vector{Float64}             # A muliplication factor applied to vertical hydraulic conductivity `kv` [-]
    soilthickness::Vector{Float64}      # Soil thickness [m]
    specific_yield::Vector{Float64}     # Specific yield (theta_s - theta_fc) [-]
    area::Vector{Float64}               # Area of cell [m²]
    top::Vector{Float64}                # Top of subsurface flow layer [m]
end

"Lateral subsurface flow model"
@with_kw struct LateralSSF{Kh, B <: SubsurfaceFlowBC} <: AbstractSubsurfaceFlowModel
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
    indices::Vector{CartesianIndex{2}},
    soil::SbmSoilParameters,
    area::Vector{Float64},
)
    elevation = ncread(
        dataset,
        config,
        "land_surface__elevation";
        optional = false,
        sel = indices,
        type = Float64,
    )
    khfrac = ncread(
        dataset,
        config,
        "subsurface_water__horizontal_to_vertical_saturated_hydraulic_conductivity_ratio";
        optional = false,
        sel = indices,
        type = Float64,
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
    n = length(zi)
    storage = @. ssf.specific_yield * (ssf.soilthickness - zi) * ssf.area
    head = @. ssf.top - zi
    variables = LateralSsfVariables(; n, zi, storage, head)
    return variables
end

"Initialize lateral subsurface flow model"
function LateralSSF(dataset::NCDataset, config::Config, domain::Domain, soil::SbmSoilModel)
    (; land, river, drain) = domain
    (; indices) = land.network
    (; area) = domain.land.parameters
    n = length(indices)
    timestepping = init_kinematic_wave_timestepping(config, n; domain = "subsurface")
    parameters = LateralSsfParameters(dataset, config, indices, soil.parameters, area)
    zi = 0.001 * soil.variables.zi
    variables = LateralSsfVariables(parameters, zi)
    recharge = Recharge(; n)
    if config.model.river_subsurface_exchange_head_based__flag
        river = GwfRiver(dataset, config, river.network.indices)
    else
        river = nothing
    end
    if config.model.drain__flag
        drain = Drainage(dataset, config, drain.network.indices)
    else
        drain = nothing
    end
    boundary_conditions = SubsurfaceFlowBC(; recharge, river, drain)
    ssf = LateralSSF(; timestepping, boundary_conditions, parameters, variables)
    return ssf
end

function update_fluxes!(model::LateralSSF, domain::Domain, dt::Float64)
    for bc in get_boundaries(model.boundary_conditions)
        indices = get_boundary_index(bc, domain)
        flux!(bc, model, indices, dt)
    end
    return nothing
end

function kinwave_subsurface_update!(
    model::LateralSSF,
    soil::SbmSoilModel,
    domain::Domain,
    dt::Float64,
)
    (; order_of_subdomains, order_subdomain, subdomain_indices, upstream_nodes) =
        domain.land.network
    (; flow_length, flow_width, area, flow_fraction_to_river, slope) =
        domain.land.parameters

    (;
        ssfin,
        ssfin_av,
        ssf,
        ssf_av,
        to_river,
        zi,
        head,
        exfiltwater,
        ssfmax,
        storage,
        q_net,
        net_flux,
    ) = model.variables
    (; specific_yield, top, soilthickness, kh_profile) = model.parameters
    (; river) = model.boundary_conditions

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
                    ssfin[v] = sum_at(
                        i -> ssf[i] * (1.0 - flow_fraction_to_river[i]),
                        upstream_nodes[n],
                    )
                    to_river[v] +=
                        sum_at(i -> ssf[i] * flow_fraction_to_river[i], upstream_nodes[n]) *
                        dt
                else
                    ssfin[v] = sum_at(i -> ssf[i], upstream_nodes[n])
                end

                ssf[v], zi[v], _exfiltwater, netflux = kinematic_wave_ssf(
                    ssfin[v],
                    ssf[v],
                    zi[v],
                    q_net[v],
                    slope[v],
                    specific_yield[v],
                    soilthickness[v],
                    dt,
                    flow_length[v],
                    flow_width[v],
                    ssfmax[v],
                    kh_profile,
                    soil,
                    v,
                )
                ssfin_av[v] += ssfin[v] * dt
                ssf_av[v] += ssf[v] * dt
                exfiltwater[v] += _exfiltwater
                net_flux[v] += netflux
                head[v] = top[v] - zi[v]
                storage[v] = specific_yield[v] * (soilthickness[v] - zi[v]) * area[v]
            end
        end
    end
end

"""
Update lateral subsurface model for a single timestep `dt`. Timestepping within `dt` is
either with a fixed timestep `dt_fixed` or adaptive.
"""
function update!(model::LateralSSF, soil::SbmSoilModel, domain::Domain, dt::Float64)
    (; ssfin_av, ssf_av, to_river, exfiltwater, net_flux) = model.variables
    (; adaptive) = model.timestepping

    ssf_av .= 0.0
    to_river .= 0.0
    ssfin_av .= 0.0
    exfiltwater .= 0.0
    net_flux .= 0.0

    set_flux_vars_bc!(model)
    t = 0.0
    while t < dt
        model.variables.q_net .= 0.0
        dt_s = adaptive ? stable_timestep(model, domain.land) : model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        update_fluxes!(model, domain, dt_s)
        kinwave_subsurface_update!(model, soil, domain, dt_s)
        t += dt_s
    end
    ssf_av ./= dt
    to_river ./= dt
    ssfin_av ./= dt
    average_flux_vars_bc!(model, dt)
    return nothing
end

"""
Compute a stable timestep size for the kinematice wave method for a lateral subsurface flow
model using a nonlinear scheme (Chow et al., 1988).

A stable time step is computed for each vector element based on the Courant timestep size
criterion. Li et al. (1975) found that the nonlinear scheme is unconditonally stable and
that a wide range of dt/dx values can be used without loss of accuracy.
"""
function stable_timestep(model::LateralSSF, domain::DomainLand)
    (; zi) = model.variables
    (; specific_yield, kh_profile) = model.parameters
    (; flow_length, slope) = domain.parameters
    (; stable_timesteps, cfl) = model.timestepping

    n = length(zi)
    stable_timesteps .= Inf
    k = 0
    for i in 1:n
        if zi[i] > 0.0
            k += 1
            c = ssf_celerity(zi[i], slope[i], specific_yield[i], kh_profile, i)
            stable_timesteps[k] = (flow_length[i] / c)
        end
    end

    dt_min = if k > 0
        minimum(@view(stable_timesteps[1:k]))
    else
        0.5
    end

    dt_min = cfl * dt_min

    return dt_min
end

# wrapper methods
get_water_depth(model::LateralSSF) = model.variables.zi
get_exfiltwater(model::LateralSSF) = model.variables.exfiltwater

function get_flux_to_river(model::LateralSSF, inds::Vector{Int})
    (; river) = model.boundary_conditions
    dt = tosecond(BASETIMESTEP) # conversion to [m³ s⁻¹]
    flux = if isnothing(river)
        model.variables.to_river[inds] ./ dt
    else
        -river.variables.flux_av ./ dt
    end
    return flux
end

get_inflow(model::LateralSSF) = model.variables.ssfin_av
get_outflow(model::LateralSSF) = model.variables.ssf_av
get_storage(model::LateralSSF) = model.variables.storage
