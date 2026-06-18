abstract type AbstractSubsurfaceFlowBC end

"""
    GroundwaterFlowModel

The upper boundary of an unconfined aquifer is the water table (the phreatic
surface).

Specific yield (or drainable porosity) represents the volumetric fraction the
aquifer will yield when all water drains and the pore volume is filled by air
instead. Specific yield will vary roughly between 0.05 (clay) and 0.45 (peat)
(Johnson, 1967).

The vertically averaged governing equation of an unconfined, inhomogeneous and
isotropic aquifer in one dimension can be written as:

S * ∂ϕ / ∂t = ∂ / ∂x * (kH * ∂ϕ / ∂x) + Q

with:
* S: storativity (or storage coefficient)
* ϕ: hydraulic head (groundwater level)
* t: time
* hydraulic_conductivity: conductivity
* H: H the (saturated) aquifer height: groundwater level - aquifer bottom
  elevation
* η: elevation of aquifer bottom
* Q: fluxes from boundary conditions (e.g. recharge or abstraction)

The simplest finite difference formulation is forward in time, central in space,
and can be written as:

Sᵢ * (ϕᵢᵗ⁺¹ - ϕᵢᵗ) / Δt = -Cᵢ₋₁ * (ϕᵢ₋₁ - ϕᵢ) -Cᵢ * (ϕᵢ₊₁ - ϕᵢ) + Qᵢ

with:
* ᵢ as cell index
* ᵗ as time index
* Δt as step size
* Cᵢ₋₁ as the intercell conductance between cell i-1 and i
* Cᵢ as the intercell conductance between cell i and i+1

Conductance is defined as:

C = kH * w / l

with:
* w the width of the cell to cell connection
* l the length of the cell to cell connection

hydraulic_conductivity and H may both vary in space; intercell conductance is therefore an average
using the properties of two cells. See the documentation below.

There is only one unknown, ϕᵢᵗ⁺¹. Reshuffling terms:

ϕᵢᵗ⁺¹ = ϕᵢᵗ + (Cᵢ₋₁ * (ϕᵢ - ϕᵢ₋₁) + Cᵢ * (ϕᵢ₊₁ - ϕᵢ) + Qᵢ) * Δt / Sᵢ

This can be generalized to two dimensions, for both regular and irregular cell
connectivity.

See this paper for more details:
     Chu, W. S., & Willis, R. (1984). An explicit finite difference model for
     unconfined aquifers. Groundwater, 22(6), 728-734.

Boundary conditions can be classified into three categories:
* specified head (Dirichlet)
* specified flux (Neumann)
* head-dependent flux (Robin)

Neumann and Robin conditions are implemented by adding to or subtracting from a
net (lumped) cell flux. Dirichlet conditions are special cased, since they
cannot (easily) be implemented via the flux, but the head is set directly
instead.
"""

@with_kw struct GroundwaterFlowVariables
    n::Int
    # hydraulic head [m]
    head::Vector{Float64}
    # conductance [m² s⁻¹]
    conductance::Vector{Float64}
    # total storage of water that can be released [m³]
    storage::Vector{Float64}
    # net flow (groundwater and boundaries) [m³ s⁻¹]
    q_net::Vector{Float64} = zeros(n)
    # cumulative net flow (groundwater and boundaries) [m³]
    q_net_cumulative::Vector{Float64} = zeros(n)
    # average net flow (groundwater and boundaries) [m³ s⁻¹]
    q_net_average::Vector{Float64} = zeros(n)
    # net flow boundaries [m³ s⁻¹]
    q_net_bnds::Vector{Float64} = zeros(n)
    # cumulative groundwater (lateral) inflow for model timestep dt [m³]
    q_in_cumulative::Vector{Float64} = zeros(n)
    # average groundwater (lateral) inflow for model timestep dt [m³ s⁻¹]
    q_in_average::Vector{Float64} = zeros(n)
    # cumulative groundwater (lateral) outflow for model timestep dt [m³]
    q_cumulative::Vector{Float64} = zeros(n)
    # average groundwater (lateral) outflow for model timestep dt [m³ s⁻¹]
    q_average::Vector{Float64} = zeros(n)
    # Cumulative exfiltration [m] (groundwater above surface level, saturated excess conditions)
    exfiltwater_cumulative::Vector{Float64} = zeros(n)
    # Average exfiltration [m s⁻¹] (groundwater above surface level, saturated excess conditions)
    exfiltwater_average::Vector{Float64} = zeros(n)
end

@with_kw struct GroundwaterFlowParameters
    # reference horizontal conductivity [m s⁻¹]
    hydraulic_conductivity::Vector{Float64}
    # top of groundwater layer [m]
    top::Vector{Float64}
    # bottom of groundwater layer [m]
    bottom::Vector{Float64}
    # area of cell [m²]
    area::Vector{Float64}
    # specific yield (theta_s - theta_fc) [m m⁻¹]
    specific_yield::Vector{Float64}
    # factor controlling the reduction of reference horizontal conductivity [-]
    hydraulic_conductivity_scale_parameter::Vector{Float64}
    # Unconfined aquifer conductance is computed with degree of saturation (only when
    # conductivity_profile is set to "exponential")
end

function GroundwaterFlowParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    top::Vector{Float64},
    bottom::Vector{Float64},
    area::Vector{Float64},
    specific_yield::Vector{Float64},
)
    hydraulic_conductivity = ncread(
        dataset,
        config,
        "subsurface_surface_water__horizontal_saturated_hydraulic_conductivity",
        Routing;
        sel = indices,
    )
    if config.model.conductivity_profile == GwfConductivityProfileType.exponential
        hydraulic_conductivity_scale_parameter = ncread(
            dataset,
            config,
            "subsurface__horizontal_saturated_hydraulic_conductivity_scale_parameter",
            Routing;
            sel = indices,
        )
    else
        hydraulic_conductivity_scale_parameter = Float64[]
    end
    parameters = GroundwaterFlowParameters(;
        hydraulic_conductivity,
        top,
        bottom,
        area,
        specific_yield,
        hydraulic_conductivity_scale_parameter,
    )
    return parameters
end

@with_kw struct ConstantHeadVariables
    # [m]
    head::Vector{Float64}
end

@with_kw struct ConstantHead
    variables::ConstantHeadVariables
    # [-]
    index::Vector{Int}
end

function ConstantHead(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    constanthead = ncread(
        dataset,
        config,
        "model_constant_boundary_condition__hydraulic_head",
        Routing;
        sel = indices,
    )
    n = length(indices)
    index_constanthead = filter(i -> !isequal(constanthead[i], MISSING_VALUE), 1:n)
    head = constanthead[index_constanthead]
    variables = ConstantHeadVariables(head)
    constant_head = ConstantHead(; variables, index = index_constanthead)
    return constant_head
end

@kwdef struct SubsurfaceFlowBC{
    Re <: Union{Nothing, AbstractSubsurfaceFlowBC},
    Ri <: Union{Nothing, AbstractSubsurfaceFlowBC},
    D <: Union{Nothing, AbstractSubsurfaceFlowBC},
    W <: Union{Nothing, AbstractSubsurfaceFlowBC},
}
    recharge::Re = nothing
    river::Ri = nothing
    drain::D = nothing
    well::W = nothing
end

get_boundaries(boundary_conditions::SubsurfaceFlowBC) = (
    boundary_conditions.recharge,
    boundary_conditions.river,
    boundary_conditions.drain,
    boundary_conditions.well,
)

@kwdef struct GroundwaterFlowModel{B <: SubsurfaceFlowBC} <: AbstractSubsurfaceFlowModel
    timestepping::TimeStepping
    parameters::GroundwaterFlowParameters
    variables::GroundwaterFlowVariables
    connectivity::Connectivity
    constanthead::ConstantHead
    boundary_conditions::B = SubsurfaceFlowBC()
    function GroundwaterFlowModel(
        timestepping::TimeStepping,
        parameters::GroundwaterFlowParameters,
        variables::GroundwaterFlowVariables,
        connectivity::Connectivity,
        constanthead::ConstantHead,
        boundary_conditions::B,
    ) where {B <: SubsurfaceFlowBC}
        initialize_conductance!(parameters, variables, connectivity)
        new{B}(
            timestepping,
            parameters,
            variables,
            connectivity,
            constanthead,
            boundary_conditions,
        )
    end
end

storativity(A::GroundwaterFlowModel) = A.parameters.specific_yield

"Initialize groundwater flow model"
function GroundwaterFlowModel(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil::SbmSoilModel,
)
    (; land, river, drain) = domain

    (; indices, reverse_indices) = land.network
    (; x_length, y_length, area) = land.parameters

    n_cells = length(indices)

    elevation = ncread(dataset, config, "land_surface__elevation", Routing; sel = indices)

    # unconfined aquifer
    if config.model.constanthead__flag
        constanthead = ConstantHead(dataset, config, indices)
    else
        variables = ConstantHeadVariables(; head = Float64[])
        constanthead = ConstantHead(; variables, index = Int64[])
    end

    connectivity = Connectivity(indices, reverse_indices, x_length, y_length)

    # cold state for groundwater head based on water table depth zi
    initial_head = elevation .- soil.variables.water_table_depth
    initial_head[river.network.land_indices] = elevation[river.network.land_indices]
    if config.model.constanthead__flag
        initial_head[constanthead.index] = constanthead.variables.head
    end
    # reset soil (cold) state and related variables based on initial_head (river cells and constanthead)
    if config.model.cold_start__flag
        (;
            water_table_depth,
            saturated_water_depth,
            unsaturated_store_capacity,
            unsaturated_layer_thickness,
            n_unsatlayers,
            total_soil_water_storage,
        ) = soil.variables
        (;
            theta_s,
            theta_r,
            soil_thickness,
            soil_water_capacity,
            cumulative_layer_depth,
            actual_layer_thickness,
        ) = soil.parameters

        @. water_table_depth = elevation - min(elevation, initial_head)
        @. saturated_water_depth =
            (soil_thickness - water_table_depth) * (theta_s - theta_r)
        @. unsaturated_store_capacity = soil_water_capacity - saturated_water_depth
        @. unsaturated_layer_thickness = set_layerthickness(
            water_table_depth,
            cumulative_layer_depth,
            actual_layer_thickness,
        )
        @. n_unsatlayers = number_of_active_layers.(unsaturated_layer_thickness)
        @. total_soil_water_storage = saturated_water_depth
    end

    bottom = elevation - soil.parameters.soil_thickness
    specific_yield =
        @. lower_bound_drainable_porosity(soil.parameters.theta_s, soil.parameters.theta_fc)
    conductance = zeros(connectivity.nconnection)
    parameters = GroundwaterFlowParameters(
        dataset,
        config,
        indices,
        elevation,
        bottom,
        area,
        specific_yield,
    )
    storage = @. (min(elevation, initial_head) - bottom) * area * parameters.specific_yield
    n = length(storage)
    variables = GroundwaterFlowVariables(; n, head = initial_head, conductance, storage)

    # river boundary of unconfined aquifer
    gwf_river_model = GwfRiverModel(dataset, config, river.network.indices)

    # recharge boundary of unconfined aquifer
    recharge_model = RechargeModel(; n = n_cells)

    # drain boundary of unconfined aquifer (optional)
    if config.model.drain__flag
        drainage_model = DrainageModel(dataset, config, drain.network.indices)
        boundary_conditions = SubsurfaceFlowBC(;
            recharge = recharge_model,
            river = gwf_river_model,
            drain = drainage_model,
        )
    else
        boundary_conditions =
            SubsurfaceFlowBC(; recharge = recharge_model, river = gwf_river_model)
    end

    alpha_coefficient = config.model.subsurface_water_flow__alpha_coefficient
    @info "Numerical stability coefficient for groundwater flow `alpha`: `$alpha_coefficient`."
    timestepping = TimeStepping(; alpha_coefficient)

    groundwaterflow_model = GroundwaterFlowModel(;
        timestepping,
        parameters,
        variables,
        connectivity,
        constanthead,
        boundary_conditions,
    )
    return groundwaterflow_model
end

"""
    harmonicmean_conductance(kH1, kH2, l1, l2, width)

The harmonic mean is the exact interblock transmissivity for steady-state
one-dimensional flow with no recharge if the transmissivity is assumed to be
spatially uniform over each finite difference block, changing abruptly at the
block interface.

Refer to:

    Goode, D. J., & Appel, C. A. (1992). Finite-Difference Interblock Transmissivity
    for Unconﬁned Aquifers and for Aquifers having Smoothly Varying Transmissivity.
    Water-resources investigations report, 92, 4124.
"""
function harmonicmean_conductance(kH1, kH2, l1, l2, width)
    if (kH1 * kH2) > 0.0
        return width * kH1 * kH2 / (kH1 * l2 + kH2 * l1)
    else
        return 0.0
    end
end

function saturated_thickness(gwf::GroundwaterFlowModel, index::Int)
    return min(gwf.parameters.top[index], gwf.variables.head[index]) -
           gwf.parameters.bottom[index]
end

function saturated_thickness(gwf::GroundwaterFlowModel)
    @. min(gwf.parameters.top, gwf.variables.head) - gwf.parameters.bottom
end

"""
    horizontal_conductance(i::Int, j::int, nzi::Int, parameters::GroundwaterFlowParameters, connectivity::Connectivity)

Compute fully saturated horizontal conductance for a single connection between two cells
(indexed with `i` and `j`). Geometry characteristics are taken from the `connectivity`
struct , using the non-zero index (nzi) of its CSC data structure.
"""
function horizontal_conductance(
    i::Int,
    j::Int,
    nzi::Int,
    parameters::GroundwaterFlowParameters,
    connectivity::Connectivity,
)
    k1 = parameters.hydraulic_conductivity[i]
    k2 = parameters.hydraulic_conductivity[j]
    H1 = parameters.top[i] - parameters.bottom[i]
    H2 = parameters.top[j] - parameters.bottom[j]
    length1 = connectivity.length1[nzi]
    length2 = connectivity.length2[nzi]
    width = connectivity.width[nzi]
    kH1 = k1 * H1
    kH2 = k2 * H2
    return harmonicmean_conductance(kH1, kH2, length1, length2, width)
end

"""
    initialize_conductance!(parameters::GroundwaterFlowParameters, variables::GroundwaterFlowVariables, connectivity::Connectivity)

Conductance is computed per timestep by multiplying by degree of saturation [0.0 - 1.0].
"""
function initialize_conductance!(
    parameters::GroundwaterFlowParameters,
    variables::GroundwaterFlowVariables,
    connectivity::Connectivity,
)
    for i in 1:(connectivity.ncell)
        # Loop over connections for cell j
        for nzi in connections(connectivity, i)
            j = connectivity.rowval[nzi]
            variables.conductance[nzi] =
                horizontal_conductance(i, j, nzi, parameters, connectivity)
        end
    end
end

"""
    conductance(gwf::GroundwaterFlowModel, i, j, nzi, conductivity_profile::GwfConductivityProfileType.T, connectivity::Connectivity)

This computes the conductance for an unconfined aquifer using the "upstream
saturated fraction" as the MODFLOW documentation calls it. In this approach, the
saturated thickness of a cell-to-cell is approximated using the cell with the
highest head. This results in a consistent overestimation of the saturated
thickness, but it avoids complexities related with cell drying and rewetting,
such as having to define a "wetting threshold" or a "wetting factor".

See the documentation for MODFLOW-NWT or MODFLOW6 for more background:
Niswonger, R.G., Panday, Sorab, and Ibaraki, Motomu, 2011, MODFLOW-NWT, A Newton
formulation for MODFLOW-2005: U.S. Geological Survey Techniques and Methods
6-A37, 44 p.

Langevin, C.D., Hughes, J.D., Banta, E.R., Niswonger, R.G., Panday, Sorab, and
Provost, A.M., 2017, Documentation for the MODFLOW 6 Groundwater Flow Model:
U.S. Geological Survey Techniques and Methods, book 6, chap. A55, 197 p.,
https://doi.org/10.3133/tm6A55.

For background on drying and rewetting, see:
McDonald, M.G., Harbaugh, A.W., Orr, B.R., and Ackerman, D.J., 1991, A method of
converting no-flow cells to variable-head cells for the U.S. Geological Survey
modular finite-difference groundwater flow model: U.S. Geological Survey
Open-File Report 91-536, 99 p
"""
function conductance(
    gwf::GroundwaterFlowModel,
    i,
    j,
    nzi,
    conductivity_profile::GwfConductivityProfileType.T,
)
    if conductivity_profile == GwfConductivityProfileType.exponential
        # Extract required variables
        zi1 = gwf.parameters.top[i] - gwf.variables.head[i]
        zi2 = gwf.parameters.top[j] - gwf.variables.head[j]
        thickness1 = gwf.parameters.top[i] - gwf.parameters.bottom[i]
        thickness2 = gwf.parameters.top[j] - gwf.parameters.bottom[j]
        # calculate conductivity values corrected for depth of water table
        kH1 =
            (
                gwf.parameters.hydraulic_conductivity[i] /
                gwf.parameters.hydraulic_conductivity_scale_parameter[i]
            ) * (
                exp(-gwf.parameters.hydraulic_conductivity_scale_parameter[i] * zi1) -
                exp(-gwf.parameters.hydraulic_conductivity_scale_parameter[i] * thickness1)
            )
        kH2 =
            (
                gwf.parameters.hydraulic_conductivity[j] /
                gwf.parameters.hydraulic_conductivity_scale_parameter[j]
            ) * (
                exp(-gwf.parameters.hydraulic_conductivity_scale_parameter[j] * zi2) -
                exp(-gwf.parameters.hydraulic_conductivity_scale_parameter[j] * thickness2)
            )
        return harmonicmean_conductance(
            kH1,
            kH2,
            gwf.connectivity.length1[nzi],
            gwf.connectivity.length2[nzi],
            gwf.connectivity.width[nzi],
        )
    elseif conductivity_profile == GwfConductivityProfileType.uniform
        head_i = gwf.variables.head[i]
        head_j = gwf.variables.head[j]
        if head_i >= head_j
            saturation =
                saturated_thickness(gwf, i) /
                (gwf.parameters.top[i] - gwf.parameters.bottom[i])
        else
            saturation =
                saturated_thickness(gwf, j) /
                (gwf.parameters.top[j] - gwf.parameters.bottom[j])
        end
        return saturation * gwf.variables.conductance[nzi]
    end
end

function flux!(
    gwf::GroundwaterFlowModel,
    conductivity_profile::GwfConductivityProfileType.T,
    dt::Float64,
)
    for i in 1:(gwf.connectivity.ncell)
        # Loop over connections for cell j
        for nzi in connections(gwf.connectivity, i)
            # connection from i -> j
            j = gwf.connectivity.rowval[nzi]
            delta_head = gwf.variables.head[i] - gwf.variables.head[j]
            cond = conductance(gwf, i, j, nzi, conductivity_profile)
            flow = cond * delta_head
            gwf.variables.q_net[i] -= flow
            if flow > 0.0
                gwf.variables.q_cumulative[i] += flow * dt
            else
                gwf.variables.q_in_cumulative[i] -= flow * dt
            end
        end
    end
    return nothing
end

"""
    stable_timestep(gwf::GroundwaterFlowModel, conductivity_profile::GwfConductivityProfileType.T, alpha_coefficient::Float64)

Compute a stable timestep size given the forward-in-time, central in space scheme.
The following criterion can be found in Chu & Willis (1984):

Δt * hydraulic_conductivity * H / (Δx * Δy * S) <= α.
"""
function stable_timestep(
    gwf::GroundwaterFlowModel,
    conductivity_profile::GwfConductivityProfileType.T,
    alpha_coefficient::Float64,
)
    dt_min = Inf
    for i in eachindex(gwf.variables.head)
        if conductivity_profile == GwfConductivityProfileType.exponential
            water_table_depth = gwf.parameters.top[i] - gwf.variables.head[i]
            thickness = gwf.parameters.top[i] - gwf.parameters.bottom[i]
            value =
                (
                    gwf.parameters.hydraulic_conductivity[i] /
                    gwf.parameters.hydraulic_conductivity_scale_parameter[i]
                ) * (
                    exp(
                        -gwf.parameters.hydraulic_conductivity_scale_parameter[i] *
                        water_table_depth,
                    ) - exp(
                        -gwf.parameters.hydraulic_conductivity_scale_parameter[i] *
                        thickness,
                    )
                )
        elseif conductivity_profile == GwfConductivityProfileType.uniform
            value = gwf.parameters.hydraulic_conductivity[i] * saturated_thickness(gwf, i)
        end
        dt = gwf.parameters.area[i] * storativity(gwf)[i] / value
        dt_min = dt < dt_min ? dt : dt_min
    end
    return dt_min * alpha_coefficient
end

minimum_head(gwf::GroundwaterFlowModel) = max.(gwf.variables.head, gwf.parameters.bottom)
maximum_head(gwf::GroundwaterFlowModel) = min.(gwf.variables.head, gwf.parameters.top)

function update_fluxes!(
    gwf::GroundwaterFlowModel,
    domain::Domain,
    conductivity_profile::GwfConductivityProfileType.T,
    dt::Float64,
)
    flux!(gwf, conductivity_profile, dt)
    for bc in get_boundaries(gwf.boundary_conditions)
        indices = get_boundary_index(bc, domain)
        flux!(bc, gwf, indices, dt)
    end
    gwf.variables.q_net .+= gwf.variables.q_net_bnds
    @. gwf.variables.q_net_cumulative += gwf.variables.q_net * dt
    return nothing
end

function update_head!(gwf::GroundwaterFlowModel, soil::SbmSoilModel, dt::Float64)
    (; head, exfiltwater_cumulative, q_net) = gwf.variables
    (; area, specific_yield) = gwf.parameters

    for i in eachindex(head)
        net_flux = q_net[i] / area[i]
        dh, exfilt = water_table_change(soil, net_flux, specific_yield[i], i, dt)
        head[i] += dh
        exfiltwater_cumulative[i] += exfilt * dt
    end
    # Set constant head (dirichlet) boundaries
    gwf.variables.head[gwf.constanthead.index] .= gwf.constanthead.variables.head
    # Make sure no heads ends up below an unconfined aquifer bottom
    gwf.variables.head .= minimum_head(gwf)
    # Adjust exfiltration rate for constant head boundaries
    exfiltwater_cumulative[gwf.constanthead.index] .= 0.0
    gwf.variables.storage .=
        saturated_thickness(gwf) .* gwf.parameters.area .* storativity(gwf)
    return nothing
end

function set_flux_vars_bc!(gwf::AbstractSubsurfaceFlowModel)
    for bc in get_boundaries(gwf.boundary_conditions)
        !isnothing(bc) && (bc.variables.flux_cumulative .= 0.0)
    end
    return nothing
end

function set_flux_vars!(gwf::AbstractSubsurfaceFlowModel)
    set_flux_vars_bc!(gwf)
    gwf.variables.exfiltwater_cumulative .= 0.0
    gwf.variables.q_in_cumulative .= 0.0
    gwf.variables.q_cumulative .= 0.0
    gwf.variables.q_net_cumulative .= 0.0
    return nothing
end

function average_flux_vars_bc!(gwf::AbstractSubsurfaceFlowModel, dt::Float64)
    for bc in get_boundaries(gwf.boundary_conditions)
        if !isnothing(bc)
            @. bc.variables.flux_average = bc.variables.flux_cumulative / dt
        end
    end
end

function average_flux_vars!(gwf::AbstractSubsurfaceFlowModel, dt::Float64)
    average_flux_vars_bc!(gwf, dt)

    @. gwf.variables.q_in_average = gwf.variables.q_in_cumulative / dt
    @. gwf.variables.q_average = gwf.variables.q_cumulative / dt
    @. gwf.variables.q_net_average = gwf.variables.q_net_cumulative / dt
    @. gwf.variables.exfiltwater_average = gwf.variables.exfiltwater_cumulative / dt
    return nothing
end

function update_subsurface_flow_model!(
    gwf_model::GroundwaterFlowModel,
    soil_model::SbmSoilModel,
    domain::Domain,
    dt::Float64,
    conductivity_profile::GwfConductivityProfileType.T,
)
    (; alpha_coefficient) = gwf_model.timestepping

    set_flux_vars!(gwf_model)
    t = 0.0
    while t < dt
        gwf_model.variables.q_net .= 0.0
        gwf_model.variables.q_net_bnds .= 0.0
        dt_s = stable_timestep(gwf_model, conductivity_profile, alpha_coefficient)
        dt_s = check_timestepsize(dt_s, t, dt)
        update_fluxes!(gwf_model, domain, conductivity_profile, dt_s)
        update_head!(gwf_model, soil_model, dt_s)
        update_ustorelayerdepth!(soil_model, gwf_model)
        t += dt_s
    end
    average_flux_vars!(gwf_model, dt)
    return nothing
end

get_water_depth(gwf_model::GroundwaterFlowModel) =
    gwf_model.parameters.top .- gwf_model.variables.head

get_flux_to_river(subsurface_flow_model::GroundwaterFlowModel, inds::Vector{Int}) =
    -subsurface_flow_model.boundary_conditions.river.variables.flux_average

function sum_boundary_fluxes(
    gwf_model::AbstractSubsurfaceFlowModel,
    domain::Domain;
    exclude = nothing,
)
    (; boundary_conditions) = gwf_model
    n = length(gwf_model.variables.storage)
    flux_in = zeros(n)
    flux_out = zeros(n)
    for bc in get_boundaries(boundary_conditions)
        isnothing(bc) && continue
        typeof(bc) == exclude && continue
        flux_av_average = bc.variables.flux_average
        indices = get_boundary_index(bc, domain)
        for (i, index) in enumerate(indices)
            flux = flux_av_average[i]
            if flux > 0.0
                flux_in[index] += flux
            else
                flux_out[index] -= flux
            end
        end
    end
    return flux_in, flux_out
end
