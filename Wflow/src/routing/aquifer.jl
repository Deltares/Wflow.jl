abstract type AbstractSubsurfaceFlowBC end

"""
    GroundwaterFlow

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
* k: conductivity
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

k and H may both vary in space; intercell conductance is therefore an average
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
    head::Vector{Float64}                          # hydraulic head [m]
    conductance::Vector{Float64}                   # conductance [m² d⁻¹]
    storage::Vector{Float64}                       # total storage of water that can be released [m³]
    q_net::Vector{Float64} = zeros(n)              # net flow (groundwater and boundaries) [m³ d⁻¹]
    q_in_av::Vector{Float64} = zeros(n)            # average groundwater (lateral) inflow for model timestep Δt [m³ d⁻¹]
    q_out_av::Vector{Float64} = zeros(n)           # average groundwater (lateral) outflow for model timestep Δt [m³ d⁻¹]
    exfiltwater::Vector{Float64} = zeros(n)        # Exfiltration [m Δt⁻¹] (groundwater above surface level, saturated excess conditions)
end

@with_kw struct GroundwaterFlowParameters
    k::Vector{Float64}                  # reference horizontal conductivity [m d⁻¹]
    top::Vector{Float64}                # top of groundwater layer [m]
    bottom::Vector{Float64}             # bottom of groundwater layer [m]
    area::Vector{Float64}               # area of cell [m²]
    specific_yield::Vector{Float64}     # specific yield (theta_s - theta_fc) [m m⁻¹]
    f::Vector{Float64}                  # factor controlling the reduction of reference horizontal conductivity [-]
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
    k = ncread(
        dataset,
        config,
        "subsurface_surface_water__horizontal_saturated_hydraulic_conductivity";
        optional = false,
        sel = indices,
        type = Float64,
    )
    if config.model.conductivity_profile == GwfConductivityProfileType.exponential
        f = ncread(
            dataset,
            config,
            "subsurface__horizontal_saturated_hydraulic_conductivity_scale_parameter";
            optional = false,
            sel = indices,
            type = Float64,
        )
    else
        f = Float64[]
    end
    parameters = GroundwaterFlowParameters(; k, top, bottom, area, specific_yield, f)
    return parameters
end

@with_kw struct ConstantHeadVariables
    head::Vector{Float64} # [m]
end

@with_kw struct ConstantHead
    variables::ConstantHeadVariables
    index::Vector{Int} # [-]
end

function ConstantHead(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    constanthead = ncread(
        dataset,
        config,
        "model_constant_boundary_condition__hydraulic_head";
        optional = false,
        sel = indices,
        type = Float64,
        fill = MISSING_VALUE,
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

@kwdef struct GroundwaterFlow{B <: SubsurfaceFlowBC} <: AbstractSubsurfaceFlowModel
    timestepping::TimeStepping
    parameters::GroundwaterFlowParameters
    variables::GroundwaterFlowVariables
    connectivity::Connectivity
    constanthead::ConstantHead
    boundary_conditions::B = SubsurfaceFlowBC()
    function GroundwaterFlow(
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

storativity(A::GroundwaterFlow) = A.parameters.specific_yield

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

function saturated_thickness(gwf::GroundwaterFlow, index::Int)
    return min(gwf.parameters.top[index], gwf.variables.head[index]) -
           gwf.parameters.bottom[index]
end

function saturated_thickness(gwf::GroundwaterFlow)
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
    k1 = parameters.k[i]
    k2 = parameters.k[j]
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
    conductance(gwf::GroundwaterFlow, i, j, nzi, conductivity_profile::GwfConductivityProfileType.T, connectivity::Connectivity)

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
    gwf::GroundwaterFlow,
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
        k1 =
            (gwf.parameters.k[i] / gwf.parameters.f[i]) *
            (exp(-gwf.parameters.f[i] * zi1) - exp(-gwf.parameters.f[i] * thickness1))
        k2 =
            (gwf.parameters.k[j] / gwf.parameters.f[j]) *
            (exp(-gwf.parameters.f[j] * zi2) - exp(-gwf.parameters.f[j] * thickness2))
        return harmonicmean_conductance(
            k1,
            k2,
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
    gwf::GroundwaterFlow,
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
                gwf.variables.q_out_av[i] += flow * dt
            else
                gwf.variables.q_in_av[i] -= flow * dt
            end
        end
    end
    return nothing
end

"""
    stable_timestep(gwf::GroundwaterFlow, conductivity_profile::GwfConductivityProfileType.T, cfl::Float64)

Compute a stable timestep size given the forward-in-time, central in space scheme.
The following criterion can be found in Chu & Willis (1984):

Δt * k * H / (Δx * Δy * S) <= cfl.
"""
function stable_timestep(
    gwf::GroundwaterFlow,
    conductivity_profile::GwfConductivityProfileType.T,
    cfl::Float64,
)
    dt_min = Inf
    for i in eachindex(gwf.variables.head)
        if conductivity_profile == GwfConductivityProfileType.exponential
            zi = gwf.parameters.top[i] - gwf.variables.head[i]
            thickness = gwf.parameters.top[i] - gwf.parameters.bottom[i]
            value =
                (gwf.parameters.k[i] / gwf.parameters.f[i]) *
                (exp(-gwf.parameters.f[i] * zi) - exp(-gwf.parameters.f[i] * thickness))
        elseif conductivity_profile == GwfConductivityProfileType.uniform
            value = gwf.parameters.k[i] * saturated_thickness(gwf, i)
        end

        dt = gwf.parameters.area[i] * storativity(gwf)[i] / value
        dt_min = dt < dt_min ? dt : dt_min
    end
    dt_min = cfl * dt_min
    return dt_min
end

minimum_head(gwf::GroundwaterFlow) = max.(gwf.variables.head, gwf.parameters.bottom)
maximum_head(gwf::GroundwaterFlow) = min.(gwf.variables.head, gwf.parameters.top)

function update_fluxes!(
    gwf::GroundwaterFlow,
    domain::Domain,
    conductivity_profile::GwfConductivityProfileType.T,
    dt::Float64,
)
    flux!(gwf, conductivity_profile, dt)
    for bc in get_boundaries(gwf.boundary_conditions)
        indices = get_boundary_index(bc, domain)
        flux!(bc, gwf, indices, dt)
    end
    return nothing
end

function update_head!(gwf::GroundwaterFlow, soil::SbmSoilModel, dt::Float64)
    (; head, exfiltwater, q_net) = gwf.variables
    (; area, specific_yield) = gwf.parameters

    for i in eachindex(head)
        net_flux = q_net[i] / area[i] * dt
        dh, exfilt = water_table_change(soil, net_flux, specific_yield[i], i)
        head[i] += dh
        exfiltwater[i] += exfilt
    end
    # Set constant head (dirichlet) boundaries
    gwf.variables.head[gwf.constanthead.index] .= gwf.constanthead.variables.head
    # Make sure no heads ends up below an unconfined aquifer bottom
    gwf.variables.head .= minimum_head(gwf)
    # Compute exfiltration rate and make sure head is not above surface for unconfined aquifer
    gwf.variables.exfiltwater .+=
        (gwf.variables.head .- maximum_head(gwf)) .* storativity(gwf)
    gwf.variables.head .= maximum_head(gwf)
    # Adjust exfiltration rate for constant head boundaries
    gwf.variables.exfiltwater[gwf.constanthead.index] .= 0.0
    gwf.variables.storage .=
        saturated_thickness(gwf) .* gwf.parameters.area .* storativity(gwf)
    return nothing
end

function set_flux_vars_bc!(gwf::AbstractSubsurfaceFlowModel)
    for bc in get_boundaries(gwf.boundary_conditions)
        !isnothing(bc) && (bc.variables.flux_av .= 0.0)
    end
    return nothing
end

function set_flux_vars!(gwf::GroundwaterFlow)
    set_flux_vars_bc!(gwf)
    gwf.variables.exfiltwater .= 0.0
    gwf.variables.q_in_av .= 0.0
    gwf.variables.q_out_av .= 0.0
    return nothing
end

function average_flux_vars_bc!(gwf::AbstractSubsurfaceFlowModel, dt::Float64)
    for bc in get_boundaries(gwf.boundary_conditions)
        !isnothing(bc) && (bc.variables.flux_av ./= dt)
    end
end

function average_flux_vars!(gwf::GroundwaterFlow, dt::Float64)
    average_flux_vars_bc!(gwf, dt)
    gwf.variables.q_in_av ./= dt
    gwf.variables.q_out_av ./= dt
    return nothing
end

function update!(
    gwf::GroundwaterFlow,
    soil::SbmSoilModel,
    domain::Domain,
    dt::Float64,
    conductivity_profile::GwfConductivityProfileType.T,
)
    (; cfl) = gwf.timestepping

    set_flux_vars!(gwf)
    t = 0.0
    while t < dt
        gwf.variables.q_net .= 0.0
        dt_s = stable_timestep(gwf, conductivity_profile, cfl)
        dt_s = check_timestepsize(dt_s, t, dt)
        update_fluxes!(gwf, domain, conductivity_profile, dt_s)
        update_head!(gwf, soil, dt_s)
        update_ustorelayerdepth!(soil, gwf)
        t += dt_s
    end
    average_flux_vars!(gwf, dt)
    return nothing
end

get_water_depth(gwf::GroundwaterFlow) = gwf.parameters.top .- gwf.variables.head

get_exfiltwater(gwf::GroundwaterFlow) = gwf.variables.exfiltwater

function get_flux_to_river(subsurface_flow::GroundwaterFlow, inds::Vector{Int})
    (; river) = subsurface_flow.boundary_conditions
    flux = -river.variables.flux_av ./ tosecond(BASETIMESTEP) # [m³ s⁻¹]
    return flux
end

function sum_boundary_fluxes(gwf::GroundwaterFlow, domain::Domain; exclude = nothing)
    (; boundary_conditions) = gwf
    n = length(gwf.variables.storage)
    flux_in = zeros(n)
    flux_out = zeros(n)
    for bc in get_boundaries(boundary_conditions)
        isnothing(bc) && continue
        typeof(bc) == exclude && continue
        indices = get_boundary_index(bc, domain)
        for (i, index) in enumerate(indices)
            flux = bc.variables.flux_av[i]
            if flux > 0.0
                flux_in[index] += flux
            else
                flux_out[index] -= flux
            end
        end
    end
    return flux_in, flux_out
end
get_inflow(gwf::GroundwaterFlow) = gwf.variables.q_in_av
get_outflow(gwf::GroundwaterFlow) = gwf.variables.q_out_av
get_storage(gwf::GroundwaterFlow) = gwf.variables.storage
