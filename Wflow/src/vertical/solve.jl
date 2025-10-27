"""
A WflowIntegrator object holds all parameters and variables
needed for numerical integration of an ODE.
"""
@kwdef mutable struct WflowIntegrator{V <: ComponentVector}
    # The fraction of the global time step that has been simulated
    progress::Float64 = 0.0
    # The sub time step
    dt_sub::Float64 = 0.0
    # The state vector
    const u::V
    # The previous value of the state vector
    const uprev::V
    # The time derivative of the state vector,
    # generally instantaneous fluxes
    const du::V
end

function WflowIntegrator(u0::V) where {V <: ComponentVector}
    return WflowIntegrator(; u = u0, uprev = copy(u0), du = zero(u0))
end

const StateType{A} = ComponentVector{Float64, Vector{Float64}, Tuple{Axis{A}}}

# Fallback method; no sub time steps
function set_sub_time_step!(model_vertical::AbstractVerticalModel, i::Integer)::Nothing
    integrator = model_vertical.integrators[i]
    integrator.dt_sub = 1.0
    integrator.progress = 1.0
    return nothing
end

# Fallback method; do nothing
dt_sub_callback!(::AbstractVerticalModel, ::LandHydrologySBM, ::Int) = nothing

# Fallback method; do nothing (TODO: Make more generic)
update!(::NoSnowModel, ::LandHydrologySBM) = nothing

function update!(
    model_vertical::M,
    land::LandHydrologySBM,
)::Nothing where {M <: AbstractVerticalModel}
    (; integrators, p) = model_vertical
    n = length(integrators)

    # Loop over the cells
    threaded_foreach(1:n; basesize = 1000) do i
        integrator = integrators[i]
        # Perform pre processing before sub time stepping;
        # things that only have to be computed once for the whole
        # global time step
        update_preamble!(model_vertical, land, i)
        (; du, u, uprev) = integrator

        # The fraction of the global time step that has passed 
        integrator.progress = 0.0

        # Perform sub time steps until the global time step has been reached
        while integrator.progress < 1.0
            # Set the instantaneous rates of changes of the states
            set_instantaneous_rates!(du, u, p, i, integrator.progress)

            # Get the sub time step 
            set_sub_time_step!(model_vertical, i)
            @assert integrator.progress <= 1.0 "Overstepped the global time step for $(nameof(M)) at cell $i."

            # Perform an Euler forward step
            @. uprev = u
            @. u += integrator.dt_sub * du
            integrator.progress += integrator.dt_sub

            # Perform post processing after sub time step
            dt_sub_callback!(model_vertical, land, i)
        end

        # Perform post processing after global time step
        update_postamble!(model_vertical, land, i)
    end
    return nothing
end