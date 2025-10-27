function update!(
    model_vertical::M,
    model::Model,
)::Nothing where {M <: AbstractVerticalModel}
    (; integrators, p) = model_vertical
    n = length(integrators)

    # Loop over the cells
    threaded_foreach(1:n; basesize = 1000) do i
        integrator = integrators[i]
        # Perform pre processing before sub time stepping;
        # things that only have to be computed once for the whole
        # global time step
        update_preamble!(model_vertical, model, i)
        (; du, u, uprev) = integrator

        # The fraction of the global time step that has passed 
        integrator.progress = 0.0

        # Perform sub time steps until the global time step has been reached
        while integrator.progress < 1.0
            # Set the instantaneous rates of changes of the states
            set_instantaneous_rates!(du, u, p, i, integrator.progress)

            # Get the sub time step 
            set_sub_time_step!(model_vertical, i)
            @assert progress <= 1.0 "Overstepped the global time step for $(nameof(M)) at cell $i."

            # Perform an Euler forward step
            @. uprev = u
            @. u += integrator.dt_sub * du
            integrator.progress += integrator.dt_sub

            # Perform post processing after sub time step
            dt_sub_callback!(model_vertical, model, i)
        end

        # Perform post processing after global time step
        update_postamble!(model_vertical, model, i)
    end
    return nothing
end