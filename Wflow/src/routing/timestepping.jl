
"Timestepping for river, overland and subsurface flow routing."
@with_kw struct TimeStepping
    stable_timesteps::Vector{Float64} = Float64[]
    dt_fixed::Float64 = 0.0
    adaptive::Bool = true
    alpha_coefficient::Float64 = 1.0
end

"Check timestep size"
function check_timestepsize(timestepsize, currenttime, endtime)
    if currenttime + timestepsize > endtime
        timestepsize = endtime - currenttime
    end
    return timestepsize
end

"Initialize timestepping for kinematic wave river, overland and lateral subsurface flow models"
function init_kinematic_wave_timestepping(config::Config, n::Int; domain::String)
    adaptive = config.model.kinematic_wave__adaptive_time_step_flag
    @info "Kinematic wave approach is used for $domain flow, adaptive timestepping = $adaptive."
    if adaptive
        stable_timesteps = zeros(n)
        if domain == "subsurface"
            alpha_coefficient = config.model.subsurface_kinematic_wave__alpha_coefficient
            @info "Numerical stability coefficient for lateral subsurface flow `alpha`: `$alpha_coefficient`."
            timestepping = TimeStepping(; stable_timesteps, adaptive, alpha_coefficient)
        else
            timestepping = TimeStepping(; stable_timesteps, adaptive)
        end
    else
        dt_fixed = getfield(config.model, Symbol("$(domain)_kinematic_wave__time_step"))
        @info "Using a fixed internal timestep (seconds) $dt_fixed for kinematic wave $domain flow."
        timestepping = TimeStepping(; dt_fixed, adaptive)
    end
    return timestepping
end
