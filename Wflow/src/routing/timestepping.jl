
"Timestepping for solving kinematic wave and local inertial river and overland flow routing."
@with_kw struct TimeStepping{T <: AbstractArray{<:AbstractFloat}, S <: AbstractFloat}
    stable_timesteps::T = Float64[]
    dt_fixed::S = 0.0
    adaptive::Bool = true
    cfl::S = 0.70
end
@adapt_structure TimeStepping

"Check timestep size"
function check_timestepsize(timestepsize, currenttime, endtime)
    if currenttime + timestepsize > endtime
        timestepsize = endtime - currenttime
    end
    return timestepsize
end