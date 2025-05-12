
"Timestepping for solving kinematic wave and local inertial river and overland flow routing."
@with_kw struct TimeStepping{T <: DenseArray{Float}}
    stable_timesteps::T = Float[]
    dt_fixed::Float = 0.0
    adaptive::Bool = true
    cfl::Float = 0.70
end

function Adapt.adapt_structure(to, from::TimeStepping)
    return TimeStepping(
        adapt(to, from.stable_timesteps),
        adapt(to, from.dt_fixed),
        adapt(to, from.adaptive),
        adapt(to, from.cfl),
    )
end

"Check timestep size"
function check_timestepsize(timestepsize, currenttime, endtime)
    if currenttime + timestepsize > endtime
        timestepsize = endtime - currenttime
    end
    return Float(timestepsize)
end