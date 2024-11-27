@with_kw struct TimeStepping{T}
    stable_timesteps::Vector{T} = Float[]
    dt_fixed::T = 0.0
    adaptive::Bool = true
    cfl::T = 0.70
end

function check_timestepsize(timestepsize, currenttime, endtime)
    if currenttime + timestepsize > endtime
        timestepsize = endtime - currenttime
    end
    return timestepsize
end