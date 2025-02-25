
"Timestepping for solving kinematic wave and local inertial river and overland flow routing."
@with_kw struct TimeStepping
    stable_timesteps::Vector{Float64} = Float64[]
    dt_fixed::Float64 = 0.0
    adaptive::Bool = true
    cfl::Float64 = 0.70
end

"Check timestep size"
function check_timestepsize(timestepsize, currenttime, endtime)
    if currenttime + timestepsize > endtime
        timestepsize = endtime - currenttime
    end
    return timestepsize
end