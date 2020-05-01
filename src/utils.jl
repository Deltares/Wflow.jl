
function lattometres(lat::Float64)
    m1 = 111132.92     #latitude calculation term 1
    m2 = -559.82       #latitude calculation term 2
    m3 = 1.175         #latitude calculation term 3
    m4 = -0.0023       #latitude calculation term 4
    p1 = 111412.84     #longitude calculation term 1
    p2 = -93.5         #longitude calculation term 2
    p3 = 0.118         #longitude calculation term 3

    #Calculate the length of a degree of latitude and longitude in meters
    latlen =
        m1 +
        (m2 * cosd(2.0 * lat)) +
        (m3 * cosd(4.0 * lat)) +
        (m4 * cosd(6.0 * lat))
    longlen = (p1 * cosd(lat)) + (p2 * cosd(3.0 * lat)) + (p3 * cosd(5.0 * lat))

    return longlen, latlen
end

"""
    readnetcdf(nc, var, inds, dpars)

Read parameter `var` from NetCDF file `nc` for indices `inds`. If `var` is not
available, a default value based on dict 'dpars' is returned.
"""
function readnetcdf(nc, var, inds, dpars)
    if haskey(nc, var)
        @info(string("read parameter ", var))
        Float64.(nc[var][:][inds])
    else
        @warn(string(var, " not found, set to default value ", dpars[var]))
        fill(dpars[var], length(inds))
    end
end

"""
    set_layerthickness(d::Float64, sl::SVector)

Calculate actual soil thickness of layers based on a reference depth (e.g. soil depth or water table depth) `d`,
and a SVector `sl` with cumulative soil depth starting at soil surface (0).
"""
function set_layerthickness(d::Float64, sl::SVector)
    act_d = sl[sl.<d]
    if d - act_d[end] > 0
        push!(act_d, d)
    end
    diff(act_d)
end
