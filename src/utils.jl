
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
function readnetcdf(nc, var, inds, dpars; transp = false)
    if haskey(nc, var)
        @info(string("read parameter ", var))
        ncvar = transp ? Float64.(permutedims(nc[var][:])[inds]) :
            Float64.(nc[var][:][inds])
    else
        @warn(string(var, " not found, set to default value ", dpars[var]))
        ncvar = fill(dpars[var], length(inds))
    end
end

"""
    set_layerthickness(d::Float64, sl::SVector)

Calculate actual soil thickness of layers based on a reference depth (e.g. soil depth or water table depth) `d`,
a SVector `sl` with cumulative soil depth starting at soil surface (0), and a SVector `tl` with actual thickness
per soil layer.
"""
function set_layerthickness(d::Float64, sl::SVector, tl::SVector)

    act_d = tl .* mv
    for i = 1:length(act_d)
        if d > sl[i+1]
            act_d = setindex(act_d, tl[i], i)
        elseif d-sl[i] > 0.0
            act_d = setindex(act_d, d-sl[i], i)
        end
    end

    nlayers = length(act_d) - sum(isnan.(act_d))
    return act_d, nlayers
end

"""
    detdrainlength(ldd, xl, yl)

Determines the drainaige length for a non square grid. Input `ldd` (drainage network), `xl` (length of cells in x direction),
`yl` (length of cells in y direction). Output is drainage length.
"""
function detdrainlength(ldd, xl, yl)
    # take into account non-square cells
    # if ldd is 8 or 2 use ylength
    # if ldd is 4 or 6 use xlength
    if ldd == 2 || ldd ==  8
        yl
    elseif ldd == 4 || ldd == 6
        xl
    else
        sqrt(xl^2 + yl^2)
    end
end

"""
    detdrainwidth(ldd, xl, yl)

Determines the drainaige width for a non square grid. Input `ldd` (drainage network), `xl` (length of cells in x direction),
`yl` (length of cells in y direction). Output is drainage width.
"""

function detdrainwidth(ldd, xl, yl)
    # take into account non-square cells
    # if ldd is 8 or 2 use xlength
    # if ldd is 4 or 6 use ylength
    slantwidth = (xl + yl) * 0.5
    if ldd == 2 || ldd == 8
        xl
    elseif ldd == 4 || ldd == 6
        yl
    else
        slantwidth
    end
end
