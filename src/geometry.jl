
@get_units @grid_loc @with_kw struct LandGeometry{T}
    # cell area [m^2]
    area::Vector{T} | "m^2"
    # drain width [m]
    width::Vector{T} | "m"
    # drain slope
    slope::Vector{T} | "-"
end

function LandGeometry(nc, config, inds)
    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc; outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cellength, sizeinmetres)
    area = xl .* yl
    ldd = ncread(nc, config, "ldd"; optional = false, sel = inds, allow_missing = true)
    drain_width = map(detdrainwidth, ldd, xl, yl)
    landslope = ncread(
        nc,
        config,
        "vertical.geometry.slope";
        optional = false,
        sel = inds,
        type = Float,
    )
    clamp!(landslope, 0.00001, Inf)

    land_geometry =
        LandGeometry{Float}(; area = area, width = drain_width, slope = landslope)
    return land_geometry
end

@get_units @grid_loc @with_kw struct RiverGeometry{T}
    # drain width [m]
    width::Vector{T} | "m"
    # drain length
    length::Vector{T} | "m"
    # slope
    slope::Vector{T} | "-"
end

function RiverGeometry(nc, config, inds)
    riverwidth = ncread(
        nc,
        config,
        "lateral.river.geometry.width";
        optional = false,
        sel = inds,
        type = Float,
    )
    riverlength = ncread(
        nc,
        config,
        "lateral.river.geometry.length";
        optional = false,
        sel = inds,
        type = Float,
    )
    riverslope = ncread(
        nc,
        config,
        "lateral.river.geometry.slope";
        optional = false,
        sel = inds,
        type = Float,
    )
    minimum(riverlength) > 0 || error("river length must be positive on river cells")
    minimum(riverwidth) > 0 || error("river width must be positive on river cells")
    clamp!(riverslope, 0.00001, Inf)

    river_geometry =
        RiverGeometry{Float}(; width = riverwidth, length = riverlength, slope = riverslope)
    return river_geometry
end