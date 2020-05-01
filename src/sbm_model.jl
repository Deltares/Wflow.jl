# timestep that the parameter units are defined in
basetimestep = Second(Day(1))
Δt = Second(Day(1))

# default parameter values (dict)
dparams = Dict(
    "Cfmax" => 3.75653 * (Δt / basetimestep),
    "TT" => 0.0,
    "TTM" => 0.0,
    "TTI" => 1.0,
    "WHC" => 0.1,
    "cf_soil" => 0.038,
    "w_soil" => 0.1125 * (Δt / basetimestep),
    "SoilThickness" => 2000.0,
    "InfiltCapSoil" => 100.0,
    "InfiltCapPath" => 10.0,
    "PathFrac" => 0.01,
    "WaterFrac" => 0.0,
    "thetaS" => 0.6,
    "thetaR" => 0.01,
    "AirEntryPressure" => 10.0,
    "KsatVer" => 3000.0 * (Δt / basetimestep),
    "MaxLeakage" => 0.0,
    "c" => 10.0,
    "M" => 300.0,
    "CapScale" => 100.0,
    "rootdistpar" => -500.0,
    "RootingDepth" => 750.0,
    "LAI" => 1.0,
    "Cmax" => 1.0,
    "CanopyGapFraction" => 0.1,
    "EoverR" => 0.1,
    "et_reftopot" => 1.0,
    "kvfrac" => 1.0,
)

"""
    initialize_sbm_model(staticmaps_path, leafarea_path)

Initial part of the SBM model concept. Reads model parameters from disk, `staticmaps_path` is the file path
of the NetCDF file with model parameters, `leafarea_path` is an optional file path for a NetCDF file with leaf
area index (LAI) values (climatology).
"""
function initialize_sbm_model(staticmaps_path, leafarea_path)

    sizeinmetres = false
    thicknesslayers = SVector(100.0, 300.0, 800.0)
    maxlayers = length(thicknesslayers) + 1 # max number of soil layers
    sumlayers = SVector(pushfirst(cumsum(thicknesslayers), 0.0))

    nc = NCDataset(staticmaps_path)

    subcatch_2d = nc["wflow_subcatch"][:]
    # indices based on catchment
    inds = Wflow.active_indices(subcatch_2d, missing)
    n = length(inds)

    altitude = Float64.(nc["wflow_dem"][:][inds])
    river = nomissing(nc["wflow_river"][:], 0)[inds]
    riverwidth = Float64.(nc["wflow_riverwidth"][:][inds])
    ldd = Float64.(nc["wflow_ldd"][:][inds])
    if "wflow_riverlength" in keys(nc)
        riverlength = Float64.(nc["wflow_riverlength"][:][inds])
    else
        @warn("wflow_riverlength not found, riverlength based on ldd...")
        # TODO calculate river based on ldd
    end

    # read x, y coordinates and calculate cell length [m]
    y_nc = "y" in keys(nc.dim) ? nomissing(nc["y"][:]) : nomissing(nc["lat"][:])
    x_nc = "x" in keys(nc.dim) ? nomissing(nc["x"][:]) : nomissing(nc["lon"][:])
    y = repeat(y_nc', outer = (length(x_nc), 1))[inds]
    cellength = abs(mean(diff(x_nc)))

    # snow parameters (also set in ini file (snow=True or False)?)
    cfmax = readnetcdf(nc, "Cfmax", inds, dparams)
    tt = readnetcdf(nc, "TT", inds, dparams)
    tti = readnetcdf(nc, "TTI", inds, dparams)
    ttm = readnetcdf(nc, "TTM", inds, dparams)
    whc = readnetcdf(nc, "WHC", inds, dparams)
    w_soil = readnetcdf(nc, "w_soil", inds, dparams)
    cf_soil = readnetcdf(nc, "cf_soil", inds, dparams)

    # soil parameters
    θₛ = readnetcdf(nc, "thetaS", inds, dparams)
    θᵣ = readnetcdf(nc, "thetaR", inds, dparams)
    kv = readnetcdf(nc, "KsatVer", inds, dparams)
    m = readnetcdf(nc, "M", inds, dparams)
    hb = readnetcdf(nc, "AirEntryPressure", inds, dparams)
    soilthickness = readnetcdf(nc, "SoilThickness", inds, dparams)
    infiltcappath = readnetcdf(nc, "InfiltCapPath", inds, dparams)
    infiltcapsoil = readnetcdf(nc, "InfiltCapSoil", inds, dparams)
    maxleakage = readnetcdf(nc, "MaxLeakage", inds, dparams)
    #TODO: store c, kvfrac in staticmaps.nc start at index 1
    c = fill(dparams["c"], (maxlayers, n))
    kvfrac = fill(dparams["kvfrac"], (maxlayers, n))
    for i in [0:1:maxlayers-1;]
        if string("c_", i) in keys(nc)
            c[i+1, :] = Float64.(nc[string("c_", i)][:][inds])
        else
            @warn(string("c_", i, " not found, set to default value ", dparams["c"]))
        end
        if string("kvfrac_", i) in keys(nc)
            kvfrac[i+1, :] = Float64.(nc[string("kvfrac_", i)][:][inds])
        else
            @warn(string(
                "kvfrac_",
                i,
                " not found, set to default value ",
                dparams["kvfrac"],
            ))
        end
    end

    # fraction open water and compacted area (land cover)
    waterfrac = readnetcdf(nc, "WaterFrac", inds, dparams)
    pathfrac = readnetcdf(nc, "PathFrac", inds, dparams)

    # vegetation parameters
    rootingdepth = readnetcdf(nc, "RootingDepth", inds, dparams)
    rootdistpar = readnetcdf(nc, "rootdistpar", inds, dparams)
    capscale = readnetcdf(nc, "CapScale", inds, dparams)
    et_reftopot = readnetcdf(nc, "et_reftopot", inds, dparams)
    # cmax, e_r, canopygapfraction only required when lai climatoly not provided
    cmax = readnetcdf(nc, "Cmax", inds, dparams)
    e_r = readnetcdf(nc, "EoverR", inds, dparams)
    canopygapfraction = readnetcdf(nc, "CanopyGapFraction", inds, dparams)

    # if lai climatology provided use sl, swood and kext to calculate cmax
    if isnothing(leafarea_path) == false
        sl = readnetcdf(nc, "Sl", inds, dparams)
        swood = readnetcdf(nc, "Swood", inds, dparams)
        kext = readnetcdf(nc, "Kext", inds, dparams)
        # set in inifile? Also type (monthly, daily, hourly) as part of netcdf variable attribute?
        # in original inifile: LAI=staticmaps/clim/LAI,monthlyclim,1.0,1
        lai_clim = NCDataset(leafarea_path) #TODO:include LAI climatology in update() vertical SBM model
    end

    sbm = Vector{SBM}(undef, n)
    for i = 1:n
        act_thickl = set_layerthickness(soilthickness[i], sumlayers)
        nlayers = length(act_thickl)
        s_layers = pushfirst(cumsum(SVector{nlayers,Float64}(act_thickl)), 0.0)

        xl = sizeinmetres ? cellength : lattometres(y[i])[1] * cellength
        yl = sizeinmetres ? cellength : lattometres(y[i])[2] * cellength
        riverfrac =
            Bool(river[i]) ? min((riverlength[i] * riverwidth[i]) / (xl * yl), 1.0) : 0.0

        sbm[i] = SBM{Float64,nlayers,nlayers + 1}(
            maxlayers = maxlayers,
            nlayers = nlayers,
            riverfrac = riverfrac,
            cfmax = cfmax[i],
            tt = tt[i],
            tti = tti[i],
            ttm = ttm[i],
            whc = whc[i],
            w_soil = w_soil[i],
            cf_soil = cf_soil[i],
            θₛ = θₛ[i],
            θᵣ = θᵣ[i],
            kv = kv[i],
            kvfrac = kvfrac[1:nlayers, i],
            m = m[i],
            hb = hb[i],
            soilthickness = soilthickness[i],
            act_thickl = act_thickl,
            sumlayers = s_layers,
            infiltcappath = infiltcappath[i],
            infiltcapsoil = infiltcapsoil[i],
            maxleakage = maxleakage[i],
            waterfrac = max(waterfrac[i] - riverfrac, 0.0),
            pathfrac = pathfrac[i],
            altitude = altitude[i],
            rootingdepth = rootingdepth[i],
            rootdistpar = rootdistpar[i],
            capscale = capscale[i],
            et_reftopot = et_reftopot[i],
            sl = sl[i],
            swood = swood[i],
            kext = kext[i],
            c = c[1:nlayers, i],
            lai = 1.0,
            cmax = cmax[i],
            canopygapfraction = canopygapfraction[i],
            e_r = e_r[i],
        )

    end

    return sbm

end
