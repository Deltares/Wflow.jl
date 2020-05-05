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
    "KsatVerFrac" => 1.0,
    "KsatHorFrac" => 1.0,
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
    if keys(nc.dim)[1] == "y" || "lat"
        trsp = true
    end

    subcatch_2d =
        trsp ? permutedims(nc["wflow_subcatch"][:]) : nc["wflow_subcatch"][:]
    # indices based on catchment
    inds = Wflow.active_indices(subcatch_2d, missing)
    n = length(inds)

    altitude = trsp ? Float64.(permutedims(nc["wflow_dem"][:])[inds]) :
        Float64.(nc["wflow_dem"][:][inds])
    river = trsp ? nomissing(permutedims(nc["wflow_river"][:])[inds], 0) :
        nomissing(nc["wflow_river"][:][inds], 0)
    riverwidth = trsp ? Float64.(permutedims(nc["wflow_riverwidth"][:])[inds]) :
        Float64.(nc["wflow_riverwidth"][:][inds])
    riverlength =
        trsp ? Float64.(permutedims(nc["wflow_riverlength"][:])[inds]) :
        Float64.(nc["wflow_riverlength"][:][inds])


    # read x, y coordinates and calculate cell length [m]
    y_nc = "y" in keys(nc.dim) ? nomissing(nc["y"][:]) : nomissing(nc["lat"][:])
    x_nc = "x" in keys(nc.dim) ? nomissing(nc["x"][:]) : nomissing(nc["lon"][:])
    y = repeat(y_nc, outer = (1, length(x_nc)))[inds]
    cellength = abs(mean(diff(x_nc)))

    # snow parameters (also set in ini file (snow=True or False)?)
    cfmax = readnetcdf(nc, "Cfmax", inds, dparams, transp = trsp)
    tt = readnetcdf(nc, "TT", inds, dparams, transp = trsp)
    tti = readnetcdf(nc, "TTI", inds, dparams, transp = trsp)
    ttm = readnetcdf(nc, "TTM", inds, dparams, transp = trsp)
    whc = readnetcdf(nc, "WHC", inds, dparams, transp = trsp)
    w_soil = readnetcdf(nc, "w_soil", inds, dparams, transp = trsp)
    cf_soil = readnetcdf(nc, "cf_soil", inds, dparams, transp = trsp)

    # soil parameters
    θₛ = readnetcdf(nc, "thetaS", inds, dparams, transp = trsp)
    θᵣ = readnetcdf(nc, "thetaR", inds, dparams, transp = trsp)
    kv₀ = readnetcdf(nc, "KsatVer", inds, dparams, transp = trsp)
    m = readnetcdf(nc, "M", inds, dparams, transp = trsp)
    hb = readnetcdf(nc, "AirEntryPressure", inds, dparams, transp = trsp)
    soilthickness =
        readnetcdf(nc, "SoilThickness", inds, dparams, transp = trsp)
    infiltcappath =
        readnetcdf(nc, "InfiltCapPath", inds, dparams, transp = trsp)
    infiltcapsoil =
        readnetcdf(nc, "InfiltCapSoil", inds, dparams, transp = trsp)
    maxleakage = readnetcdf(nc, "MaxLeakage", inds, dparams, transp = trsp)
    #TODO: store c, kvfrac in staticmaps.nc start at index 1
    c = fill(dparams["c"], (maxlayers, n))
    kvfrac = fill(dparams["KsatVerFrac"], (maxlayers, n))
    for i in [0:1:maxlayers-1;]
        if string("c_", i) in keys(nc)
            c[i+1, :] =
                trsp ? Float64.(permutedims(nc[string("c_", i)][:])[inds]) :
                Float64.(nc[string("c_", i)][:][inds])
        else
            @warn(string(
                "c_",
                i,
                " not found, set to default value ",
                dparams["c"],
            ))
        end
        if string("KsatVerFrac_", i) in keys(nc)
            kvfrac[i+1, :] = trsp ?
                Float64.(permutedims(nc[string("KsatVerFrac_", i)][:])[inds]) :
                Float64.(nc[string("KsatVerFrac_", i)][:][inds])
        else
            @warn(string(
                "KsatVerFrac_",
                i,
                " not found, set to default value ",
                dparams["KsatVerFrac"],
            ))
        end
    end

    # fraction open water and compacted area (land cover)
    waterfrac = readnetcdf(nc, "WaterFrac", inds, dparams, transp = trsp)
    pathfrac = readnetcdf(nc, "PathFrac", inds, dparams, transp = trsp)

    # vegetation parameters
    rootingdepth = readnetcdf(nc, "RootingDepth", inds, dparams, transp = trsp)
    rootdistpar = readnetcdf(nc, "rootdistpar", inds, dparams, transp = trsp)
    capscale = readnetcdf(nc, "CapScale", inds, dparams, transp = trsp)
    et_reftopot = readnetcdf(nc, "et_reftopot", inds, dparams, transp = trsp)
    # cmax, e_r, canopygapfraction only required when lai climatoly not provided
    cmax = readnetcdf(nc, "Cmax", inds, dparams, transp = trsp)
    e_r = readnetcdf(nc, "EoverR", inds, dparams, transp = trsp)
    canopygapfraction =
        readnetcdf(nc, "CanopyGapFraction", inds, dparams, transp = trsp)

    # if lai climatology provided use sl, swood and kext to calculate cmax
    if isnothing(leafarea_path) == false
        sl = readnetcdf(nc, "Sl", inds, dparams, transp = trsp)
        swood = readnetcdf(nc, "Swood", inds, dparams, transp = trsp)
        kext = readnetcdf(nc, "Kext", inds, dparams, transp = trsp)
        # set in inifile? Also type (monthly, daily, hourly) as part of netcdf variable attribute?
        # in original inifile: LAI=staticmaps/clim/LAI,monthlyclim,1.0,1
        lai_clim = NCDataset(leafarea_path) #TODO:include LAI climatology in update() vertical SBM model
    end

    xl = fill(mv, n)
    yl = fill(mv, n)

    sbm = Vector{SBM}(undef, n)
    for i = 1:n
        act_thickl = set_layerthickness(soilthickness[i], sumlayers)
        nlayers = length(act_thickl)
        s_layers = pushfirst(cumsum(SVector{nlayers,Float64}(act_thickl)), 0.0)

        xl[i] = sizeinmetres ? cellength : lattometres(y[i])[1] * cellength
        yl[i] = sizeinmetres ? cellength : lattometres(y[i])[2] * cellength
        riverfrac = Bool(river[i]) ?
            min((riverlength[i] * riverwidth[i]) / (xl[i] * yl[i]), 1.0) : 0.0

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
            kv₀ = kv₀[i],
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

    # lateral part sbm
    khfrac = readnetcdf(nc, "KsatHorFrac", inds, dparams, transp = trsp)
    βₗ = trsp ? Float64.(permutedims(nc["Slope"][:])[inds]) :
        Float64.(nc["Slope"][:][inds])
    ldd =
        trsp ? permutedims(nc["wflow_ldd"][:])[inds] : nc["wflow_ldd"][:][inds]
    kh₀ = khfrac .* kv₀
    dl = fill(mv, n)
    dw = fill(mv, n)

    for i = 1:n
        dl[i] = detdrainlength(ldd[i], xl[i], yl[i])
        dw[i] = detdrainwidth(ldd[i], xl[i], yl[i])
    end

    ssf = LateralSSF{Float64,n}(
        kh₀ = kh₀,
        f = getfield.(sbm, :f),
        zi = getfield.(sbm, :zi),
        soilthickness = soilthickness,
        θₑ = θₛ - θᵣ,
        Δt = Δt,
        βₗ = βₗ,
        dl = dl * 1000.0,
        dw = dw * 1000.0,
    )

    dag = flowgraph(ldd, inds, Wflow.pcrdir)
    starttime = DateTime(2000, 1, 1)

    # create a Model
    model = Model(dag, ssf, sbm, Clock(starttime, 1, Δt), nothing, nothing)

end
