# timestep that the parameter units are defined in
const basetimestep = Second(Day(1))
const Δt = Second(Day(1))

# default parameter values (dict)
const dparams = Dict(
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
    thicknesslayers = SVector(100.0, 300.0, 800.0, mv)
    maxlayers = length(thicknesslayers) # max number of soil layers
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
    # TODO: store c, kvfrac in staticmaps.nc start at index 1
    c = fill(dparams["c"], (maxlayers, n))
    kvfrac = fill(dparams["KsatVerFrac"], (maxlayers, n))
    for i in [0:1:maxlayers - 1;]
        if string("c_", i) in keys(nc)
            c[i + 1, :] =
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
            kvfrac[i + 1, :] = trsp ?
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
        lai_clim = NCDataset(leafarea_path) # TODO:include LAI climatology in update() vertical SBM model
    end

    # these are filled in the loop below
    # TODO see if we can replace this approach
    nlayers = zeros(Int, n)
    act_thickl = zeros(Float64, maxlayers, n)
    s_layers = zeros(Float64, maxlayers+1, n)
    xl = fill(mv, n)
    yl = fill(mv, n)
    riverfrac = fill(mv, n)

    for i = 1:n
        act_thickl_, nlayers_ = set_layerthickness(soilthickness[i], sumlayers, thicknesslayers)
        s_layers_ = pushfirst(cumsum(act_thickl_), 0.0)

        xl[i] = sizeinmetres ? cellength : lattometres(y[i])[1] * cellength
        yl[i] = sizeinmetres ? cellength : lattometres(y[i])[2] * cellength
        riverfrac[i] = Bool(river[i]) ?
            min((riverlength[i] * riverwidth[i]) / (xl[i] * yl[i]), 1.0) : 0.0

        nlayers[i] = nlayers_
        act_thickl[:, i] = act_thickl_
        s_layers[:, i] = s_layers_
    end

    # needed for derived parameters below
    act_thickl = svectorscopy(act_thickl, Val{maxlayers}())
    θₑ = θₛ .- θᵣ
    soilwatercapacity = soilthickness .* θₑ
    satwaterdepth = 0.85 .* soilwatercapacity

    # copied to array of sarray below
    vwc = fill(mv, maxlayers, n)
    vwc_perc = fill(mv, maxlayers, n)

    sbm = Table(
        maxlayers = Fill(maxlayers, n), # Fill{Float64}
        nlayers = nlayers,  # Vector{Float64}
        riverfrac = riverfrac, # Vector{Float64}
        cfmax = cfmax,  # Vector{Float64}
        tt = tt,  # Vector{Float64}
        tti = tti,  # Vector{Float64}
        ttm = ttm,  # Vector{Float64}
        whc = whc,  # Vector{Float64}
        w_soil = w_soil,  # Vector{Float64}
        cf_soil = cf_soil,  # Vector{Float64}
        θₛ = θₛ,  # Vector{Float64}
        θᵣ = θᵣ,  # Vector{Float64}
        kv₀ = kv₀,  # Vector{Float64}
        kvfrac = svectorscopy(kvfrac, Val{maxlayers}()),  # Vector{SVector{maxlayers,Float64}}
        m = m,  # Vector{Float64}
        hb = hb,  # Vector{Float64}
        soilthickness = soilthickness,  # Vector{Float64}
        act_thickl = act_thickl,  # Vector{SVector{maxlayers,Float64}}
        sumlayers = svectorscopy(s_layers, Val{maxlayers+1}()),  # Vector{SVector{maxlayers,Float64}}
        infiltcappath = infiltcappath,  # Vector{Float64}
        infiltcapsoil = infiltcapsoil,  # Vector{Float64}
        maxleakage = maxleakage,  # Vector{Float64}
        waterfrac = max.(waterfrac .- riverfrac, 0.0),  # Vector{Float64}
        pathfrac = pathfrac,  # Vector{Float64}
        altitude = altitude,  # Vector{Float64}
        rootingdepth = rootingdepth,  # Vector{Float64}
        rootdistpar = rootdistpar,  # Vector{Float64}
        capscale = capscale,  # Vector{Float64}
        et_reftopot = et_reftopot,  # Vector{Float64}
        sl = sl,  # Vector{Float64}
        swood = swood,  # Vector{Float64}
        kext = kext,  # Vector{Float64}
        c = svectorscopy(c, Val{maxlayers}()),  # Vector{SVector{maxlayers,Float64}}
        lai = Fill(1.0, n), # Fill{Float64}
        cmax = cmax,  # Vector{Float64}
        canopygapfraction = canopygapfraction,  # Vector{Float64}
        e_r = e_r,  # Vector{Float64}

        # filled in by SBM struct defaults
        f = θₑ ./ m,
        ustorelayerdepth = act_thickl .* 0.0,
        satwaterdepth = satwaterdepth,
        zi = max.(0.0, soilthickness .- satwaterdepth ./ θₑ),
        soilwatercapacity = soilwatercapacity,
        snow = fill(0.0, n),
        snowwater = fill(0.0, n),
        tsoil = Fill(10.0, n),
        canopystorage = fill(0.0, n),

        # preallocated but set to mv
        # TODO check if we can use FillArrays for some
        precipitation = fill(mv, n),       # Precipitation [mm]
        temperature = fill(mv, n),         # Temperature [ᵒC]
        potevap = fill(mv, n),             # Potential evapotranspiration [mm]
        pottrans_soil = fill(mv, n),       # Potential transpiration, open water, river and soil evaporation (after subtracting interception from potevap)
        transpiration = fill(mv, n),       # Transpiration [mm]
        ae_ustore = fill(mv, n),           # Actual evaporation from unsaturated store [mm]
        ae_sat = fill(mv, n),              # Actual evaporation from saturated store [mm]
        interception = fill(mv, n),        # Interception [mm]
        soilevap = fill(mv, n),            # Soil evaporation [mm]
        actevapsat = fill(mv, n),          # Actual evaporation from saturated store (transpiration and soil evaporation) [mm]
        actevap = fill(mv, n),             # Total actual evaporation (transpiration + soil evapation + open water evaporation) [mm]
        ae_openw_l = fill(mv, n),          # Actual evaporation from open water (land) [mm]
        ae_openw_r = fill(mv, n),          # Actual evaporation from river [mm]
        avail_forinfilt = fill(mv, n),     # Water available for infiltration [mm]
        actinfilt = fill(mv, n),           # Actual infiltration into the unsaturated zone [mm]
        actinfiltsoil = fill(mv, n),       # Actual infiltration non-compacted fraction [mm]
        actinfiltpath = fill(mv, n),       # Actual infiltration compacted fraction [mm]
        infiltexcess = fill(mv, n),        # Infiltration excess water [mm]
        excesswater = fill(mv, n),         # Water that cannot infiltrate due to saturated soil (saturation excess) [mm]
        exfiltsatwater = fill(mv, n),      # Water exfiltrating during saturation excess conditions [mm]
        exfiltustore = fill(mv, n),        # Water exfiltrating from unsaturated store because of change in water table [mm]
        excesswatersoil = fill(mv, n),     # Excess water for non-compacted fraction [mm]
        excesswaterpath = fill(mv, n),     # Excess water for compacted fraction [mm]
        runoff = fill(mv, n),              # Total surface runoff from infiltration and saturation excess [mm]
        vwc = svectorscopy(vwc, Val{maxlayers}()),           # Volumetric water content [mm mm⁻¹] per soil layer (including θᵣ and saturated zone)
        vwc_perc = svectorscopy(vwc_perc, Val{maxlayers}()),      # Volumetric water content [%] per soil layer (including θᵣ and saturated zone)
        rootstore = fill(mv, n),           # Root water storage [mm] in unsaturated and saturated zone (excluding θᵣ)
        vwc_root = fill(mv, n),            # Volumetric water content [mm mm⁻¹] in root zone (including θᵣ and saturated zone)
        vwc_percroot = fill(mv, n),        # Volumetric water content [%] in root zone (including θᵣ and saturated zone)
        ustoredepth = fill(mv, n),         # Amount of available water in the unsaturated zone [mm]
        transfer = fill(mv, n),            # Downward flux from unsaturated to saturated zone [mm]
        capflux = fill(mv, n),             # Capilary rise [mm]
        recharge = fill(mv, n),            # Net recharge to saturated store [mm]

    )

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

    ssf = LateralSSF{Float64}(
        kh₀ = kh₀,
        f = sbm.f,
        zi = sbm.zi,
        soilthickness = soilthickness,
        θₑ = θₛ .- θᵣ,
        Δt = 1.0,
        βₗ = βₗ,
        dl = dl .* 1000.0,
        dw = dw .* 1000.0,
    )

    dag = flowgraph(ldd, inds, Wflow.pcrdir)
    starttime = DateTime(2000, 1, 1)

    model = Model(dag, ssf, sbm, Clock(starttime, 1, Δt), nothing, nothing)
    return model
end
