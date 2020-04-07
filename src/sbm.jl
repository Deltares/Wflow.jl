const mv = NaN

# from ini file
# [run]
# timestepsecs = 86400
timestepsecs = 86400.0

Base.@kwdef struct SBM{N}
    cfmax::Float64
    tt::Float64
    ttm::Float64
    tti::Float64
    whc::Float64
    cf_soil::Float64
    w_soil::Float64
    soilthickness::Float64
    infiltcapsoil::Float64
    infiltcappath::Float64
    pathfrac::Float64
    waterfrac::Float64
    θₛ::Float64
    θᵣ::Float64
    hb ::Float64
    kv ::Float64
    maxleakage::Float64
    c::SVector{N,Float64}
    m::Float64
    f::Float64 = (θₛ - θᵣ) / m
    capscale::Float64
    rootdistpar::Float64
    rootingdepth::Float64
    lai::Float64
    sl::Float64
    kext::Float64
    swood::Float64
    et_reftopot::Float64
    altitude::Float64
    precipitation::Float64 = mv
    temperature::Float64 = mv
    potevap::Float64 = mv
    pottrans_soil::Float64 = mv
    transpiration::Float64 = mv
    ae_ustore::Float64 = mv
    ae_sat::Float64 = mv
    interception::Float64 = mv
    ae::Float64 = mv
    ae_openw_l::Float64 = mv
    ae_openw_r::Float64 = mv
    avail_forinfilt::Float64 = mv
    zi::Float64 = mv
    ustorelayerdepth::SVector{N,Float64} = fill(mv, SVector{N,Float64}) #TODO:define nLayers per grid cell
    ustoredepth::Float64 = mv
    transfer::Float64 = mv
    capflux::Float64 = mv
    recharge::Float64 = mv
    soilwatercapacity::Float64 = soilthickness * (θₛ - θᵣ)
    satwaterdepth::Float64 = 0.85 * soilwatercapacity
    snow::Float64 = 0.0
    snowwater::Float64 = 0.0
    tsoil::Float64 = 10.0
    canopystorage::Float64 = 0.0
end

function readnetcdf(nc, var, inds, dpars)
    if haskey(nc, var)
        @info(string("read parameter ", var))
        Float64.(nc[var][:][inds])
    else
        @warn(string(var, " not found, set to default value ", dpars[var]))
        fill(dpars[var], length(inds))
    end
end

function statenames()

    # depends on ini file settings (optional: glaciers, snow, irrigation)
    states = [
        :satwaterdepth,
        :snow,
        :tsoil,
        :ustorelayerdepth,
        :snowwater,
        :canopystorage,
    ]
    #TODO: (warm) states read from netcdf file or cold state (reinit=1, setting in ini file)

end

"Initial part of the model. Reads model parameters from disk"
function initialize(staticmaps_path, leafarea_path)

    basetimestep = 86400.0
    nlayers = 4 #max number of soil layers (ini file setting)
    #TODO: generate layer thickness for each grid cell (input SoilThicknes and vector with UStoreLayerThickness)

    #default parameter values (dict)
    dparams = Dict(
        "Cfmax" => 3.75653 * (timestepsecs / basetimestep),
        "TT" => 0.0,
        "TTM" => 0.0,
        "TTI" => 1.0,
        "WHC" => 0.1,
        "cf_soil" => 0.038,
        "w_soil" => 0.1125 * (timestepsecs / basetimestep),
        "SoilThickness" => 2000.0,
        "InfiltCapSoil" => 100.0,
        "InfiltCapPath" => 10.0,
        "PathFrac" => 0.01,
        "WaterFrac" => 0.0,
        "thetaS" => 0.6,
        "thetaR" => 0.01,
        "AirEntryPressure" => 10.0,
        "KsatVer" => 3000.0 * (timestepsecs / basetimestep),
        "MaxLeakage" => 0.0,
        "c" => 10.0,
        "M" => 300.0,
        "CapScale" => 100.0,
        "rootdistpar" => -500.0,
        "RootingDepth" => 750.0,
        "LAI" => 1.0,
        "et_reftopot" => 1.0,
    )

    nc = NCDataset(staticmaps_path)

    subcatch_2d = "wflow_subcatch" in keys(nc) ? nc["wflow_subcatch"][:] :
        @error("wflow_subcatch not found")
    # indices based on catchment
    inds = Wflow.active_indices(subcatch_2d, missing)
    n = length(inds)

    altitude =
        "wflow_dem" in keys(nc) ? readnetcdf(nc, "wflow_dem", inds, dparams) :
        @error("wflow_dem not found")

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

    waterfrac = readnetcdf(nc, "WaterFrac", inds, dparams)
    pathfrac = readnetcdf(nc, "PathFrac", inds, dparams)

    rootingdepth = readnetcdf(nc, "RootingDepth", inds, dparams)
    rootdistpar = readnetcdf(nc, "rootdistpar", inds, dparams)
    capscale = readnetcdf(nc, "CapScale", inds, dparams)
    et_reftopot = readnetcdf(nc, "et_reftopot", inds, dparams)

    # only read Sl, Swood and Kext if LAI attribute (ini file)
    sl = readnetcdf(nc, "Sl", inds, dparams)
    swood = readnetcdf(nc, "Swood", inds, dparams)
    kext = readnetcdf(nc, "Kext", inds, dparams)
    # otherwise SBM needs EoverR, Cmax and CanopyGapFraction

    #TODO: store c in staticmaps.nc as Array (4,:)
    c = zeros(nlayers, n)
    for i in [0:1:nlayers-1;]
        c[i+1, :] = readnetcdf(nc, string("c_", i), inds, dparams)
    end

    # set in inifile? Also type (monthly, daily, hourly) as part of netcdf variable attribute?
    # in original inifile: LAI=staticmaps/clim/LAI,monthlyclim,1.0,1
    lai_clim = NCDataset(leafarea_path) #TODO:include LAI climatology in update() vertical SBM model

    sbm = Vector{SBM}(undef, n)
    for i = 1:n
        sbm[i] = SBM{nlayers}(
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
            m = m[i],
            hb = hb[i],
            soilthickness = soilthickness[i],
            infiltcappath = infiltcappath[i],
            infiltcapsoil = infiltcapsoil[i],
            maxleakage = maxleakage[i],
            waterfrac = waterfrac[i],
            pathfrac = pathfrac[i],
            altitude = altitude[i],
            rootingdepth = rootingdepth[i],
            rootdistpar = rootdistpar[i],
            capscale = capscale[i],
            et_reftopot = et_reftopot[i],
            sl = sl[i],
            swood = swood[i],
            kext = kext[i],
            c = c[:, i],
            lai = 1.0,
        )

    end

    return sbm

end
