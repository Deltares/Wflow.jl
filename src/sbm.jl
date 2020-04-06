const mv = NaN

# from ini file
# [run]
# timestepsecs = 86400
timestepsecs = 86400.0

Base.@kwdef struct SBM{N}
    Cfmax::Float64
    TT::Float64
    TTM::Float64
    TTI::Float64
    WHC::Float64
    cf_soil::Float64
    w_soil::Float64
    SoilThickness::Float64
    InfiltCapSoil::Float64
    InfiltCapPath::Float64
    PathFrac::Float64
    WaterFrac::Float64
    thetaS::Float64
    thetaR::Float64
    AirEntryPressure::Float64
    KsatVer::Float64
    MaxLeakage::Float64
    c::SVector{N,Float64}
    M::Float64
    f::Float64 = (thetaS - thetaR) / M
    CapScale::Float64
    rootdistpar::Float64
    RootingDepth::Float64
    LAI::Float64
    Sl::Float64
    Kext::Float64
    Swood::Float64
    et_RefToPot::Float64
    Altitude::Float64
    Precipitation::Float64 = mv
    Temperature::Float64 = mv
    PotenEvap::Float64 = mv
    PotTransSoil::Float64 = mv
    Transpiration::Float64 = mv
    ActEvapUStore::Float64 = mv
    ActEvapSat::Float64 = mv
    Interception::Float64 = mv
    ActEvap::Float64 = mv
    ActEvapOpenWaterLand::Float64 = mv
    ActEvapOpenWaterRiver::Float64 = mv
    AvailableForInfiltration::Float64 = mv
    zi::Float64 = mv
    UStoreLayerDepth::SVector{N,Float64} = fill(mv, SVector{N,Float64}) #TODO:define nLayers per grid cell
    UstoreDepth::Float64 = mv
    Transfer::Float64 = mv
    CapFlux::Float64 = mv
    Recharge::Float64 = mv
    SoilWaterCapacity::Float64 = SoilThickness * (thetaS - thetaR)
    SatWaterDepth::Float64 = 0.85 * SoilWaterCapacity
    Snow::Float64 = 0.0
    SnowWater::Float64 = 0.0
    TSoil::Float64 = 10.0
    CanopyStorage::Float64 = 0.0
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

function StateVariables()

    # depends on ini file settings (optional: glaciers, snow, irrigation)
    states = [
        :SatWaterDepth,
        :Snow,
        :Tsoil,
        :UStoreLayerDepth,
        :SnowWater,
        :CanopyStorage,
    ]
    #TODO: (warm) states read from netcdf file or cold state (reinit=1, setting in ini file)

end

"Initial part of the model. Reads model parameters from disk"
function initialize(staticmaps_path, leafarea_path)

    basetimestep = 86400.0
    nLayers = 4 #max number of soil layers (ini file setting)
    #TODO: generate layer thickness for each grid cell (input SoilThicknes and vector with UStoreLayerThickness)

    #default parameter values (dict)
    dParams = Dict(
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
        "et_RefToPot" => 1.0,
    )

    nc = NCDataset(staticmaps_path)

    subcatch_2d = "wflow_subcatch" in keys(nc) ? nc["wflow_subcatch"][:] :
        @error("wflow_subcatch not found")
    # indices based on catchment
    inds = Wflow.active_indices(subcatch_2d, missing)
    n = length(inds)

    Altitude =
        "wflow_dem" in keys(nc) ? readnetcdf(nc, "wflow_dem", inds, dParams) :
        @error("wflow_dem not found")

    # snow parameters (also set in ini file (snow=True or False)?)
    Cfmax = readnetcdf(nc, "Cfmax", inds, dParams)
    TT = readnetcdf(nc, "TT", inds, dParams)
    TTI = readnetcdf(nc, "TTI", inds, dParams)
    TTM = readnetcdf(nc, "TTM", inds, dParams)
    WHC = readnetcdf(nc, "WHC", inds, dParams)
    w_soil = readnetcdf(nc, "w_soil", inds, dParams)
    cf_soil = readnetcdf(nc, "cf_soil", inds, dParams)

    # soil parameters
    thetaS = readnetcdf(nc, "thetaS", inds, dParams)
    thetaR = readnetcdf(nc, "thetaR", inds, dParams)
    KsatVer = readnetcdf(nc, "KsatVer", inds, dParams)
    M = readnetcdf(nc, "M", inds, dParams)
    AirEntryPressure = readnetcdf(nc, "AirEntryPressure", inds, dParams)
    SoilThickness = readnetcdf(nc, "SoilThickness", inds, dParams)
    InfiltCapPath = readnetcdf(nc, "InfiltCapPath", inds, dParams)
    InfiltCapSoil = readnetcdf(nc, "InfiltCapSoil", inds, dParams)
    MaxLeakage = readnetcdf(nc, "MaxLeakage", inds, dParams)

    WaterFrac = readnetcdf(nc, "WaterFrac", inds, dParams)
    PathFrac = readnetcdf(nc, "PathFrac", inds, dParams)

    RootingDepth = readnetcdf(nc, "RootingDepth", inds, dParams)
    rootdistpar = readnetcdf(nc, "rootdistpar", inds, dParams)
    CapScale = readnetcdf(nc, "CapScale", inds, dParams)
    et_RefToPot = readnetcdf(nc, "et_RefToPot", inds, dParams)

    # only read Sl, Swood and Kext if LAI attribute (ini file)
    Sl = readnetcdf(nc, "Sl", inds, dParams)
    Swood = readnetcdf(nc, "Swood", inds, dParams)
    Kext = readnetcdf(nc, "Kext", inds, dParams)
    # otherwise SBM needs EoverR, Cmax and CanopyGapFraction

    #TODO: store c in staticmaps.nc as Array (4,:)
    c = zeros(nLayers, n)
    for i in [0:1:nLayers-1;]
        c[i+1, :] = readnetcdf(nc, string("c_", i), inds, dParams)
    end

    # set in inifile? Also type (monthly, daily, hourly) as part of netcdf variable attribute?
    # in original inifile: LAI=staticmaps/clim/LAI,monthlyclim,1.0,1
    lai_clim = NCDataset(leafarea_path) #TODO:include LAI climatology in update() vertical SBM model

    sbm = Vector{SBM}(undef, n)
    for i = 1:n
        sbm[i] = SBM{nLayers}(
            Cfmax = Cfmax[i],
            TT = TT[i],
            TTI = TTI[i],
            TTM = TTM[i],
            WHC = WHC[i],
            w_soil = w_soil[i],
            cf_soil = cf_soil[i],
            thetaS = thetaS[i],
            thetaR = thetaR[i],
            KsatVer = KsatVer[i],
            M = M[i],
            AirEntryPressure = AirEntryPressure[i],
            SoilThickness = SoilThickness[i],
            InfiltCapPath = InfiltCapPath[i],
            InfiltCapSoil = InfiltCapSoil[i],
            MaxLeakage = MaxLeakage[i],
            WaterFrac = WaterFrac[i],
            PathFrac = PathFrac[i],
            Altitude = Altitude[i],
            RootingDepth = RootingDepth[i],
            rootdistpar = rootdistpar[i],
            CapScale = CapScale[i],
            et_RefToPot = et_RefToPot[i],
            Sl = Sl[i],
            Swood = Swood[i],
            Kext = Kext[i],
            c = c[:, i],
            LAI = 1.0,
        )

    end

    return sbm

end
