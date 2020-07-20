"""
    initialize_sbm_model(config::Config)

Initial part of the SBM model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_sbm_model(config::Config)

    # unpack the paths to the NetCDF files
    tomldir = dirname(config)
    staticmaps_path = joinpath(tomldir, config.input.staticmaps)
    forcing_path = joinpath(tomldir, config.input.forcing)
    instate_path = joinpath(tomldir, config.state.input.path)
    output_path = joinpath(tomldir, config.output.path)

    Δt = Second(config.timestepsecs)
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
        "N" => 0.072,
        "NRiver" => 0.036,
    )

    sizeinmetres = Bool(get(config.model, "sizeinmetres", 0))
    reinit = Bool(get(config.model, "reinit", 0))
    config_thicknesslayers = get(config.model, "thicknesslayers", 0.0)
    if length(config_thicknesslayers) > 0
        thicknesslayers = SVector(Tuple(push!(Float64.(config_thicknesslayers), mv)))
        sumlayers = pushfirst(cumsum(thicknesslayers), 0.0)
        maxlayers = length(thicknesslayers) # max number of soil layers
    else
        maxlayers = 1
    end
    do_reservoirs = Bool(get(config.model, "reservoirs", 0))
    do_lakes = Bool(get(config.model, "lakes", 0))
    do_snow = Bool(get(config.model, "modelsnow", 0))

    nc = NCDataset(staticmaps_path)
    dims = dimnames(nc["wflow_subcatch"])

    # There is no need to permute the dimensions of the data, since the active indices are
    # correctly calculated in both ways.
    # The dimension order only needs to be known for interpreting the LDD directions
    # and creating the coordinate maps.
    trsp = dims[2] in ("y", "lat")

    subcatch_2d = ncread(nc, "wflow_subcatch", allow_missing = true)
    # indices based on catchment
    inds = Wflow.active_indices(subcatch_2d, missing)
    n = length(inds)

    altitude = ncread(nc, "wflow_dem"; sel = inds, type = Float64)
    river_2d = ncread(nc, "wflow_river"; type = Bool, fill = false)
    river = river_2d[inds]
    riverwidth_2d = ncread(nc, "wflow_riverwidth"; type = Float64, fill = 0)
    riverwidth = riverwidth_2d[inds]
    riverlength_2d = ncread(nc, "wflow_riverlength"; type = Float64, fill = 0)
    riverlength = riverlength_2d[inds]

    # read x, y coordinates and calculate cell length [m]
    y_nc = "y" in keys(nc.dim) ? ncread(nc, "y") : ncread(nc, "lat")
    x_nc = "x" in keys(nc.dim) ? ncread(nc, "x") : ncread(nc, "lon")
    if trsp
        y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    else
        y = repeat(y_nc, outer = (1, length(x_nc)))[inds]
    end
    cellength = abs(mean(diff(x_nc)))


    cfmax = ncread(nc, "Cfmax"; sel = inds, defaults = dparams, type = Float64)
    tt = ncread(nc, "TT"; sel = inds, defaults = dparams, type = Float64)
    tti = ncread(nc, "TTI"; sel = inds, defaults = dparams, type = Float64)
    ttm = ncread(nc, "TTM"; sel = inds, defaults = dparams, type = Float64)
    whc = ncread(nc, "WHC"; sel = inds, defaults = dparams, type = Float64)
    w_soil = ncread(nc, "w_soil"; sel = inds, defaults = dparams, type = Float64)
    cf_soil = ncread(nc, "cf_soil"; sel = inds, defaults = dparams, type = Float64)


    # soil parameters
    θₛ = ncread(nc, "thetaS"; sel = inds, defaults = dparams, type = Float64)
    θᵣ = ncread(nc, "thetaR"; sel = inds, defaults = dparams, type = Float64)
    kv₀ = ncread(nc, "KsatVer"; sel = inds, defaults = dparams, type = Float64)
    m = ncread(nc, "M"; sel = inds, defaults = dparams, type = Float64)
    hb = ncread(nc, "AirEntryPressure"; sel = inds, defaults = dparams, type = Float64)
    soilthickness =
        ncread(nc, "SoilThickness"; sel = inds, defaults = dparams, type = Float64)
    infiltcappath =
        ncread(nc, "InfiltCapPath"; sel = inds, defaults = dparams, type = Float64)
    infiltcapsoil =
        ncread(nc, "InfiltCapSoil"; sel = inds, defaults = dparams, type = Float64)
    maxleakage = ncread(nc, "MaxLeakage"; sel = inds, defaults = dparams, type = Float64)
    # TODO: store c, kvfrac in staticmaps.nc start at index 1
    c = fill(dparams["c"], (maxlayers, n))
    kvfrac = fill(dparams["KsatVerFrac"], (maxlayers, n))
    for i in [0:1:maxlayers-1;]
        if string("c_", i) in keys(nc)
            c[i+1, :] = ncread(nc, string("c_", i); sel = inds, type = Float64)
        else
            @warn(string("c_", i, " not found, set to default value ", dparams["c"]))
        end
        if string("KsatVerFrac_", i) in keys(nc)
            kvfrac[i+1, :] =
                ncread(nc, string("KsatVerFrac_", i); sel = inds, type = Float64)
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
    waterfrac = ncread(nc, "WaterFrac"; sel = inds, defaults = dparams, type = Float64)
    pathfrac = ncread(nc, "PathFrac"; sel = inds, defaults = dparams, type = Float64)

    # vegetation parameters
    rootingdepth =
        ncread(nc, "RootingDepth"; sel = inds, defaults = dparams, type = Float64)
    rootdistpar = ncread(nc, "rootdistpar"; sel = inds, defaults = dparams, type = Float64)
    capscale = ncread(nc, "CapScale"; sel = inds, defaults = dparams, type = Float64)
    et_reftopot = ncread(nc, "et_reftopot"; sel = inds, defaults = dparams, type = Float64)

    # if lai climatology provided use sl, swood and kext to calculate cmax, e_r and canopygapfraction
    # TODO replace by something else
    if isnothing(true)
        # cmax, e_r, canopygapfraction only required when lai climatoly not provided
        cmax = ncread(nc, "Cmax"; sel = inds, defaults = dparams, type = Float64)
        e_r = ncread(nc, "EoverR"; sel = inds, defaults = dparams, type = Float64)
        canopygapfraction =
            ncread(nc, "CanopyGapFraction"; sel = inds, defaults = dparams, type = Float64)
        sl = fill(mv, n)
        swood = fill(mv,n)
        kext = fill(mv,n)
    else
        # TODO confirm if lai climatology is present in the NetCDF
        sl = ncread(nc, "Sl"; sel = inds, defaults = dparams, type = Float64)
        swood = ncread(nc, "Swood"; sel = inds, defaults = dparams, type = Float64)
        kext = ncread(nc, "Kext"; sel = inds, defaults = dparams, type = Float64)
        cmax = fill(mv, n)
        e_r = fill(mv, n)
        canopygapfraction = fill(mv, n)
    end


    # these are filled in the loop below
    # TODO see if we can replace this approach
    nlayers = zeros(Int, n)
    act_thickl = zeros(Float64, maxlayers, n)
    s_layers = zeros(Float64, maxlayers + 1, n)
    xl = fill(mv, n)
    yl = fill(mv, n)
    riverfrac = fill(mv, n)

    for i = 1:n
        xl[i] = sizeinmetres ? cellength : lattometres(y[i])[1] * cellength
        yl[i] = sizeinmetres ? cellength : lattometres(y[i])[2] * cellength
        riverfrac[i] =
            river[i] ? min((riverlength[i] * riverwidth[i]) / (xl[i] * yl[i]), 1.0) : 0.0

        if length(config_thicknesslayers) > 0
            act_thickl_, nlayers_ =
                set_layerthickness(soilthickness[i], sumlayers, thicknesslayers)
            s_layers_ = pushfirst(cumsum(act_thickl_), 0.0)

            nlayers[i] = nlayers_
            act_thickl[:, i] = act_thickl_
            s_layers[:, i] = s_layers_
        else
            nlayers[i] = 1
            act_thickl[:, i] = SVector(soilthickness[i])
            s_layers[:, i] = pushfirst(cumsum(SVector(soilthickness[i])), 0.0)
        end
    end

    # needed for derived parameters below
    act_thickl = svectorscopy(act_thickl, Val{maxlayers}())
    θₑ = θₛ .- θᵣ
    soilwatercapacity = soilthickness .* θₑ
    satwaterdepth = 0.85 .* soilwatercapacity # cold state value for satwaterdepth

    # copied to array of sarray below
    vwc = fill(mv, maxlayers, n)
    vwc_perc = fill(mv, maxlayers, n)

    # states sbm concept
    if do_snow
        statenames =("satwaterdepth", "snow", "tsoil", "ustorelayerdepth", "snowwater", "canopystorage")
    else
        statenames = ("satwaterdepth", "ustorelayerdepth", "canopystorage")
    end

    sbm = SBM{Float64, maxlayers, maxlayers+1}(
        maxlayers = maxlayers,
        n = n,
        nlayers = nlayers,
        yl = yl,
        xl = xl,
        riverfrac = riverfrac,
        θₛ = θₛ,
        θᵣ = θᵣ,
        θₑ = θₑ,
        kv₀ = kv₀,
        kvfrac = svectorscopy(kvfrac, Val{maxlayers}()),
        m = m,
        hb = hb,
        soilthickness = soilthickness,
        act_thickl = act_thickl,
        sumlayers = svectorscopy(s_layers, Val{maxlayers + 1}()),
        infiltcappath = infiltcappath,
        infiltcapsoil = infiltcapsoil,
        maxleakage = maxleakage,
        waterfrac = max.(waterfrac .- riverfrac, 0.0),
        pathfrac = pathfrac,
        altitude = altitude,
        rootingdepth = rootingdepth,
        rootdistpar = rootdistpar,
        capscale = capscale,
        et_reftopot = et_reftopot,
        c = svectorscopy(c, Val{maxlayers}()),
        stemflow = fill(mv, n),
        throughfall = fill(mv, n),
        f = θₑ ./ m,
        ustorelayerdepth = act_thickl .* 0.0,
        satwaterdepth = satwaterdepth,
        zi = max.(0.0, soilthickness .- satwaterdepth ./ θₑ),
        soilwatercapacity = soilwatercapacity,
        canopystorage = fill(0.0, n),
        cmax = cmax,
        canopygapfraction = canopygapfraction,
        e_r = e_r,
        precipitation = fill(mv, n),
        temperature = fill(mv, n),
        potevap = fill(mv, n),
        pottrans_soil = fill(mv, n),
        transpiration = fill(mv, n),
        ae_ustore = fill(mv, n),
        ae_sat = fill(mv, n),
        interception = fill(mv, n),
        soilevap = fill(mv, n),
        actevapsat = fill(mv, n),
        actevap = fill(mv, n),
        runoff_river = fill(mv, n),
        runoff_land = fill(mv, n),
        ae_openw_l = fill(mv, n),
        ae_openw_r = fill(mv, n),
        avail_forinfilt = fill(mv, n),
        actinfilt = fill(mv, n),
        actinfiltsoil = fill(mv, n),
        actinfiltpath = fill(mv, n),
        infiltexcess = fill(mv, n),
        excesswater = fill(mv, n),
        exfiltsatwater = fill(mv, n),
        exfiltustore = fill(mv, n),
        excesswatersoil = fill(mv, n),
        excesswaterpath = fill(mv, n),
        runoff = fill(mv, n),
        vwc = svectorscopy(vwc, Val{maxlayers}()),
        vwc_perc = svectorscopy(vwc_perc, Val{maxlayers}()),
        rootstore = fill(mv, n),
        vwc_root = fill(mv, n),
        vwc_percroot = fill(mv, n),
        ustoredepth = fill(mv, n),
        transfer = fill(mv, n),
        capflux = fill(mv, n),
        recharge = fill(mv, n),
        # snow parameters
        cfmax = cfmax,
        tt = tt,
        tti = tti,
        ttm = ttm,
        whc = whc,
        w_soil = w_soil,
        cf_soil = cf_soil,
        snow = fill(0.0, n),
        snowwater = fill(0.0, n),
        rainfallplusmelt = fill(mv, n),
        tsoil = fill(10.0, n),
        # Interception related to climatology (LAI)
        sl = sl,
        swood = swood,
        kext = kext,
        lai = fill(mv, n),
        )

    inds_riv = filter(i -> !isequal(river_2d[i], 0), inds)
    # reservoirs
    pits = zeros(Bool, n)
    if do_reservoirs
        # read only reservoir data if reservoirs true
        reslocs_2d = ncread(nc, "wflow_reservoirlocs"; type = Int, fill = 0)
        # allow reservoirs only in river cells
        inds_res = filter(i -> reslocs_2d[inds][i] > 0 && isequal(river[i], 1), 1:n)
        resdemand = ncread(nc, "ResDemand"; sel = inds_riv, type = Float64, fill = 0)
        resmaxrelease =
            ncread(nc, "ResMaxRelease"; sel = inds_riv, type = Float64, fill = 0)
        resmaxvolume = ncread(nc, "ResMaxVolume"; sel = inds_riv, type = Float64, fill = 0)
        resarea = ncread(nc, "ResSimpleArea"; sel = inds_riv, type = Float64, fill = 0)
        res_targetfullfrac =
            ncread(nc, "ResTargetFullFrac"; sel = inds_riv, type = Float64, fill = 0)
        res_targetminfrac =
            ncread(nc, "ResTargetMinFrac"; sel = inds_riv, type = Float64, fill = 0)
        reslocs = reslocs_2d[inds_riv]

        # for surface water routing reservoir locations are considered pits in the flow network
        # all upstream flow goes to the river and flows into the reservoir
        pits[inds_res] .= true

        reservoirs = [
            reslocs[i] == 0 ? missing :
                SimpleReservoir{Float64}(
                demand = resdemand[i],
                maxrelease = resmaxrelease[i],
                maxvolume = resmaxvolume[i],
                area = resarea[i],
                targetfullfrac = res_targetfullfrac[i],
                targetminfrac = res_targetminfrac[i],
            ) for i = 1:length(reslocs)
        ]
        statenames =  (statenames..., "volume_reservoir")
    end

    # lakes
    if do_lakes
        # read only lake data if lakes true
        lakelocs_2d = ncread(nc, "wflow_lakelocs"; type = Int, fill = 0)
        # allow lakes only in river cells
        inds_lakes = filter(i -> lakelocs_2d[inds][i] > 0 && isequal(river[i], 1), 1:n)
        lakearea = ncread(nc, "LakeArea"; sel = inds_riv, type = Float64, fill = 0)
        lake_b = ncread(nc, "Lake_b"; sel = inds_riv, type = Float64, fill = 0)
        lake_e = ncread(nc, "Lake_e"; sel = inds_riv, type = Float64, fill = 0)
        lake_threshold =
            ncread(nc, "LakeThreshold"; sel = inds_riv, type = Float64, fill = 0)
        linked_lakelocs = ncread(nc, "LinkedLakeLocs"; sel = inds_riv, type = Int, fill = 0)
        lake_storfunc = ncread(nc, "LakeStorFunc"; sel = inds_riv, type = Int, fill = 0)
        lake_outflowfunc =
            ncread(nc, "LakeOutflowFunc"; sel = inds_riv, type = Int, fill = 0)
        lake_avglevel = ncread(nc, "LakeAvgLevel"; sel = inds_riv, type = Float64, fill = 0)
        lakelocs = lakelocs_2d[inds_riv]

        # for surface water routing lake locations are considered pits in the flow network
        # all upstream flow goes to the river and flows into the lake
        pits[inds_lakes] .= true

        n_lakes = length(lakelocs)
        lakes = Vector{Union{Missing,NaturalLake{Float64}}}(missing, n_lakes)
        for i = 1:length(n_lakes)
            if lakelocs[i] > 0
                lakes[i] = NaturalLake(
                    loc_id = lakelocs[i],
                    lowerlake_ind = linked_lakelocs[i] > 0 ? i : 0,
                    area = lakearea[i],
                    threshold = lake_threshold[i],
                    storfunc = lake_storfunc[i],
                    outflowfunc = lake_outflowfunc[i],
                    b = lake_b[i],
                    e = lake_e[i],
                    avg_waterlevel = lake_avglevel[i],
                    sh = lake_storfunc[i] == 2 ?
                             CSV.read(("Lake_SH_", string(lakelocs[i])), type = Float64) :
                             DataFrame(),
                    hq = lake_outflowfunc[i] == 1 ?
                             CSV.read(("Lake_HQ_", string(lakelocs[i])), type = Float64) :
                             DataFrame(),
                )

                if lake_outflowfunc[i] == 3 && lake_storfunc[i] != 1
                    @warn("For the modified pulse approach (LakeOutflowFunc = 3) the LakeStorFunc should be 1")
                end
            end
        end
        statenames = (statenames..., "waterlevel_lake")
    end
    # lateral part sbm
    khfrac = ncread(nc, "KsatHorFrac"; sel = inds, defaults = dparams, type = Float64)
    βₗ = ncread(nc, "Slope"; sel = inds, type = Float64)
    clamp!(βₗ, 0.00001, Inf)
    ldd_2d = ncread(nc, "wflow_ldd", allow_missing = true)
    ldd = ldd_2d[inds]
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
        wb_pit = pits,
    )

    n_land = ncread(nc, "N"; sel = inds, defaults = dparams, type = Float64)

    olf = SurfaceFlow(
        sl = βₗ,
        n = n_land,
        dl = dl,
        Δt = Float64(Δt.value),
        width = dw,
        wb_pit = pits,
    )

    pcr_dir = trsp ? permute_indices(Wflow.pcrdir) : Wflow.pcrdir
    dag = flowgraph(ldd, inds, pcr_dir)


    riverslope = ncread(nc, "RiverSlope"; sel = inds_riv, type = Float64)
    clamp!(riverslope, 0.00001, Inf)
    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]
    n_river = ncread(nc, "N_River"; sel = inds_riv, defaults = dparams, type = Float64)
    ldd_riv = ldd_2d[inds_riv]
    dag_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)
    nr = length(inds_riv)

    rf = SurfaceFlow(
        sl = riverslope,
        n = n_river,
        dl = riverlength,
        Δt = Float64(Δt.value),
        width = riverwidth,
        reservoir = do_reservoirs ? reservoirs : fill(missing, nr),
        lake = do_lakes ? lakes : fill(missing, nr),
        rivercells = river,
    )

    statenames = (statenames...,"ssf", "q_river", "h_river", "q_land", "h_land")

    reader = prepare_reader(forcing_path, staticmaps_path, inds, inds_riv, config)

    modelmap = (vertical = sbm, subsurface = ssf, land = olf, river = rf)
    writer = prepare_writer(config, reader, output_path, modelmap, maxlayers)

    model = Model(
        (land = dag, river = dag_riv),
        (subsurface = ssf, land = olf, river = rf),
        sbm,
        Clock(config.starttime, 1, Δt),
        reader,
        writer,
    )

    if reinit
        ds = NCDataset(instate_path)
        for state in statenames
            if "layer" in dimnames(ds[state])
                var = Float64.(ds[state][inds,:,1])
                model = set(model, paramap[state], [SVector{ds.dim["layer"]}(var[i,:]) for i=1:n])
            else
                var = Float64.(nomissing(ds[state][inds,1],mv))
                model = set(model, paramap[state], var)
            end
        end
    end

    # make sure the forcing is already loaded
    # it's fine to run twice, and may help catching errors earlier
    update_forcing!(model)
    update_cyclic!(model)
    return model
end

function update(
    model,
    config,
    toposort_land,
    toposort_river,
    frac_toriver,
    index_river,
    nl,
    nr,
)
    @unpack lateral, vertical, network, clock = model

    update_forcing!(model)
    update_cyclic!(model)

    update_until_snow(vertical, config)

    if Bool(get(config.model, "masswasting", 0))
        snowflux_frac =
            min.(0.5, lateral.land.sl ./ 5.67) .* min.(1.0, vertical.snow ./ 10000.0)
        maxflux = snowflux_frac .* vertical.snow
        vertical.snow .=
            accucapacityflux(network.land, toposort_land, vertical.snow, maxflux)
    end

    update_until_recharge(vertical, config)
    lateral.subsurface.recharge .= vertical.recharge
    lateral.subsurface.recharge .*= lateral.subsurface.dl
    lateral.subsurface.zi .= vertical.zi

    update(
        lateral.subsurface,
        network.land,
        toposort_land,
        frac_toriver,
        lateral.river.rivercells,
    )

    update_after_lateralflow(
        vertical,
        lateral.subsurface.zi,
        lateral.subsurface.exfiltwater,
    )

    lateral.land.qlat .=
        (vertical.runoff .* vertical.xl .* vertical.yl .* 0.001) ./ 86400.0 ./
        lateral.land.dl

    update(
        lateral.land,
        network.land,
        toposort_land,
        nl,
        frac_toriver = frac_toriver,
        river = lateral.river.rivercells,
        do_iter = true,
    )

    lateral.river.qlat .=
        (
            lateral.subsurface.to_river[index_river] ./ 1.0e9 ./ lateral.river.Δt .+
            lateral.land.to_river[index_river]
        ) ./ lateral.river.dl

    update(
        lateral.river,
        network.river,
        toposort_river,
        nr,
        do_iter = true,
        doy = dayofyear(clock.time),
    )

    write_output(model, model.writer)

    # update the clock
    clock.iteration += 1
    clock.time += clock.Δt

    return model
end
