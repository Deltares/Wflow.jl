### Soil erosion ###
@get_units @with_kw struct LandSediment{T}
    # number of cells
    n::Int | "-"
    ### Soil erosion part ###
    # length of cells in y direction [m]
    yl::Vector{T} | "m"
    # length of cells in x direction [m]
    xl::Vector{T} | "m"
    # Fraction of river [-]
    riverfrac::Vector{T} | "-"
    # Waterbodies areas [-]
    wbcover::Vector{T} | "-"
    # Depth of overland flow [m]
    h_land::Vector{T} | "m"
    # Canopy interception [mm]
    interception::Vector{T} | "mm"
    # Precipitation [mm]
    precipitation::Vector{T} | "mm"
    # Overland flow [m3/s]
    q_land::Vector{T} | "m3 s-1"
    # Canopy height [m]
    canopyheight::Vector{T} | "m"
    # Canopy gap fraction [mm]
    canopygapfraction::Vector{T} | "mm"
    # Coefficient for EUROSEM rainfall erosion [-]
    erosk::Vector{T} | "-"
    # Exponent for EUROSEM rainfall erosion [-]
    erosspl::Vector{T} | "-"
    # Coefficient for ANSWERS overland flow erosion [-]
    erosov::Vector{T} | "-"
    # Fraction of impervious area per grid cell [-]
    pathfrac::Vector{T} | "-"
    # Land slope [-]
    slope::Vector{T} | "-"
    # USLE crop management factor [-]
    usleC::Vector{T} | "-"
    # USLE soil erodibility factor [-]
    usleK::Vector{T} | "-"
    # Sediment eroded by rainfall [ton]
    sedspl::Vector{T} | "t"
    # Sediment eroded by overland flow [ton]
    sedov::Vector{T} | "t"
    # Total eroded soil [ton]
    soilloss::Vector{T} | "t"
    # Eroded soil per particle class [ton]
    erosclay::Vector{T} | "t"  # clay
    erossilt::Vector{T} | "t"  # silt
    erossand::Vector{T} | "t"  # sand
    erossagg::Vector{T} | "t"  # small aggregates
    eroslagg::Vector{T} | "t"  # large aggregates
    ## Interception related to leaf_area_index climatology ###
    # Specific leaf storage [mm]
    sl::Vector{T} | "mm"
    # Storage woody part of vegetation [mm]
    swood::Vector{T} | "mm"
    # Extinction coefficient [-] (to calculate canopy gap fraction)
    kext::Vector{T} | "-"
    # Leaf area index [m² m⁻²]
    leaf_area_index::Vector{T} | "m2 m-2"
    ### Transport capacity part ###                          
    # Drain length [m]
    dl::Vector{T} | "m"
    # Flow width [m]                          
    width::Vector{T} | "m"
    # Govers transport capacity coefficients [-]
    cGovers::Vector{T} | "-"
    nGovers::Vector{T} | "-"
    # Median particle diameter of the topsoil [mm]
    D50::Vector{T} | "mm"
    # Median diameter per particle class [um]
    dmclay::Vector{T} | "µm"
    dmsilt::Vector{T} | "µm"
    dmsand::Vector{T} | "µm"
    dmsagg::Vector{T} | "µm"
    dmlagg::Vector{T} | "µm"
    # Fraction of the different particle class [-]
    fclay::Vector{T} | "-"
    fsilt::Vector{T} | "-"
    fsand::Vector{T} | "-"
    fsagg::Vector{T} | "-"
    flagg::Vector{T} | "-"
    # Density of sediment [kg/m3]
    rhos::Vector{T} | "kg m-3"
    # Filter with river cells
    rivcell::Vector{T} | "-"
    # Total transport capacity of overland flow [ton]          
    TCsed::Vector{T} | "t"
    # Transport capacity of overland flow per particle class [ton]          
    TCclay::Vector{T} | "t"
    TCsilt::Vector{T} | "t"
    TCsand::Vector{T} | "t"
    TCsagg::Vector{T} | "t"
    TClagg::Vector{T} | "t"

    function LandSediment{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

statevars(::LandSediment) = ()

function initialize_landsed(nc, config, river, riverfrac, xl, yl, inds)
    # Initialize parameters for the soil loss part
    n = length(inds)
    do_river = get(config.model, "runrivermodel", false)::Bool
    # Reservoir / lake
    do_reservoirs = get(config.model, "doreservoir", false)::Bool
    do_lakes = get(config.model, "dolake", false)::Bool
    # Rainfall erosion equation: ["answers", "eurosem"]
    rainerosmethod = get(config.model, "rainerosmethod", "answers")::String
    # Overland flow transport capacity method: ["yalinpart", "govers", "yalin"]
    landtransportmethod = get(config.model, "landtransportmethod", "yalinpart")::String

    altitude =
        ncread(nc, param(config, "input.vertical.altitude"); sel = inds, type = Float64)
    canopyheight = ncread(
        nc,
        param(config, "input.vertical.canopyheight", nothing);
        sel = inds,
        defaults = 3.0,
        type = Float64,
    )
    erosk = ncread(
        nc,
        param(config, "input.vertical.erosk", nothing);
        sel = inds,
        defaults = 0.6,
        type = Float64,
    )
    erosspl = ncread(
        nc,
        param(config, "input.vertical.erosspl", nothing);
        sel = inds,
        defaults = 2.0,
        type = Float64,
    )
    erosov = ncread(
        nc,
        param(config, "input.vertical.erosov", nothing);
        sel = inds,
        defaults = 0.9,
        type = Float64,
    )
    pathfrac = ncread(
        nc,
        param(config, "input.vertical.pathfrac", nothing);
        sel = inds,
        defaults = 0.01,
        type = Float64,
    )
    slope = ncread(
        nc,
        param(config, "input.vertical.slope", nothing);
        sel = inds,
        defaults = 0.01,
        type = Float64,
    )
    usleC = ncread(
        nc,
        param(config, "input.vertical.usleC", nothing);
        sel = inds,
        defaults = 0.01,
        type = Float64,
    )
    usleK = ncread(
        nc,
        param(config, "input.vertical.usleK", nothing);
        sel = inds,
        defaults = 0.1,
        type = Float64,
    )

    cmax, _, canopygapfraction, sl, swood, kext = initialize_canopy(nc, config, inds)

    # Initialise parameters for the transport capacity part
    βₗ = slope
    clamp!(βₗ, 0.00001, Inf)
    ldd_2d = ncread(nc, param(config, "input.ldd"); allow_missing = true)
    ldd = ldd_2d[inds]
    dl = fill(mv, n)
    dw = fill(mv, n)

    for i = 1:n
        dl[i] = detdrainlength(ldd[i], xl[i], yl[i])
        dw[i] = detdrainwidth(ldd[i], xl[i], yl[i])
    end
    dmclay = ncread(
        nc,
        param(config, "input.vertical.dmclay", nothing);
        sel = inds,
        defaults = 2.0,
        type = Float64,
    )
    dmsilt = ncread(
        nc,
        param(config, "input.vertical.dmsilt", nothing);
        sel = inds,
        defaults = 10.0,
        type = Float64,
    )
    dmsand = ncread(
        nc,
        param(config, "input.vertical.dmsand", nothing);
        sel = inds,
        defaults = 200.0,
        type = Float64,
    )
    dmsagg = ncread(
        nc,
        param(config, "input.vertical.dmsagg", nothing);
        sel = inds,
        defaults = 30.0,
        type = Float64,
    )
    dmlagg = ncread(
        nc,
        param(config, "input.vertical.dmlagg", nothing);
        sel = inds,
        defaults = 500.0,
        type = Float64,
    )
    pclay = ncread(
        nc,
        param(config, "input.vertical.pclay", nothing);
        sel = inds,
        defaults = 0.1,
        type = Float64,
    )
    psilt = ncread(
        nc,
        param(config, "input.vertical.psilt", nothing);
        sel = inds,
        defaults = 0.1,
        type = Float64,
    )
    rhos = ncread(
        nc,
        param(config, "input.vertical.rhosed", nothing);
        sel = inds,
        defaults = 2650.0,
        type = Float64,
    )

    ### Initialize transport capacity variables ###
    rivcell = float(river)
    # Percent Sand
    psand = 100 .- pclay .- psilt
    # Govers coefficient for transport capacity
    if landtransportmethod != "yalinpart"
        # Calculation of D50 and fraction of fine and very fine sand (fvfs) from Fooladmand et al, 2006
        psand999 = psand .* ((999 - 25) / (1000 - 25))
        vd50 = log.((1 ./ (0.01 .* (pclay .+ psilt)) .- 1) ./ (1 ./ (0.01 .* pclay) .- 1))
        wd50 =
            log.(
                (1 ./ (0.01 .* (pclay .+ psilt .+ psand999)) .- 1) ./
                (1 ./ (0.01 .* pclay) .- 1),
            )
        ad50 = 1 / log((25 - 1) / (999 - 1))
        bd50 = ad50 ./ log((25 - 1) / 1)
        cd50 = ad50 .* log.(vd50 ./ wd50)
        ud50 = (.-vd50) .^ (1 .- bd50) ./ ((.-wd50) .^ (.-bd50))
        D50 = 1 .+ (-1 ./ ud50 .* log.(1 ./ (1 ./ (0.01 .* pclay) .- 1))) .^ (1 ./ cd50) #[um]
        D50 = D50 ./ 1000 # [mm]
    else
        D50 = fill(mv, n)
    end
    if landtransportmethod == "govers"
        cGovers = ((D50 .* 1000 .+ 5) ./ 0.32) .^ (-0.6)
        nGovers = ((D50 .* 1000 .+ 5) ./ 300) .^ (0.25)
    else
        cGovers = fill(mv, n)
        nGovers = fill(mv, n)
    end
    if do_river || landtransportmethod == "yalinpart"
        # Determine sediment size distribution, estimated from primary particle size distribution (Foster et al., 1980)
        fclay = 0.20 .* pclay ./ 100
        fsilt = 0.13 .* psilt ./ 100
        fsand = 0.01 .* psand .* (1 .- 0.01 .* pclay) .^ (2.4)
        fsagg = 0.28 .* (0.01 .* pclay .- 0.25) .+ 0.5
        for i = 1:n
            if pclay[i] > 50.0
                fsagg[i] = 0.57
            elseif pclay[i] < 25
                fsagg[i] = 2.0 * 0.01 * pclay[i]
            end
        end
        flagg = 1.0 .- fclay .- fsilt .- fsand .- fsagg
    else
        fclay = fill(mv, n)
        fsilt = fill(mv, n)
        fsand = fill(mv, n)
        fsagg = fill(mv, n)
        flagg = fill(mv, n)
    end

    # Reservoir and lakes
    wbcover = fill(0.0, n)
    if do_reservoirs
        rescoverage_2d = ncread(
            nc,
            param(config, "input.vertical.resareas");
            sel = inds,
            type = Float64,
            fill = 0.0,
        )
        wbcover = wbcover .+ rescoverage_2d
    end
    if do_lakes
        lakecoverage_2d = ncread(
            nc,
            param(config, "input.vertical.lakeareas");
            sel = inds,
            type = Float64,
            fill = 0.0,
        )
        wbcover = wbcover .+ lakecoverage_2d
    end

    eros = LandSediment{Float64}(
        n = n,
        yl = yl,
        xl = xl,
        riverfrac = riverfrac,
        wbcover = wbcover,
        ### Soil erosion part ###
        # Forcing
        interception = fill(mv, n),
        h_land = fill(mv, n),
        precipitation = fill(mv, n),
        q_land = fill(mv, n),
        # Parameters
        canopyheight = canopyheight,
        canopygapfraction = canopygapfraction,
        erosk = erosk,
        erosspl = erosspl,
        erosov = erosov,
        pathfrac = pathfrac,
        slope = slope,
        usleC = usleC,
        usleK = usleK,
        # Interception related to climatology (leaf_area_index)
        sl = sl,
        swood = swood,
        kext = kext,
        leaf_area_index = fill(mv, n),
        # Outputs
        sedspl = fill(mv, n),
        sedov = fill(mv, n),
        soilloss = fill(mv, n),
        erosclay = fill(mv, n),
        erossilt = fill(mv, n),
        erossand = fill(mv, n),
        erossagg = fill(mv, n),
        eroslagg = fill(mv, n),
        ### Transport capacity part ###
        # Parameters
        dl = dl,
        width = dw,
        cGovers = cGovers,
        D50 = D50,
        dmclay = dmclay,
        dmsilt = dmsilt,
        dmsand = dmsand,
        dmsagg = dmsagg,
        dmlagg = dmlagg,
        fclay = fclay,
        fsilt = fsilt,
        fsand = fsand,
        fsagg = fsagg,
        flagg = flagg,
        nGovers = nGovers,
        rhos = rhos,
        rivcell = rivcell,
        # Outputs
        TCsed = fill(mv, n),
        TCclay = fill(mv, n),
        TCsilt = fill(mv, n),
        TCsand = fill(mv, n),
        TCsagg = fill(mv, n),
        TClagg = fill(mv, n),
    )

    return eros
end

# Soil erosion
function update_until_ols(eros::LandSediment, config, network)
    # Options from config
    do_lai = haskey(config.input.vertical, "leaf_area_index")
    rainerosmethod = get(config.model, "rainerosmethod", "answers")::String
    Δt = Second(config.timestepsecs)
    ts = Float64(Δt.value)

    for i = 1:eros.n

        ### Splash / Rainfall erosion ###
        # ANSWERS method
        if rainerosmethod == "answers"
            # calculate rainfall intensity [mm/min]
            rintnsty = eros.precipitation[i] / (ts / 60)
            # splash erosion [kg/min]
            sedspl =
                0.108 * eros.usleC[i] * eros.usleK[i] * eros.xl[i] * eros.yl[i] * rintnsty^2
            # [ton/timestep]
            sedspl = sedspl * (ts / 60) * 10^(-3)
        end
        # TODO check eurosem method output values (too high!)
        if rainerosmethod == "eurosem"
            if do_lai
                cmax = eros.sl[i] * eros.leaf_area_index[i] + eros.swood[i]
                canopygapfraction = exp(-eros.kext[i] * eros.leaf_area_index[i])
            end
            # calculate rainfall intensity [mm/h]
            rintnsty = eros.precipitation[i] / (ts / 3600)
            # Kinetic energy of direct throughfall [J/m2/mm]
            # kedir = max(11.87 + 8.73 * log10(max(0.0001, rintnsty)),0.0) #basis used in USLE
            kedir = max(8.95 + 8.44 * log10(max(0.0001, rintnsty)), 0.0) #variant used in most distributed mdoels
            # Kinetic energy of leaf drainage [J/m2/mm]
            pheff = 0.5 * eros.canopyheight[i]
            keleaf = max((15.8 * pheff^0.5) - 5.87, 0.0)

            #Depths of rainfall (total, leaf drianage, direct) [mm]
            rdtot = eros.precipitation[i]
            rdleaf = rdtot * 0.1 * canopygapfraction #stemflow
            rddir = max(rdtot - rdleaf - eros.interception[i], 0.0) #throughfall

            #Total kinetic energy by rainfall [J/m2]
            ketot = (rddir * kedir + rdleaf * keleaf) * 0.001
            # Rainfall / splash erosion [g/m2]
            sedspl = eros.erosk[i] * ketot * exp(-eros.erosspl[i] * eros.h_land[i])
            sedspl = sedspl * eros.xl[i] * eros.yl[i] * 10^(6) # ton/cell
        end

        # Remove the impervious area
        sedspl = sedspl * (1.0 - eros.pathfrac[i])

        ### Overland flow erosion ###
        # ANWERS method
        # Overland flow rate [m2/min]
        qr_land = eros.q_land[i] * 60 / ((eros.xl[i] + eros.yl[i]) / 2)
        # Sine of the slope
        sinslope = sin(atan(eros.slope[i]))

        # Overland flow erosion [kg/min]
        # For a wide range of slope, it is better to use the sine of slope rather than tangeant
        sedov =
            eros.erosov[i] *
            eros.usleC[i] *
            eros.usleK[i] *
            eros.xl[i] *
            eros.yl[i] *
            sinslope *
            qr_land
        # [ton/timestep]
        sedov = sedov * (ts / 60) * 10^(-3)
        # Remove the impervious area
        sedov = sedov * (1.0 - eros.pathfrac[i])

        # Assume no erosion in reservoir/lake areas
        if eros.wbcover[i] > 0
            sedspl = 0.0
            sedov = 0.0
        end

        # Total soil loss [ton/cell/timestep]
        soilloss = sedspl + sedov

        # update the outputs and states
        eros.sedspl[i] = sedspl
        eros.sedov[i] = sedov
        eros.soilloss[i] = soilloss
        # Eroded amount per particle class
        eros.erosclay[i] = soilloss * eros.fclay[i]
        eros.erossilt[i] = soilloss * eros.fsilt[i]
        eros.erossand[i] = soilloss * eros.fsand[i]
        eros.erossagg[i] = soilloss * eros.fsagg[i]
        eros.eroslagg[i] = soilloss * eros.flagg[i]
    end

end

### Sediment transport capacity in overland flow ###
function update_until_oltransport(ols::LandSediment, config, network)

    do_river = get(config.model, "runrivermodel", false)::Bool
    tcmethod = get(config.model, "landtransportmethod", "yalinpart")::String
    Δt = Second(config.timestepsecs)
    ts = Float64(Δt.value)

    for i = 1:ols.n
        sinslope = sin(atan(ols.slope[i]))

        if do_river != true
            # Total transport capacity without particle differentiation
            if tcmethod == "govers"
                # Transport capacity from govers 1990
                # Unit stream power
                if ols.h_land[i] > 0.0
                    velocity = ols.q_land[i] / (ols.dl[i] * ols.h_land[i])
                else
                    velocity = 0.0
                end
                omega = 10 * sinslope * 100 * velocity #cm/s
                if omega > 0.4
                    TCf = ols.cGovers[i] * (omega - 0.4)^(ols.nGovers[i]) * 2650 #kg/m3
                else
                    TCf = 0.0
                end
                TC = TCf * ols.q_land[i] * ts * 10^(-3) #[ton/cell]
            end

            if tcmethod == "yalin"
                # Transport capacity from Yalin without particle differentiation
                delta = max(
                    (
                        ols.h_land[i] * sinslope /
                        (ols.D50[i] * 0.001 * (ols.rhos[i] / 1000 - 1)) / 0.06 - 1
                    ),
                    0.0,
                )
                alphay = delta * 2.45 / (0.001 * ols.rhos[i])^0.4 * 0.06^(0.5)
                if ols.q_land[i] > 0.0 && alphay != 0.0
                    TC = (
                        ols.dl[i] / ols.q_land[i] *
                        (ols.rhos[i] - 1000) *
                        ols.D50[i] *
                        0.001 *
                        (9.81 * ols.h_land[i] * sinslope) *
                        0.635 *
                        delta *
                        (1 - log(1 + alphay) / (alphay))
                    ) # [kg/m3]
                    TC = TC * ols.q_land[i] * ts * 10^(-3) #[ton]
                else
                    TC = 0.0
                end
            end

            # Filter TC land for river cells (0 in order for sediment from land to stop when entering the river)
            if ols.rivcell[i] == 1.0
                TC = 0.0
            end

            # Set particle TC to 0
            TCclay = 0.0
            TCsilt = 0.0
            TCsand = 0.0
            TCsagg = 0.0
            TClagg = 0.0
        end

        if do_river || tcmethod == "yalinpart"
            # Transport capacity from Yalin with particle differentiation
            # Delta parameter of Yalin for each particle class
            delta = ols.h_land[i] * sinslope / (10^(-6) * (ols.rhos[i] / 1000 - 1) / 0.06)
            dclay = max(1 / ols.dmclay[i] * delta - 1, 0.0)
            dsilt = max(1 / ols.dmsilt[i] * delta - 1, 0.0)
            dsand = max(1 / ols.dmsand[i] * delta - 1, 0.0)
            dsagg = max(1 / ols.dmsagg[i] * delta - 1, 0.0)
            dlagg = max(1 / ols.dmlagg[i] * delta - 1, 0.0)
            # Total transportability
            dtot = dclay + dsilt + dsand + dsagg + dlagg

            # Yalin transport capacity of overland flow for each particle class
            if ols.q_land[i] > 0.0
                TCa =
                    ols.dl[i] / ols.q_land[i] *
                    (ols.rhos[i] - 1000) *
                    10^(-6) *
                    (9.81 * ols.h_land[i] * sinslope)
            else
                TCa = 0.0
            end
            TCb = 2.45 / (ols.rhos[i] / 1000)^0.4 * 0.06^0.5
            if dtot != 0.0 && dclay != 0.0
                TCclay =
                    TCa * ols.dmclay[i] * dclay / dtot *
                    0.635 *
                    dclay *
                    (1 - log(1 + dclay * TCb) / dclay * TCb) # [kg/m3]
                TCclay = TCclay * ols.q_land[i] * ts * 10^(-3) # [ton]
            else
                TCclay = 0.0
            end
            if dtot != 0.0 && dsilt != 0.0
                TCsilt =
                    TCa * ols.dmsilt[i] * dsilt / dtot *
                    0.635 *
                    dsilt *
                    (1 - log(1 + dsilt * TCb) / dsilt * TCb) # [kg/m3]
                TCsilt = TCsilt * ols.q_land[i] * ts * 10^(-3) # [ton]
            else
                TCsilt = 0.0
            end
            if dtot != 0.0 && dsand != 0.0
                TCsand =
                    TCa * ols.dmsand[i] * dsand / dtot *
                    0.635 *
                    dsand *
                    (1 - log(1 + dsand * TCb) / dsand * TCb) # [kg/m3]
                TCsand = TCsand * ols.q_land[i] * ts * 10^(-3) # [ton]
            else
                TCsand = 0.0
            end
            if dtot != 0.0 && dsagg != 0.0
                TCsagg =
                    TCa * ols.dmsagg[i] * dsagg / dtot *
                    0.635 *
                    dsagg *
                    (1 - log(1 + dsagg * TCb) / dsagg * TCb) # [kg/m3]
                TCsagg = TCsagg * ols.q_land[i] * ts * 10^(-3) # [ton]
            else
                TCsagg = 0.0
            end
            if dtot != 0.0 && dlagg != 0.0
                TClagg =
                    TCa * ols.dmlagg[i] * dlagg / dtot *
                    0.635 *
                    dlagg *
                    (1 - log(1 + dlagg * TCb) / dlagg * TCb) # [kg/m3]
                TClagg = TClagg * ols.q_land[i] * ts * 10^(-3) # [ton]
            else
                TClagg = 0.0
            end

            # Assume that ols all reach the river in reservoir/lake areas (very high TC)
            if ols.wbcover[i] > 0
                TC = 10^9
                TCclay = 10^9
                TCsilt = 10^9
                TCsand = 10^9
                TCsagg = 10^9
                TClagg = 10^9
            end

            # Filter TC land for river cells (0 in order for sediment from land to stop when entering the river)
            if ols.rivcell[i] == 1.0
                TCclay = 0.0
                TCsilt = 0.0
                TCsand = 0.0
                TCsagg = 0.0
                TClagg = 0.0
            end

            # Set total TC to 0
            TC = 0.0

        end

        # update the outputs and states
        ols.TCsed[i] = TC
        ols.TCclay[i] = TCclay
        ols.TCsilt[i] = TCsilt
        ols.TCsand[i] = TCsand
        ols.TCsagg[i] = TCsagg
        ols.TClagg[i] = TClagg

    end

end

### Sediment transport in overland flow ###
@get_units @with_kw struct OverlandFlowSediment{T}
    # number of cells
    n::Int | "-"
    # Filter with river cells
    rivcell::Vector{T} | "-"
    # Forcing for the river
    h_riv::Vector{T} | "m"
    q_riv::Vector{T} | "m3 s-1"
    # Total eroded soil [ton]
    soilloss::Vector{T} | "t"
    # Eroded soil per particle class [ton]
    erosclay::Vector{T} | "t"
    erossilt::Vector{T} | "t"
    erossand::Vector{T} | "t"
    erossagg::Vector{T} | "t"
    eroslagg::Vector{T} | "t"
    # Total transport capacity of overland flow [ton]          
    TCsed::Vector{T} | "t"
    # Transport capacity of overland flow per particle class [ton]          
    TCclay::Vector{T} | "t"
    TCsilt::Vector{T} | "t"
    TCsand::Vector{T} | "t"
    TCsagg::Vector{T} | "t"
    TClagg::Vector{T} | "t"
    # Sediment flux in overland flow [ton]
    olsed::Vector{T} | "t"
    olclay::Vector{T} | "t"
    olsilt::Vector{T} | "t"
    olsand::Vector{T} | "t"
    olsagg::Vector{T} | "t"
    ollagg::Vector{T} | "t"
    # Sediment reaching the river with overland flow [ton]
    inlandsed::Vector{T} | "t"
    inlandclay::Vector{T} | "t"
    inlandsilt::Vector{T} | "t"
    inlandsand::Vector{T} | "t"
    inlandsagg::Vector{T} | "t"
    inlandlagg::Vector{T} | "t"


    function OverlandFlowSediment{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

statevars(::OverlandFlowSediment) = ()

function partial_update!(inland, rivcell, eroded)
    no_erosion = zero(eltype(eroded))
    for i in eachindex(inland)
        inland[i] = rivcell[i] == 1 ? eroded[i] : no_erosion
    end
    return inland
end

function update(ols::OverlandFlowSediment, network, config)
    do_river = get(config.model, "runrivermodel", false)::Bool
    tcmethod = get(config.model, "landtransportmethod", "yalinpart")::String
    zeroarr = fill(0.0, ols.n)

    # transport sediment down the network
    if do_river || tcmethod == "yalinpart"
        # clay
        accucapacityflux!(ols.olclay, ols.erosclay, network, ols.TCclay)
        partial_update!(ols.inlandclay, ols.rivcell, ols.erosclay)
        # silt
        accucapacityflux!(ols.olsilt, ols.erossilt, network, ols.TCsilt)
        partial_update!(ols.inlandsilt, ols.rivcell, ols.erossilt)
        # sand
        accucapacityflux!(ols.olsand, ols.erossand, network, ols.TCsand)
        partial_update!(ols.inlandsand, ols.rivcell, ols.erossand)
        # small aggregates
        accucapacityflux!(ols.olsagg, ols.erossagg, network, ols.TCsagg)
        partial_update!(ols.inlandsagg, ols.rivcell, ols.erossagg)
        # large aggregates
        accucapacityflux!(ols.ollagg, ols.eroslagg, network, ols.TClagg)
        partial_update!(ols.inlandlagg, ols.rivcell, ols.eroslagg)

        # total sediment, all particle classes
        ols.olsed .= ols.olclay .+ ols.olsilt .+ ols.olsand .+ ols.olsagg .+ ols.ollagg
        ols.inlandsed .=
            ols.inlandclay .+ ols.inlandsilt .+ ols.inlandsand .+ ols.inlandsagg .+
            ols.inlandlagg
    else
        accucapacityflux!(ols.olsed, ols.soilloss, network, ols.TCsed)
        ols.inlandsed .= ifelse.(ols.rivcell .== 1, ols.soilloss, zeroarr)
    end
end

### River transport and processes ###
@get_units @with_kw struct RiverSediment{T}
    # number of cells
    n::Int | "-"
    # Timestep [s]
    Δt::T | "s"
    # River geometry (slope [-], length [m], width [m])
    sl::Vector{T} | "m"
    dl::Vector{T} | "m"
    width::Vector{T} | "m"
    # Sediment mean diameter in the river bed [mm]
    d50::Vector{T}
    # Particle mean diameter [mm]
    dmclay::Vector{T}
    dmsilt::Vector{T}
    dmsand::Vector{T}
    dmsagg::Vector{T}
    dmlagg::Vector{T}
    dmgrav::Vector{T}
    # River bed and bank particle fraction composition [-]
    fclayriv::Vector{T} | "-"
    fsiltriv::Vector{T} | "-"
    fsandriv::Vector{T} | "-"
    fgravriv::Vector{T} | "-"
    # Sediment mean diameter for Engelund and Hansen transport equation [mm]
    d50engelund::Vector{T}
    # Parameters for Bagnold transport equation
    cbagnold::Vector{T} | "-"
    ebagnold::Vector{T} | "-"
    # Parameters for Kodatie transport equation
    ak::Vector{T} | "-"
    bk::Vector{T} | "-"
    ck::Vector{T} | "-"
    dk::Vector{T} | "-"
    # Critical bed and bank shear stress [N/m2]
    TCrbank::Vector{T} | "N m-2"
    TCrbed::Vector{T} | "N m-2"
    # Bed and bank erodibilities [m3/N.s]
    kdbank::Vector{T} | "m3 N-1 s-1"
    kdbed::Vector{T} | "m3 N-1 s-1"
    # Sediment density [kg/m3]
    rhos::Vector{T} | "kg m-3"
    # River water level [m]
    h_riv::Vector{T} | "m"
    # River discharge [m3/s]
    q_riv::Vector{T} | "m3 s-1"
    # Sediment input from land erosion [ton]
    inlandclay::Vector{T} | "t"
    inlandsilt::Vector{T} | "t"
    inlandsand::Vector{T} | "t"
    inlandsagg::Vector{T} | "t"
    inlandlagg::Vector{T} | "t"
    # Sediment / particle left in the cell [ton]
    sedload::Vector{T} | "t"
    clayload::Vector{T} | "t"
    siltload::Vector{T} | "t"
    sandload::Vector{T} | "t"
    saggload::Vector{T} | "t"
    laggload::Vector{T} | "t"
    gravload::Vector{T} | "t"
    # Sediment / particle stored on the river bed after deposition [ton]
    sedstore::Vector{T} | "t"
    claystore::Vector{T} | "t"
    siltstore::Vector{T} | "t"
    sandstore::Vector{T} | "t"
    saggstore::Vector{T} | "t"
    laggstore::Vector{T} | "t"
    gravstore::Vector{T} | "t"
    # Sediment / particle flux [ton]
    outsed::Vector{T} | "t"
    outclay::Vector{T} | "t"
    outsilt::Vector{T} | "t"
    outsand::Vector{T} | "t"
    outsagg::Vector{T} | "t"
    outlagg::Vector{T} | "t"
    outgrav::Vector{T} | "t"
    # Sediment concentrations [kg/m3]
    Sedconc::Vector{T} | "kg m-3"
    SSconc::Vector{T} | "kg m-3"
    Bedconc::Vector{T} | "kg m-3"
    # Reservoir and lakes
    wbcover::Vector{T} | "-"
    wblocs::Vector{T} | "-"
    wbarea::Vector{T} | "m2"

    # function RiverSediment{T}(args...) where {T}
    #     equal_size_vectors(args)
    #     return new(args...)
    # end
end

statevars(::RiverSediment) = (
    :clayload,
    :siltload,
    :sandload,
    :saggload,
    :laggload,
    :gravload,
    :claystore,
    :siltstore,
    :sandstore,
    :saggstore,
    :laggstore,
    :gravstore,
    :outclay,
    :outsilt,
    :outsand,
    :outsagg,
    :outlagg,
    :outgrav,
)

function initialize_riversed(nc, config, riverwidth, riverlength, inds_riv)
    # Initialize river parameters
    nriv = length(inds_riv)
    # River flow transport capacity method: ["bagnold", "engelund", "yang", "kodatie", "molinas"]
    tcmethodriv = get(config.model, "rivtransportmethod", "bagnold")::String
    Δt = Second(config.timestepsecs)
    # Reservoir / lakes
    do_reservoirs = get(config.model, "doreservoir", false)::Bool
    do_lakes = get(config.model, "dolake", false)::Bool
    wbcover = fill(0.0, nriv)
    wblocs = fill(0.0, nriv)
    wbarea = fill(0.0, nriv)

    if do_reservoirs
        reslocs = ncread(
            nc,
            param(config, "input.lateral.river.reslocs");
            sel = inds_riv,
            type = Float64,
            fill = 0,
        )
        rescoverage_2d = ncread(
            nc,
            param(config, "input.lateral.river.resareas");
            sel = inds_riv,
            type = Float64,
            fill = 0,
        )
        resarea = ncread(
            nc,
            param(config, "input.lateral.river.resarea");
            sel = inds_riv,
            type = Float64,
            fill = 0,
        )

        wbcover = wbcover .+ rescoverage_2d
        wblocs = wblocs .+ reslocs
        wbarea = wbarea .+ resarea
    end

    if do_lakes
        lakelocs = ncread(
            nc,
            param(config, "input.lateral.river.lakelocs");
            sel = inds_riv,
            type = Float64,
            fill = 0,
        )
        lakecoverage_2d = ncread(
            nc,
            param(config, "input.lateral.river.lakeareas");
            sel = inds_riv,
            type = Float64,
            fill = 0,
        )
        lakearea = ncread(
            nc,
            param(config, "input.lateral.river.lakearea");
            sel = inds_riv,
            type = Float64,
            fill = 0,
        )

        wbcover = wbcover .+ lakecoverage_2d
        wblocs = wblocs .+ lakelocs
        wbarea = wbarea .+ lakearea
    end

    riverslope = ncread(
        nc,
        param(config, "input.lateral.river.slope");
        sel = inds_riv,
        type = Float64,
    )
    clamp!(riverslope, 0.00001, Inf)
    rhos = ncread(
        nc,
        param(config, "input.lateral.river.rhosed", nothing);
        sel = inds_riv,
        defaults = 2650.0,
        type = Float64,
    )
    dmclay = ncread(
        nc,
        param(config, "input.lateral.river.dmclay", nothing);
        sel = inds_riv,
        defaults = 2.0,
        type = Float64,
    )
    dmsilt = ncread(
        nc,
        param(config, "input.lateral.river.dmsilt", nothing);
        sel = inds_riv,
        defaults = 10.0,
        type = Float64,
    )
    dmsand = ncread(
        nc,
        param(config, "input.lateral.river.dmsand", nothing);
        sel = inds_riv,
        defaults = 200.0,
        type = Float64,
    )
    dmsagg = ncread(
        nc,
        param(config, "input.lateral.river.dmsagg", nothing);
        sel = inds_riv,
        defaults = 30.0,
        type = Float64,
    )
    dmlagg = ncread(
        nc,
        param(config, "input.lateral.river.dmlagg", nothing);
        sel = inds_riv,
        defaults = 500.0,
        type = Float64,
    )
    dmgrav = ncread(
        nc,
        param(config, "input.lateral.river.dmgrav", nothing);
        sel = inds_riv,
        defaults = 2000.0,
        type = Float64,
    )
    fclayriv = ncread(
        nc,
        param(config, "input.lateral.river.fclayriv");
        sel = inds_riv,
        type = Float64,
    )
    fsiltriv = ncread(
        nc,
        param(config, "input.lateral.river.fsiltriv");
        sel = inds_riv,
        type = Float64,
    )
    fsandriv = ncread(
        nc,
        param(config, "input.lateral.river.fsandriv");
        sel = inds_riv,
        type = Float64,
    )
    fgravriv = ncread(
        nc,
        param(config, "input.lateral.river.fgravriv");
        sel = inds_riv,
        type = Float64,
    )
    d50riv =
        ncread(nc, param(config, "input.lateral.river.d50"); sel = inds_riv, type = Float64)
    d50engelund = ncread(
        nc,
        param(config, "input.lateral.river.d50engelund");
        sel = inds_riv,
        type = Float64,
    )
    cbagnold = ncread(
        nc,
        param(config, "input.lateral.river.cbagnold");
        sel = inds_riv,
        type = Float64,
    )
    ebagnold = ncread(
        nc,
        param(config, "input.lateral.river.ebagnold");
        sel = inds_riv,
        type = Float64,
    )

    # Initialisation of parameters for Kodatie transport capacity
    ak = fill(0.0, nriv)
    bk = fill(0.0, nriv)
    ck = fill(0.0, nriv)
    dk = fill(0.0, nriv)
    if tcmethodriv == "kodatie"
        for i = 1:nriv
            if d50riv[i] <= 0.05
                ak[i] = 281.4
                bk[i] = 2.622
                ck[i] = 0.182
                dk[i] = 0.0
            elseif d50riv[i] <= 0.25
                ak[i] = 2829.6
                bk[i] = 3.646
                ck[i] = 0.406
                dk[i] = 0.412
            elseif d50riv[i] <= 2.0
                ak[i] = 2123.4
                bk[i] = 3.3
                ck[i] = 0.468
                dk[i] = 0.613
            else
                ak[i] = 431884.8
                bk[i] = 1.0
                ck[i] = 1.0
                dk[i] = 2.0
            end
        end
    end
    # Initialisation of parameters for river erosion
    # Bed and Bank from Shields diagram, Da Silva & Yalin (2017)
    E_ = (2.65 - 1) * 9.81
    E = (E_ .* (d50riv .* 10^(-3)) .^ 3 ./ 10^(-12)) .^ 0.33
    TCrbed = (
        E_ .* d50riv .* (
            0.13 .* E .^ (-0.392) .* exp.(-0.015 .* E .^ 2) .+
            0.045 .* (1 .- exp.(-0.068 .* E))
        )
    )
    TCrbank = TCrbed
    # kd from Hanson & Simon 2001
    kdbank = 0.2 .* TCrbank .^ (-0.5) .* 10^(-6)
    kdbed = 0.2 .* TCrbed .^ (-0.5) .* 10^(-6)

    rs = RiverSediment(
        n = nriv,
        Δt = Float64(Δt.value),
        # Parameters
        sl = riverslope,
        dl = riverlength,
        width = riverwidth,
        dmclay = dmclay,
        dmsilt = dmsilt,
        dmsand = dmsand,
        dmsagg = dmsagg,
        dmlagg = dmlagg,
        dmgrav = dmgrav,
        fclayriv = fclayriv,
        fsiltriv = fsiltriv,
        fsandriv = fsandriv,
        fgravriv = fgravriv,
        d50 = d50riv,
        d50engelund = d50engelund,
        cbagnold = cbagnold,
        ebagnold = ebagnold,
        ak = ak,
        bk = bk,
        ck = ck,
        dk = dk,
        kdbank = kdbank,
        kdbed = kdbed,
        TCrbank = TCrbank,
        TCrbed = TCrbed,
        rhos = rhos,
        # Forcing
        h_riv = fill(mv, nriv),
        q_riv = fill(mv, nriv),
        # Input from land
        inlandclay = fill(0.0, nriv),
        inlandsilt = fill(0.0, nriv),
        inlandsand = fill(0.0, nriv),
        inlandsagg = fill(0.0, nriv),
        inlandlagg = fill(0.0, nriv),
        # Outputs and states
        sedload = fill(0.0, nriv),
        clayload = fill(0.0, nriv),
        siltload = fill(0.0, nriv),
        sandload = fill(0.0, nriv),
        saggload = fill(0.0, nriv),
        laggload = fill(0.0, nriv),
        gravload = fill(0.0, nriv),
        sedstore = fill(0.0, nriv),
        claystore = fill(0.0, nriv),
        siltstore = fill(0.0, nriv),
        sandstore = fill(0.0, nriv),
        saggstore = fill(0.0, nriv),
        laggstore = fill(0.0, nriv),
        gravstore = fill(0.0, nriv),
        outsed = fill(0.0, nriv),
        outclay = fill(0.0, nriv),
        outsilt = fill(0.0, nriv),
        outsand = fill(0.0, nriv),
        outsagg = fill(0.0, nriv),
        outlagg = fill(0.0, nriv),
        outgrav = fill(0.0, nriv),
        Sedconc = fill(0.0, nriv),
        SSconc = fill(0.0, nriv),
        Bedconc = fill(0.0, nriv),
        # Reservoir / lake
        wbcover = wbcover,
        wblocs = wblocs,
        wbarea = wbarea,
    )

    return rs
end

function update(rs::RiverSediment, network, config)
    @unpack graph, order = network
    tcmethod = get(config.model, "rivtransportmethod", "bagnold")::String

    # River sediment loads are separated into different particle class.
    # Clay, silt and sand can both come from land, resuspension or river channel erosion.
    # Small and large aggregates only come from land erosion or resuspension.
    # Gravel only comes from resuspension or river channel erosion.

    for v in order
        ### Sediment input in the cell (left from previous timestep + from land + from upstream outflux) ###
        upstream_nodes = inneighbors(graph, v)

        inrivclay = 0.0
        inrivsilt = 0.0
        inrivsand = 0.0
        inrivsagg = 0.0
        inrivlagg = 0.0
        inrivgrav = 0.0
        if !isempty(upstream_nodes)
            for i in upstream_nodes
                if rs.outclay[i] >= 0.0 # avoid NaN from upstream non-river cells
                    inrivclay += rs.outclay[i]
                    inrivsilt += rs.outsilt[i]
                    inrivsand += rs.outsand[i]
                    inrivsagg += rs.outsagg[i]
                    inrivlagg += rs.outlagg[i]
                    inrivgrav += rs.outgrav[i]
                end
            end
        end

        inclay = rs.clayload[v] + rs.inlandclay[v] + inrivclay
        insilt = rs.siltload[v] + rs.inlandsilt[v] + inrivsilt
        insand = rs.sandload[v] + rs.inlandsand[v] + inrivsand
        insagg = rs.saggload[v] + rs.inlandsagg[v] + inrivsagg
        inlagg = rs.laggload[v] + rs.inlandlagg[v] + inrivlagg
        ingrav = rs.gravload[v] + inrivgrav

        insed = inclay + insilt + insand + insagg + inlagg + ingrav

        ### Transport capacity of the flow ###
        # Hydraulic radius of the river [m] (rectangular channel)
        hydrad = rs.h_riv[v] * rs.width[v] / (rs.width[v] + 2 * rs.h_riv[v])

        # Engelund and Hansen transport formula
        if tcmethod == "engelund"
            vmean =
                ifelse(rs.h_riv[v] > 0.0, rs.q_riv[v] / (rs.width[v] * rs.h_riv[v]), 0.0)
            vshear = (9.81 * hydrad * rs.sl[v])^0.5
            # Concentration by weight
            cw = ifelse(
                hydrad > 0.0,
                (
                    rs.rhos[v] / 1000 * 0.5 * vmean * vshear^3 /
                    ((rs.rhos[v] / 1000 - 1)^2 * 9.81^2 * rs.d50engelund[v] * hydrad)
                ),
                0.0,
            )
            cw = min(1.0, cw)
            # Transport capacity [tons/m3]
            maxsed = max(cw / (cw + (1 - cw) * rs.rhos[v] / 1000) * rs.rhos[v] / 1000, 0.0)
        elseif tcmethod == "bagnold"
            maxsed =
                rs.cbagnold[v] * (rs.q_riv[v] / (rs.h_riv[v] * rs.width[v]))^rs.ebagnold[v]
        elseif tcmethod == "kodatie"
            vmean =
                ifelse(rs.h_riv[v] > 0.0, rs.q_riv[v] / (rs.width[v] * rs.h_riv[v]), 0.0)
            maxsed = rs.ak[v] * vmean^rs.bk[v] * rs.h_riv[v]^rs.ck[v] * rs.sl[v]^rs.dk[v]
            # Transport capacity [tons/m3]
            maxsed =
                ifelse(rs.q_riv[v] > 0.0, maxsed * rs.width[v] / (rs.q_riv[v] * rs.Δt), 0.0)
        elseif tcmethod == "yang"
            ws = 411 * rs.d50[v]^2 / 3600
            vshear = (9.81 * hydrad * rs.sl[v])^0.5
            var1 = vshear * rs.d50[v] / 1000 / (1.16 * 10^(-6))
            var2 = ws * rs.d50[v] / 1000 / (1.16 * 10^(-6))
            vcr = min(
                0.0,
                ifelse(var1 >= 70.0, 2.05 * ws, ws * (2.5 / (log10(var1) - 0.06) + 0.66)),
            )
            # Sand equation
            if (rs.width[v] * rs.h_riv[v]) >= vcr && rs.d50[v] < 2.0
                logcppm = (
                    5.435 - 0.286 * log10(var2) - 0.457 * log10(vshear / ws) + 1.799 -
                    0.409 * log10(var2) -
                    0.314 *
                    log10(vshear / ws) *
                    log10(
                        (rs.q_riv[v] / (rs.width[v] * rs.h_riv[v]) - vcr) * rs.sl[v] / ws,
                    )
                )
                # Gravel equation
            elseif (rs.width[v] * rs.h_riv[v]) >= vcr && rs.d50[v] < 2.0
                logcppm = (
                    6.681 - 0.633 * log10(var2) - 4.816 * log10(vshear / ws) + 2.784 -
                    0.305 * log10(var2) -
                    0.282 *
                    log10(vshear / ws) *
                    log10(
                        (rs.q_riv[v] / (rs.width[v] * rs.h_riv[v]) - vcr) * rs.sl[v] / ws,
                    )
                )
            else
                logcppm = 0.0
            end
            # Sediment concentration by weight
            cw = 10^(logcppm) * 10^(-6)
            # Transport capacity [ton/m3]
            maxsed = max(cw / (cw + (1 - cw) * rs.rhos[v] / 1000) * rs.rhos[v] / 1000, 0.0)
        elseif tcmethod == "molinas"
            ws = 411 * rs.d50[v]^2 / 3600
            vmean =
                ifelse(rs.h_riv[v] > 0.0, rs.q_riv[v] / (rs.width[v] * rs.h_riv[v]), 0.0)
            if rs.h_riv[v] > 0.0
                psi = (
                    vmean^3 / (
                        (rs.rhos[v] / 1000 - 1) *
                        9.81 *
                        rs.h_riv[v] *
                        ws *
                        log10(1000 * rs.h_riv[v] / rs.d50[v])^2
                    )
                )
            else
                psi = 0.0
            end
            # Concentration by weight
            cw = 1430 * (0.86 + psi^0.5) * psi^1.5 / (0.016 + psi) * 10^(-6)
            # Transport capacity [ton/m3]
            maxsed = max(cw / (cw + (1 - cw) * rs.rhos[v] / 1000) * rs.rhos[v] / 1000, 0.0)
        end

        # 1285 g/L: boundary between streamflow and debris flow (Costa, 1988)
        maxsed = min(maxsed, 1.285)
        # Transport capacity [ton]
        maxsed = maxsed * (rs.h_riv[v] * rs.width[v] * rs.dl[v] + rs.q_riv[v] * rs.Δt)

        ### River erosion ###
        # Erosion only if the load is below the transport capacity of the flow.
        sedex = max(maxsed - insed, 0.0)
        # No erosion in lake and reservoir cells
        if rs.wbcover[v] > 0.0
            sedex = 0.0
        end
        # Bed and bank are eroded only if the previously deposited material is not enough
        rs.sedstore[v] =
            rs.claystore[v] +
            rs.siltstore[v] +
            rs.sandstore[v] +
            rs.saggstore[v] +
            rs.laggstore[v] +
            rs.gravstore[v]
        if sedex > 0.0 && sedex > rs.sedstore[v]
            # Effective sediment needed fom river bed and bank erosion [ton]
            effsedex = sedex - rs.sedstore[v]

            # Repartition of the effective shear stress between the bank and the bed from Knight et al. 1984 [%]
            SFbank = ifelse(
                rs.h_riv[v] > 0.0,
                exp(-3.23 * log10(rs.width[v] / rs.h_riv[v] + 3) + 6.146),
                0.0,
            )
            # Effective shear stress on river bed and banks [N/m2]
            TEffbank = ifelse(
                rs.h_riv[v] > 0.0,
                1000 * 9.81 * hydrad * rs.sl[v] * SFbank / 100 *
                (1 + rs.width[v] / (2 * rs.h_riv[v])),
                0.0,
            )
            TEffbed =
                1000 *
                9.81 *
                hydrad *
                rs.sl[v] *
                (1 - SFbank / 100) *
                (1 + 2 * rs.h_riv[v] / rs.width[v])
            # Potential erosion rates of the bed and bank [t/cell/timestep] 
            #(assuming only one bank is eroding)
            Tex = max(TEffbank - rs.TCrbank[v], 0.0)
            # 1.4 is bank default bulk density
            ERbank = max(0.0, rs.kdbank[v] * Tex * rs.dl[v] * rs.h_riv[v] * 1.4 * rs.Δt)
            # 1.5 is bed default bulk density
            ERbed = max(
                0.0,
                rs.kdbed[v] *
                (TEffbed - rs.TCrbed[v]) *
                rs.dl[v] *
                rs.width[v] *
                1.5 *
                rs.Δt,
            )
            # Relative potential erosion rates of the bed and the bank [-]
            RTEbank = ifelse(ERbank + ERbed > 0.0, ERbank / (ERbank + ERbed), 0.0)
            RTEbed = 1.0 - RTEbank

            # Bank erosion (difference between effective and potential erosion) [ton]
            sedbank = ifelse(effsedex * RTEbank <= ERbank, effsedex * RTEbank, ERbank)
            claybank = rs.fclayriv[v] * sedbank
            siltbank = rs.fsiltriv[v] * sedbank
            sandbank = rs.fsandriv[v] * sedbank
            gravbank = rs.fgravriv[v] * sedbank

            # Bed erosion [ton]
            sedbed = ifelse(effsedex * RTEbed <= ERbed, effsedex * RTEbed, ERbed)
            claybed = rs.fclayriv[v] * sedbed
            siltbed = rs.fsiltriv[v] * sedbed
            sandbed = rs.fsandriv[v] * sedbed
            gravbed = rs.fgravriv[v] * sedbed
        else
            sedbank = 0.0
            claybank = 0.0
            siltbank = 0.0
            sandbank = 0.0
            gravbank = 0.0

            sedbed = 0.0
            claybed = 0.0
            siltbed = 0.0
            sandbed = 0.0
            gravbed = 0.0
        end

        # Erosion/degradation of the previously deposited sediment (from clay to gravel) [ton]
        if sedex > 0.0
            degstoreclay = ifelse(rs.claystore[v] >= sedex, sedex, rs.claystore[v])
            rs.claystore[v] = rs.claystore[v] - degstoreclay
            #Update amount of sediment that need to be degraded
            sedex = sedex - degstoreclay
            degstoresilt = ifelse(rs.siltstore[v] >= sedex, sedex, rs.siltstore[v])
            rs.siltstore[v] = rs.siltstore[v] - degstoresilt
            sedex = max(0.0, sedex - degstoresilt)
            degstoresagg = ifelse(rs.saggstore[v] >= sedex, sedex, rs.saggstore[v])
            rs.saggstore[v] = rs.saggstore[v] - degstoresagg
            sedex = max(0.0, sedex - degstoresagg)
            degstoresand = ifelse(rs.sandstore[v] >= sedex, sedex, rs.sandstore[v])
            rs.sandstore[v] = rs.sandstore[v] - degstoresand
            sedex = max(0.0, sedex - degstoresand)
            degstorelagg = ifelse(rs.laggstore[v] >= sedex, sedex, rs.laggstore[v])
            rs.laggstore[v] = rs.laggstore[v] - degstorelagg
            sedex = max(0.0, sedex - degstorelagg)
            degstoregrav = ifelse(rs.gravstore[v] >= sedex, sedex, rs.gravstore[v])
            rs.gravstore[v] = rs.gravstore[v] - degstoregrav
            sedex = max(0.0, sedex - degstoregrav)
            degstoresed =
                degstoreclay +
                degstoresilt +
                degstoresagg +
                degstoresand +
                degstorelagg +
                degstoregrav
        else
            degstoreclay = 0.0
            degstoresilt = 0.0
            degstoresagg = 0.0
            degstoresand = 0.0
            degstorelagg = 0.0
            degstoregrav = 0.0
            degstoresed = 0.0
        end

        # Sum all erosion sources per particle class
        erodsed = sedbank + sedbed + degstoresed
        erodclay = claybank + claybed + degstoreclay
        erodsilt = siltbank + siltbed + degstoresilt
        erodsand = sandbank + sandbed + degstoresand
        erodsagg = degstoresagg
        erodlagg = degstorelagg
        erodgrav = gravbank + gravbed + degstoregrav


        ### Deposition / settling ###
        # Fractions of deposited particles in river cells from the Einstein formula [-]
        # Particle fall velocity [m/s] from Stokes
        xs = ifelse(rs.q_riv[v] > 0.0, 1.055 * rs.dl[v] / (rs.q_riv[v] / rs.width[v]), 0.0)
        xclay = min(1.0, 1.0 - 1.0 / exp(xs * (411.0 * (rs.dmclay[v] / 1000)^2 / 3600)))
        xsilt = min(1.0, 1.0 - 1.0 / exp(xs * (411.0 * (rs.dmsilt[v] / 1000)^2 / 3600)))
        xsand = min(1.0, 1.0 - 1.0 / exp(xs * (411.0 * (rs.dmsand[v] / 1000)^2 / 3600)))
        xsagg = min(1.0, 1.0 - 1.0 / exp(xs * (411.0 * (rs.dmsagg[v] / 1000)^2 / 3600)))
        xlagg = min(1.0, 1.0 - 1.0 / exp(xs * (411.0 * (rs.dmlagg[v] / 1000)^2 / 3600)))
        xgrav = min(1.0, 1.0 - 1.0 / exp(xs * (411.0 * (rs.dmgrav[v] / 1000)^2 / 3600)))

        # Sediment deposited in the channel [ton]
        # From natural settling with Einstein formula (density controlled).
        settclay = xclay * (inclay + erodclay)
        settsilt = xsilt * (insilt + erodsilt)
        settsand = xsand * (insand + erodsand)
        settsagg = xsagg * (insagg + erodsagg)
        settlagg = xlagg * (inlagg + erodlagg)
        settgrav = xgrav * (ingrav + erodgrav)

        # Sediment deposited in the channel (from gravel to clay) [ton]
        # From transport capacity exceedance (insed > maxsed)
        insedex = min(insed - maxsed, 0.0)
        if insedex > 0.0
            depgrav = ifelse(ingrav >= insedex, insedex, ingrav)
            insedex = max(insedex - depgrav, 0.0)
            deplagg = ifelse(inlagg >= insedex, insedex, inlagg)
            insedex = max(insedex - deplagg, 0.0)
            depsand = ifelse(insand >= insedex, insedex, insand)
            insedex = max(insedex - depsand, 0.0)
            depsagg = ifelse(insagg >= insedex, insedex, insagg)
            insedex = max(insedex - depsagg, 0.0)
            depsilt = ifelse(insilt >= insedex, insedex, insilt)
            insedex = max(insedex - depsilt, 0.0)
            depclay = ifelse(inclay >= insedex, insedex, inclay)
            insedex = max(insedex - depclay, 0.0)
        else
            depclay = settclay
            depsilt = settsilt
            depsand = settsand
            depsagg = settsagg
            deplagg = settlagg
            depgrav = settgrav
        end

        # No deposition in regular lake and reservoir cells, only at the outlet
        # Deposition in lake/reservoir from Camp 1945
        if rs.wbcover[v] > 0.0 && rs.wblocs[v] > 0.0
            # Compute deposition
            vcres = rs.q_riv[v] / rs.wbarea[v]
            DCres = 411 / 3600 / vcres
            depclay = (inclay + erodclay) * min(1.0, (DCres * (rs.dmclay[v] / 1000)^2))
            depsilt = (insilt + erodsilt) * min(1.0, (DCres * (rs.dmsilt[v] / 1000)^2))
            depsand = (insand + erodsand) * min(1.0, (DCres * (rs.dmsand[v] / 1000)^2))
            depsagg = (insagg + erodsagg) * min(1.0, (DCres * (rs.dmsagg[v] / 1000)^2))
            deplagg = (inlagg + erodlagg) * min(1.0, (DCres * (rs.dmlagg[v] / 1000)^2))
            depgrav = (ingrav + erodgrav) * min(1.0, (DCres * (rs.dmgrav[v] / 1000)^2))
        elseif rs.wbcover[v] > 0.0
            depsed = 0.0
            depclay = 0.0
            depsilt = 0.0
            depsand = 0.0
            depsagg = 0.0
            deplagg = 0.0
            depgrav = 0.0
        end

        depsed = depclay + depsilt + depsand + depsagg + deplagg + depgrav

        # Update the river deposited sediment storage
        rs.sedstore[v] = rs.sedstore[v] + depsed
        rs.claystore[v] = rs.claystore[v] + depclay
        rs.siltstore[v] = rs.siltstore[v] + depsilt
        rs.sandstore[v] = rs.sandstore[v] + depsand
        rs.saggstore[v] = rs.saggstore[v] + depsagg
        rs.laggstore[v] = rs.laggstore[v] + deplagg
        rs.gravstore[v] = rs.gravstore[v] + depgrav

        ### Ouput loads ###
        # Sediment transported out of the cell during the timestep [ton]
        # 0 in case all sediment are deposited in the cell
        # Reduce the fraction so that there is still some sediment staying in the river cell
        fwaterout = min(rs.q_riv[v] * rs.Δt / (rs.h_riv[v] * rs.width[v] * rs.dl[v]), 1.0)
        rs.outsed[v] = fwaterout * (insed + erodsed - depsed)
        rs.outclay[v] = fwaterout * (inclay + erodclay - depclay)
        rs.outsilt[v] = fwaterout * (insilt + erodsilt - depsilt)
        rs.outsand[v] = fwaterout * (insand + erodsand - depsand)
        rs.outsagg[v] = fwaterout * (insagg + erodsagg - depsagg)
        rs.outlagg[v] = fwaterout * (inlagg + erodlagg - deplagg)
        rs.outgrav[v] = fwaterout * (ingrav + erodgrav - depgrav)

        ### Mass balance ###
        # Sediment left in the cell [ton]
        rs.sedload[v] = insed + erodsed - depsed - rs.outsed[v]
        rs.clayload[v] = inclay + erodclay - depclay - rs.outclay[v]
        rs.siltload[v] = insilt + erodsilt - depsilt - rs.outsilt[v]
        rs.sandload[v] = insand + erodsand - depsand - rs.outsand[v]
        rs.saggload[v] = insagg + erodsagg - depsagg - rs.outsagg[v]
        rs.gravload[v] = ingrav + erodgrav - depgrav - rs.outgrav[v]

        ### Concentrations and suspended sediments ###
        # Conversion from load [ton] to concentration for rivers [mg/L]
        toconc = ifelse(rs.q_riv[v] > 0.0, 10^6 / (rs.q_riv[v] * rs.Δt), 0.0)
        rs.Sedconc[v] = rs.outsed[v] * toconc

        # Differentiation of bed and suspended load using Rouse number for suspension
        # threshold diameter between bed load and mixed load using Rouse number
        dbedf = 10^3 * (2.5 * 3600 * 0.41 / 411 * (9.81 * rs.h_riv[v] * rs.sl[v])^0.5)^0.5
        # threshold diameter between suspended load and mixed load using Rouse number
        dsuspf = 10^3 * (1.2 * 3600 * 0.41 / 411 * (9.81 * rs.h_riv[v] * rs.sl[v])^0.5)^0.5
        # Rouse with diameter
        SSclay = ifelse(
            rs.dmclay[v] <= dsuspf,
            rs.outclay[v],
            ifelse(rs.dmclay[v] <= dbedf, rs.outclay[v] / 2, 0.0),
        )
        SSsilt = ifelse(
            rs.dmsilt[v] <= dsuspf,
            rs.outsilt[v],
            ifelse(rs.dmsilt[v] <= dbedf, rs.outsilt[v] / 2, 0.0),
        )
        SSsagg = ifelse(
            rs.dmsagg[v] <= dsuspf,
            rs.outsagg[v],
            ifelse(rs.dmsagg[v] <= dbedf, rs.outsagg[v] / 2, 0.0),
        )
        SSsand = ifelse(
            rs.dmsand[v] <= dsuspf,
            rs.outsand[v],
            ifelse(rs.dmsand[v] <= dbedf, rs.outsand[v] / 2, 0.0),
        )
        SSlagg = ifelse(
            rs.dmlagg[v] <= dsuspf,
            rs.outlagg[v],
            ifelse(rs.dmlagg[v] <= dbedf, rs.outlagg[v] / 2, 0.0),
        )
        SSgrav = ifelse(
            rs.dmgrav[v] <= dsuspf,
            rs.outgrav[v],
            ifelse(rs.dmgrav[v] <= dbedf, rs.outgrav[v] / 2, 0.0),
        )

        SS = SSclay + SSsilt + SSsagg + SSsand + SSlagg + SSgrav
        Bed = rs.outsed[v] - SS

        rs.SSconc[v] = SS * toconc
        rs.Bedconc[v] = Bed * toconc

    end

end
