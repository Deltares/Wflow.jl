### Soil erosion ###
Base.@kwdef struct LandSed{T}
    # number of cells
    n::Int
    ### Soil erosion part ###
    # length of cells in y direction [m]
    yl::Vector{T}
    # length of cells in x direction [m]
    xl::Vector{T}
    # Fraction of river [-]
    riverfrac::Vector{T}
    # Depth of overland flow [m]
    h_land::Vector{T} = fill(mv, n)
    # Canopy interception [mm]
    interception::Vector{T} = fill(mv, n)
    # Precipitation [mm]
    precipitation::Vector{T} = fill(mv, n)
    # Overland flow [m3/s]
    q_land::Vector{T} = fill(mv, n)
    # Canopy height [m]
    canopyheight::Vector{T}
    # Canopy gap fraction [mm]
    canopygapfraction::Vector{T}
    # Coefficient for EUROSEM rainfall erosion [-]
    erosk::Vector{T}
    # Exponent for EUROSEM rainfall erosion [-]
    erosspl::Vector{T}
    # Coefficient for ANSWERS overland flow erosion [-]
    erosov::Vector{T}
    # Fraction of impervious area per grid cell [-]
    pathfrac::Vector{T}
    # Land slope [-]
    slope::Vector{T}
    # USLE crop management factor [-]
    usleC::Vector{T}
    # USLE soil erodibility factor [-]
    usleK::Vector{T}
    # Sediment eroded by rainfall [ton]
    sedspl::Vector{T} = fill(mv, n)
    # Sediment eroded by overland flow [ton]
    sedov::Vector{T} = fill(mv, n)
    # Total eroded soil [ton]
    soilloss::Vector{T} = fill(mv, n)
    # Eroded soil per particle class [ton]
    erosclay::Vector{T} = fill(mv, n)
    erossilt::Vector{T} = fill(mv, n)
    erossand::Vector{T} = fill(mv, n)
    erossagg::Vector{T} = fill(mv, n)
    eroslagg::Vector{T} = fill(mv, n)
    ## Interception related to leaf_area_index climatology ###
    # Specific leaf storage [mm]
    sl::Vector{T} = fill(mv, n)
    # Storage woody part of vegetation [mm]
    swood::Vector{T} = fill(mv, n)
    # Extinction coefficient [-] (to calculate canopy gap fraction)
    kext::Vector{T} = fill(mv, n)
    # Leaf area index [m² m⁻²]
    leaf_area_index::Vector{T} = fill(mv, n)
    ### Transport capacity part ###                          
    # Drain length [m]
    dl::Vector{T}
    # Flow width [m]                          
    width::Vector{T}
    # Govers transport capacity coefficients [-]
    cGovers::Vector{T}
    nGovers::Vector{T}
    # Median particle diameter of the topsoil [mm]
    D50::Vector{T}
    # Median diameter per particle class [um]
    dmclay::Vector{T}
    dmsilt::Vector{T}
    dmsand::Vector{T}
    dmsagg::Vector{T}
    dmlagg::Vector{T}
    # Fraction of the different particle class [-]
    fclay::Vector{T}
    fsilt::Vector{T}
    fsand::Vector{T}
    fsagg::Vector{T}
    flagg::Vector{T}  
    # Density of sediment [kg/m3]
    rhos::Vector{T}  
    # Filter with river cells
    rivcell::Vector{T}  
    # Total transport capacity of overland flow [ton]          
    TCsed::Vector{T} = fill(mv, n)
    # Transport capacity of overland flow per particle class [ton]          
    TCclay::Vector{T} = fill(mv, n)
    TCsilt::Vector{T} = fill(mv, n)
    TCsand::Vector{T} = fill(mv, n)
    TCsagg::Vector{T} = fill(mv, n)
    TClagg::Vector{T} = fill(mv, n)

    function LandSed{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

statevars(::LandSed) = ()

# Soil erosion
function update_until_ols(eros::LandSed, config)
 
    # # start dummy variables (should be generated from model reader and from Config.jl TOML)
    do_lai = haskey(config.input.vertical, "leaf_area_index")
    rainerosmethod = get(config.model, "rainerosmethod", "answers")
    #precipitation = 3.0
    #q_land = 0.01
    Δt = Second(config.timestepsecs)
    ts = Float64(Δt.value)
    #basetimestep = Second(Day(1))
    # end dummpy variables

    for i = 1:eros.n
                
        ### Splash / Rainfall erosion ###
        # ANSWERS method
        if rainerosmethod == "answers"
            # calculate rainfall intensity [mm/min]
            rintnsty = eros.precipitation[i] / (ts/60)
            # splash erosion [kg/min]
            sedspl = 0.108 * eros.usleC[i] * eros.usleK[i] * eros.xl[i] * eros.yl[i] * rintnsty^2
            # [ton/timestep]
            sedspl = sedspl * (ts/60) * 10^(-3)
        end
        # TODO check eurosem method output values (too high!)
        if rainerosmethod == "eurosem"
            if do_lai
                cmax = eros.sl[i] * eros.leaf_area_index[i] + eros.swood[i]
                canopygapfraction = exp(-eros.kext[i] * eros.leaf_area_index[i])
            end
            # calculate rainfall intensity [mm/h]
            rintnsty = eros.precipitation[i] / (ts/3600)
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
        qr_land = eros.q_land[i] * 60 / (( eros.xl[i] + eros.yl[i] ) / 2)
        # Sine of the slope
        sinslope = sin(atan(eros.slope[i]))

        # Overland flow erosion [kg/min]
        # For a wide range of slope, it is better to use the sine of slope rather than tangeant
        sedov = eros.erosov[i] * eros.usleC[i] * eros.usleK[i] * eros.xl[i] * eros.yl[i] * sinslope * qr_land
        # [ton/timestep]
        sedov = sedov * (ts/60) * 10^(-3)
        # Remove the impervious area
        sedov = sedov * (1.0 - eros.pathfrac[i])

        # Total soil loss [ton/cell/timestep]
        soilloss = sedspl + sedov

        # Eroded amount per particle class
        erosclay = soilloss * eros.fclay[i]
        erossilt = soilloss * eros.fsilt[i]
        erossand = soilloss * eros.fsand[i]
        erossagg = soilloss * eros.fsagg[i]
        eroslagg = soilloss * eros.flagg[i]
        

        # update the outputs and states
        eros.sedspl[i] = sedspl
        eros.sedov[i] = sedov
        eros.soilloss[i] = soilloss
        eros.erosclay[i] = erosclay
        eros.erossilt[i] = erossilt
        eros.erossand[i] = erossand
        eros.erossagg[i] = erossagg
        eros.eroslagg[i] = eroslagg
    end
end

### Sediment transport capacity in overland flow ###
function update_until_oltransport(ols::LandSed, config)

    # # start dummy variables (should be generated from model reader and from Config.jl TOML)
    do_river = get(config.model, "runrivermodel", false)
    tcmethod = get(config.model, "landtransportmethod", "yalinpart")
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
                delta = max((ols.h_land[i] * sinslope / (ols.D50[i] * 0.001 * (ols.rhos[i]/1000 -1))/0.06 - 1), 0.0)
                alphay = delta*2.45/(0.001*ols.rhos[i])^0.4 * 0.06^(0.5)
                if ols.q_land[i] > 0.0 && alphay != 0.0
                    TC = (ols.dl[i] / ols.q_land[i] * (ols.rhos[i] - 1000) * ols.D50[i] * 0.001 
                    * (9.81 * ols.h_land[i] * sinslope) * 0.635 * delta
                    * (1 - log(1+alphay)/(alphay))) # [kg/m3]
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
            delta = ols.h_land[i] * sinslope / (10^(-6) * (ols.rhos[i]/1000 - 1) / 0.06)
            dclay = max(1 / ols.dmclay[i] * delta - 1, 0.0)
            dsilt = max(1 / ols.dmsilt[i] * delta - 1, 0.0)
            dsand = max(1 / ols.dmsand[i] * delta - 1, 0.0)
            dsagg = max(1 / ols.dmsagg[i] * delta - 1, 0.0)
            dlagg = max(1 / ols.dmlagg[i] * delta - 1, 0.0)
            # Total transportability
            dtot = dclay + dsilt + dsand + dsagg + dlagg

            # Yalin transport capacity of overland flow for each particle class
            if ols.q_land[i] > 0.0
                TCa = ols.dl[i] / ols.q_land[i] * (ols.rhos[i] - 1000) * 10^(-6) * (9.81 * ols.h_land[i] * sinslope)
            else
                TCa = 0.0
            end
            TCb = 2.45 / (ols.rhos[i] / 1000)^0.4 * 0.06^0.5
            if dtot != 0.0 && dclay != 0.0
                TCclay = TCa * ols.dmclay[i] * dclay / dtot * 0.635 * dclay * (1 - log(1 + dclay * TCb) / dclay * TCb) # [kg/m3]
                TCclay = TCclay * ols.q_land[i] *  ts * 10^(-3) # [ton]
            else
                TCclay = 0.0
            end
            if dtot != 0.0 && dsilt != 0.0
                TCsilt = TCa * ols.dmsilt[i] * dsilt / dtot * 0.635 * dsilt * (1 - log(1 + dsilt * TCb) / dsilt * TCb) # [kg/m3]
                TCsilt = TCsilt * ols.q_land[i] *  ts * 10^(-3) # [ton]
            else
                TCsilt = 0.0
            end
            if dtot != 0.0 && dsand != 0.0
                TCsand = TCa * ols.dmsand[i] * dsand / dtot * 0.635 * dsand* (1 - log(1 + dsand * TCb) / dsand * TCb) # [kg/m3]
                TCsand = TCsand * ols.q_land[i] *  ts * 10^(-3) # [ton]
            else
                TCsand = 0.0
            end
            if dtot != 0.0 && dsagg != 0.0
                TCsagg = TCa * ols.dmsagg[i] * dsagg / dtot * 0.635 * dsagg * (1 - log(1 + dsagg * TCb) / dsagg * TCb) # [kg/m3]
                TCsagg = TCsagg * ols.q_land[i] *  ts * 10^(-3) # [ton]
            else
                TCsagg = 0.0
            end
            if dtot != 0.0 && dlagg != 0.0
                TClagg = TCa * ols.dmlagg[i] * dlagg / dtot * 0.635 * dlagg * (1 - log(1 + dlagg * TCb) / dlagg * TCb) # [kg/m3]
                TClagg = TClagg * ols.q_land[i] *  ts * 10^(-3) # [ton]
            else
                TClagg = 0.0
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
Base.@kwdef struct OLFSed{T}
    # number of cells
    n::Int
    # Filter with river cells
    rivcell::Vector{T}
    # Forcing for the river
    h_riv::Vector{T} = fill(mv, n)
    q_riv::Vector{T} = fill(mv, n)
    # Total eroded soil [ton]
    soilloss::Vector{T} = fill(mv, n)
    # Eroded soil per particle class [ton]
    erosclay::Vector{T} = fill(mv, n)
    erossilt::Vector{T} = fill(mv, n)
    erossand::Vector{T} = fill(mv, n)
    erossagg::Vector{T} = fill(mv, n)
    eroslagg::Vector{T} = fill(mv, n)
    # Total transport capacity of overland flow [ton]          
    TCsed::Vector{T} = fill(mv, n)
    # Transport capacity of overland flow per particle class [ton]          
    TCclay::Vector{T} = fill(mv, n)
    TCsilt::Vector{T} = fill(mv, n)
    TCsand::Vector{T} = fill(mv, n)
    TCsagg::Vector{T} = fill(mv, n)
    TClagg::Vector{T} = fill(mv, n)
    # Sediment flux in overland flow [ton]
    olsed::Vector{T} = fill(mv, n)
    olclay::Vector{T} = fill(mv, n)
    olsilt::Vector{T} = fill(mv, n)
    olsand::Vector{T} = fill(mv, n)
    olsagg::Vector{T} = fill(mv, n)
    ollagg::Vector{T} = fill(mv, n)
    # Sediment reaching the river with overland flow [ton]
    inlandsed::Vector{T} = fill(mv, n)
    inlandclay::Vector{T} = fill(mv, n)
    inlandsilt::Vector{T} = fill(mv, n)
    inlandsand::Vector{T} = fill(mv, n)
    inlandsagg::Vector{T} = fill(mv, n)
    inlandlagg::Vector{T} = fill(mv, n)


    function OLFSed{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

statevars(::OLFSed) = ()

function update(ols::OLFSed, network, config)
    do_river = get(config.model, "runrivermodel", false)
    tcmethod = get(config.model, "landtransportmethod", "yalinpart")
    zeroarr = fill(0.0, ols.n)

    if do_river || tcmethod == "yalinpart"
        ols.olclay .= accucapacityflux(network, ols.erosclay, ols.TCclay)
        depclay = accucapacitystate(network, ols.erosclay, ols.TCclay)
        ols.inlandclay .= ifelse.(ols.rivcell .== 1, depclay, zeroarr)
        ols.olsilt .= accucapacityflux(network, ols.erossilt, ols.TCsilt)
        depsilt = accucapacitystate(network, ols.erossilt, ols.TCsilt)
        ols.inlandsilt .= ifelse.(ols.rivcell .== 1, depsilt, zeroarr)
        ols.olsand .= accucapacityflux(network, ols.erossand, ols.TCsand)
        depsand = accucapacitystate(network, ols.erossand, ols.TCsand)
        ols.inlandsand .= ifelse.(ols.rivcell .== 1, depsand, zeroarr)
        ols.olsagg .= accucapacityflux(network, ols.erossagg, ols.TCsagg)
        depsagg = accucapacitystate(network, ols.erossagg, ols.TCsagg)
        ols.inlandsagg .= ifelse.(ols.rivcell .== 1, depsagg, zeroarr)
        ols.ollagg .= accucapacityflux(network, ols.eroslagg, ols.TClagg)
        deplagg = accucapacitystate(network, ols.eroslagg, ols.TClagg)
        ols.inlandlagg .= ifelse.(ols.rivcell .== 1, deplagg, zeroarr)

        ols.olsed .= ols.olclay .+ ols.olsilt .+ ols.olsand .+ ols.olsagg .+ ols.ollagg
        ols.inlandsed .= ols.inlandclay .+ ols.inlandsilt .+ ols.inlandsand .+ ols.inlandsagg .+ ols.inlandlagg
    else
        ols.olsed .= accucapacityflux(network, ols.soilloss, ols.TCsed)
        depsed = accucapacitystate(network, ols.soilloss, ols.TCsed)
        ols.inlandsed .= ifelse.(ols.rivcell .== 1, depsed, zeroarr)
    end
end

### River transport and processes ###
Base.@kwdef struct RiverSed{T}
    # number of cells
    n::Int
    # Timestep [s]
    Δt::T
    # River geometry (slope [-], length [m], width [m])
    sl::Vector{T}
    dl::Vector{T}
    width::Vector{T}
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
    fclayriv::Vector{T}
    fsiltriv::Vector{T}
    fsandriv::Vector{T}
    fgravriv::Vector{T}
    # Sediment mean diameter for Engelund and Hansen transport equation [mm]
    d50engelund::Vector{T}
    # Parameters for Bagnold transport equation
    cbagnold::Vector{T}
    ebagnold::Vector{T}
    # Parameters for Kodatie transport equation
    ak::Vector{T}
    bk::Vector{T}
    ck::Vector{T}
    dk::Vector{T}
    # Critical bed and bank shear stress [N/m2]
    TCrbank::Vector{T}
    TCrbed::Vector{T}
    # Bed and bank erodibilities [m3/N.s]
    kdbank::Vector{T}
    kdbed::Vector{T}
    # Sediment density [kg/m3]
    rhos::Vector{T}
    # River water level [m]
    h_riv::Vector{T} = fill(mv, n)
    # River discharge [m3/s]
    q_riv::Vector{T} = fill(mv, n)
    # Sediment input from land erosion [ton]
    inlandclay::Vector{T} = fill(0.0, n)
    inlandsilt::Vector{T} = fill(0.0, n)
    inlandsand::Vector{T} = fill(0.0, n)
    inlandsagg::Vector{T} = fill(0.0, n)
    inlandlagg::Vector{T} = fill(0.0, n)
    # Sediment / particle left in the cell [ton]
    sedload::Vector{T} = fill(0.0, n)
    clayload::Vector{T} = fill(0.0, n)
    siltload::Vector{T} = fill(0.0, n)
    sandload::Vector{T} = fill(0.0, n)
    saggload::Vector{T} = fill(0.0, n)
    laggload::Vector{T} = fill(0.0, n)
    gravload::Vector{T} = fill(0.0, n)
    # Sediment / particle stored on the river bed after deposition [ton]
    sedstore::Vector{T} = fill(0.0, n)
    claystore::Vector{T} = fill(0.0, n)
    siltstore::Vector{T} = fill(0.0, n)
    sandstore::Vector{T} = fill(0.0, n)
    saggstore::Vector{T} = fill(0.0, n)
    laggstore::Vector{T} = fill(0.0, n)
    gravstore::Vector{T} = fill(0.0, n)
    # Sediment / particle flux [ton]
    outsed::Vector{T} = fill(0.0, n)
    outclay::Vector{T} = fill(0.0, n)
    outsilt::Vector{T} = fill(0.0, n)
    outsand::Vector{T} = fill(0.0, n)
    outsagg::Vector{T} = fill(0.0, n)
    outlagg::Vector{T} = fill(0.0, n)
    outgrav::Vector{T} = fill(0.0, n)
    # Sediment concentrations [mg/L]
    Sedconc::Vector{T} = fill(0.0, n)
    SSconc::Vector{T} = fill(0.0, n)
    Bedconc::Vector{T} = fill(0.0, n)
    

    function RiverSed{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

statevars(::RiverSed) = (
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

function update(rs::RiverSed, network, config)
    @unpack graph, order = network
    tcmethod = get(config.model, "rivtransportmethod", "bagnold")

    # River sediment loads are separated into different particle class.
    # Clay, silt and sand can both come from land, resuspension or river channel erosion.
    # Small and large aggregates only come from land erosion or resuspension.
    # Gravel only comes from resuspension or river channel erosion.

    for v in order
        ### Sediment input in the cell (left from previous timestep + from land + from upstream outflux) ###
        upstream_nodes = inneighbors(graph, v)
        if !isempty(upstream_nodes)
            inrivclay = sum(rs.outclay[i] for i in upstream_nodes)
            inrivsilt = sum(rs.outsilt[i] for i in upstream_nodes)
            inrivsand = sum(rs.outsand[i] for i in upstream_nodes)
            inrivsagg = sum(rs.outsagg[i] for i in upstream_nodes)
            inrivlagg = sum(rs.outlagg[i] for i in upstream_nodes)
            inrivgrav = sum(rs.outgrav[i] for i in upstream_nodes)
        else
            inrivclay = 0.0
            inrivsilt = 0.0
            inrivsand = 0.0
            inrivsagg = 0.0
            inrivlagg = 0.0
            inrivgrav = 0.0
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
            vmean = ifelse(rs.h_riv[v] > 0.0, rs.q_riv[v] / (rs.width[v] * rs.h_riv[v]), 0.0)
            vshear = (9.81 * hydrad * rs.sl[v])^0.5
            # Concentration by weight
            cw  = ifelse(
                hydrad > 0.0,
                (rs.rhos[v]/1000 * 0.5 * vmean * vshear^3 /
                ((rs.rhos[v]/1000 - 1)^2 * 9.81^2 * rs.d50engelund[v] * hydrad)),
                0.0
            )
            cw = min(1.0, cw)
            # Transport capacity [tons/m3]
            maxsed = max(cw / (cw + (1 - cw) * rs.rhos[v]/1000) * rs.rhos[v]/1000, 0.0)
        elseif tcmethod == "bagnold"
            maxsed = rs.cbagnold[v] * (rs.q_riv[v] / (rs.h_riv[v] * rs.width[v]))^rs.ebagnold[v]
        elseif tcmethod == "kodatie"
            vmean = ifelse(rs.h_riv[v] > 0.0, rs.q_riv[v] / (rs.width[v] * rs.h_riv[v]), 0.0)
            maxsed = rs.ak[v] * vmean^rs.bk[v] * rs.h_riv[v]^rs.ck[v] * rs.sl[v]^rs.dk[v]
            # Transport capacity [tons/m3]
            maxsed = ifelse(
                rs.q_riv[v] > 0.0, 
                maxsed * rs.width[v] / (rs.q_riv[v] * rs.Δt),
                0.0
            )
        elseif tcmethod == "yang"
            ws = 411 * rs.d50[v]^2 / 3600
            vshear = (9.81 * hydrad * rs.sl[v])^0.5
            var1 = vshear * rs.d50[v] / 1000 / (1.16 * 10^(-6))
            var2 = ws * rs.d50[v] / 1000 / (1.16 * 10^(-6))
            vcr = min(0.0, ifelse(var1 >= 70.0, 2.05 * ws, ws*(2.5 / (log10(var1) - 0.06) + 0.66)))
            # Sand equation
            if (rs.width[v] * rs.h_riv[v]) >= vcr && rs.d50[v] < 2.0
                logcppm = (
                    5.435 - 0.286 * log10(var2) - 0.457 * log10(vshear / ws)
                    + 1.799 - 0.409 * log10(var2) - 0.314 * log10(vshear / ws)
                    * log10((rs.q_riv[v] / (rs.width[v] * rs.h_riv[v]) - vcr) * rs.sl[v] / ws)
                )
            # Gravel equation
            elseif (rs.width[v] * rs.h_riv[v]) >= vcr && rs.d50[v] < 2.0
                logcppm = (
                    6.681 - 0.633 * log10(var2) - 4.816 * log10(vshear / ws)
                    + 2.784 - 0.305 * log10(var2) - 0.282 * log10(vshear / ws)
                    * log10((rs.q_riv[v] / (rs.width[v] * rs.h_riv[v]) - vcr) * rs.sl[v] / ws)
                )
            else
                logcppm = 0.0
            end
            # Sediment concentration by weight
            cw = 10^(logcppm) * 10^(-6)
            # Transport capacity [ton/m3]
            maxsed = max(cw / (cw + (1 - cw) * rs.rhos[v]/1000) * rs.rhos[v]/1000, 0.0)
        elseif tcmethod == "molinas"
            ws = 411 * rs.d50[v]^2 / 3600
            vmean = ifelse(rs.h_riv[v] > 0.0, rs.q_riv[v] / (rs.width[v] * rs.h_riv[v]), 0.0)
            if rs.h_riv[v] > 0.0
                psi = (
                    vmean^3 / (
                        (rs.rhos[v] / 1000 - 1)
                        * 9.81 * rs.h_riv[v] * ws
                        * log10(1000 * rs.h_riv[v] / rs.d50[v])^2
                    )
                )
            else
                psi = 0.0
            end
            # Concentration by weight
            cw = 1430 * (0.86 + psi^0.5) * psi^1.5 / (0.016 + psi) * 10^(-6)
            # Transport capacity [ton/m3]
            maxsed = max(cw / (cw + (1 - cw) * rs.rhos[v]/1000) * rs.rhos[v]/1000, 0.0)
        end
        
        # 1285 g/L: boundary between streamflow and debris flow (Costa, 1988)
        maxsed = min(maxsed, 1.285)
        # Transport capacity [ton]
        maxsed = maxsed * (rs.h_riv[v] * rs.width[v] * rs.dl[v] + rs.q_riv[v] * rs.Δt)

        ### River erosion ###
        # Erosion only if the load is below the transport capacity of the flow.
        sedex = max(maxsed - insed, 0.0)
        # Bed and bank are eroded only if the previously deposited material is not enough
        rs.sedstore[v] = rs.claystore[v] + rs.siltstore[v] + rs.sandstore[v] + rs.saggstore[v] + rs.laggstore[v] + rs.gravstore[v]
        if sedex > 0.0 && sedex > rs.sedstore[v]
            # Effective sediment needed fom river bed and bank erosion [ton]
            effsedex = sedex - rs.sedstore[v]

            # Repartition of the effective shear stress between the bank and the bed from Knight et al. 1984 [%]
            SFbank = ifelse(
                rs.h_riv[v] > 0.0,
                exp(-3.23 * log10(rs.width[v] / rs.h_riv[v] + 3) + 6.146),
                0.0 
            )
            # Effective shear stress on river bed and banks [N/m2]
            TEffbank = ifelse(
                rs.h_riv[v] > 0.0,
                1000 * 9.81 * hydrad * rs.sl[v] * SFbank / 100 * (1 + rs.width[v] / (2 * rs.h_riv[v])),
                0.0
            )
            TEffbed = 1000 * 9.81 * hydrad * rs.sl[v] * (1 - SFbank / 100) * (1 + 2 * rs.h_riv[v] / rs.width[v])
            # Potential erosion rates of the bed and bank [t/cell/timestep] 
            #(assuming only one bank is eroding)
            Tex = max(TEffbank - rs.TCrbank[v], 0.0)
            # 1.4 is bank default bulk density
            ERbank = max(0.0, rs.kdbank[v] * Tex * rs.dl[v] * rs.h_riv[v] * 1.4 * rs.Δt)
            # 1.5 is bed default bulk density
            ERbed = max(0.0, rs.kdbed[v] * (TEffbed - rs.TCrbed[v]) * rs.dl[v] * rs.width[v] * 1.5 * rs.Δt)
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
            degstoresed = degstoreclay + degstoresilt + degstoresagg + degstoresand + degstorelagg + degstoregrav
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
            ifelse(rs.dmclay[v] <= dbedf, rs.outclay[v] / 2, 0.0)
        )
        SSsilt = ifelse(
            rs.dmsilt[v] <= dsuspf,
            rs.outsilt[v],
            ifelse(rs.dmsilt[v] <= dbedf, rs.outsilt[v] / 2, 0.0)
        )
        SSsagg = ifelse(
            rs.dmsagg[v] <= dsuspf,
            rs.outsagg[v],
            ifelse(rs.dmsagg[v] <= dbedf, rs.outsagg[v] / 2, 0.0)
        )
        SSsand = ifelse(
            rs.dmsand[v] <= dsuspf,
            rs.outsand[v],
            ifelse(rs.dmsand[v] <= dbedf, rs.outsand[v] / 2, 0.0)
        )
        SSlagg = ifelse(
            rs.dmlagg[v] <= dsuspf,
            rs.outlagg[v],
            ifelse(rs.dmlagg[v] <= dbedf, rs.outlagg[v] / 2, 0.0)
        )
        SSgrav = ifelse(
            rs.dmgrav[v] <= dsuspf,
            rs.outgrav[v],
            ifelse(rs.dmgrav[v] <= dbedf, rs.outgrav[v] / 2, 0.0)
        )

        SS = SSclay + SSsilt + SSsagg + SSsand + SSlagg + SSgrav
        Bed = rs.outsed[v] - SS

        rs.SSconc[v] = SS * toconc
        rs.Bedconc[v] = Bed * toconc

    end


end