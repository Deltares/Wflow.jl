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
    TCsed::Vector{T} = fill(0.0, n)
    # Transport capacity of overland flow per particle class [ton]          
    TCclay::Vector{T} = fill(0.0, n)
    TCsilt::Vector{T} = fill(0.0, n)
    TCsand::Vector{T} = fill(0.0, n)
    TCsagg::Vector{T} = fill(0.0, n)
    TClagg::Vector{T} = fill(0.0, n)

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
        eros.canopygapfraction[i] = canopygapfraction
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