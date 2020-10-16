Base.@kwdef struct EROS{T}
    # number of cells
    n::Int
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
    ## Interception related to leaf_area_index climatology ###
    # Specific leaf storage [mm]
    sl::Vector{T} = fill(mv, n)
    # Storage woody part of vegetation [mm]
    swood::Vector{T} = fill(mv, n)
    # Extinction coefficient [-] (to calculate canopy gap fraction)
    kext::Vector{T} = fill(mv, n)
    # Leaf area index [m² m⁻²]
    leaf_area_index::Vector{T} = fill(mv, n)

    function EROS{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

statevars(::EROS) = ()

function update_until_ols(eros::EROS, config)

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
        

        # update the outputs and states
        eros.canopygapfraction[i] = canopygapfraction
        eros.sedspl[i] = sedspl
        eros.sedov[i] = sedov
        eros.soilloss[i] = soilloss
    end
end
