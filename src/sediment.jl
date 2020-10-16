Base.@kwdef struct EROS{T}
    # number of cells
    n::Int
    # length of cells in y direction [m]
    yl::Vector{T}
    # length of cells in x direction [m]
    xl::Vector{T}
    # Fraction of river [-]
    riverfrac::Vector{T}
    # Precipitation [mm]
    precipitation::Vector{T} = fill(mv, n)
    # Overland flow [m3/s]
    q_land::Vector{T} = fill(mv, n)
    # Coefficient for ANSWERS overland flow erosion [-]
    erosov::Vector{T}
    # Fraction of impervious area per grid cell [-]
    pathfrac::Vector{T}
    #
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

    function EROS{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

statevars(::EROS) = ()

function update_until_ols(eros::EROS, config)

    # # start dummy variables (should be generated from model reader and from Config.jl TOML)
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
        eros.sedspl[i] = sedspl
        eros.sedov[i] = sedov
        eros.soilloss[i] = soilloss
    end
end
