module wflow_trim

using Wflow
using Wflow: ModelType, RoutingType, VerticalConductivityProfile

"""
Build profile: the model options that are supported by this binary. These are compile time
constants, only code for these options is compiled. A TOML with other options gives an
error.
"""
const PROFILE = (
    soil_layers = 3,
    gash_interception = true,
    snow = false,
    glacier = false,
    land_surface_temperature = false,
    water_demand = false,
    reservoir = false,
    floodplain_1d = false,
    drain = false,
    river_subsurface_exchange_head_based = false,
    water_mass_balance = false,
    land_routing = RoutingType.kinematic_wave,
    river_routing = RoutingType.kinematic_wave,
    saturated_hydraulic_conductivity_profile = VerticalConductivityProfile.exponential,
)

"Error message for common exception types, `showerror` needs dynamic dispatch."
function error_message(@nospecialize(e))::String
    msg = if e isa ErrorException
        e.msg
    elseif e isa ArgumentError
        e.msg
    else
        nothing
    end
    return if msg isa String
        msg
    elseif e isa Wflow.NC.NetCDFError
        string("NetCDF error code ", e.code, ": ", e.msg)
    elseif e isa SystemError
        string(e.prefix, ": errno ", e.errnum)
    else
        string("exception of type ", nameof(typeof(e)))
    end
end

function (@main)(args::Vector{String})::Cint
    if length(args) != 1
        println(Core.stdout, "Usage: wflow_trim path/to/config.toml")
        return 1
    end
    toml_path = only(args)
    # Logging macros are disabled in trimmed binaries, and the Wflow logger setup
    # (TerminalLoggers, LoggingExtras) is not trim compatible, so it is bypassed here.
    println(Core.stdout, "Wflow v", Wflow.VERSION, " (trimmed build), running ", toml_path)
    start = time()
    try
        config = Wflow.Config(toml_path)
        Wflow.run(config, (Wflow.SbmModel{PROFILE}(),))
    catch e
        println(Core.stderr, "Wflow simulation failed: ", error_message(e))
        return 1
    end
    println(Core.stdout, "Simulation finished in ", round(time() - start; digits = 2), " s")
    return 0
end

end # module
