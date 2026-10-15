module wflow_trim

using Wflow
using Wflow: ModelType, RoutingType, VerticalConductivityProfile

"""
Build profile: the model options that are supported by this binary. These are compile time
constants, only code for these options is compiled. A TOML with other options gives an
error. The profile can be overridden at build time with the `WFLOW_TRIM_PROFILE`
environment variable, containing a NamedTuple expression.
"""
const DEFAULT_PROFILE = (
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

const PROFILE = if haskey(ENV, "WFLOW_TRIM_PROFILE")
    merge(DEFAULT_PROFILE, eval(Meta.parse(ENV["WFLOW_TRIM_PROFILE"])))
else
    DEFAULT_PROFILE
end

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

include("trim_logging.jl")

function (@main)(args::Vector{String})::Cint
    if length(args) != 1
        println(Core.stdout, "Usage: wflow_trim path/to/config.toml")
        return 1
    end
    toml_path = only(args)
    # The Wflow logger setup (TerminalLoggers, LoggingExtras) is not trim compatible, it is
    # replaced by the `TrimLogger`, see trim_logging.jl.
    start = time()
    try
        config = Wflow.Config(toml_path)
        init_trim_logger(config)
        @info "Wflow version `v$(Wflow.VERSION)` (trimmed build)"
        Wflow.run(config, (Wflow.SbmModel{PROFILE}(),))
    catch e
        message = error_message(e)
        @error "Wflow simulation failed" exception = message _id = :wflow_run
        println(Core.stderr, "Wflow simulation failed: ", message)
        close_trim_logger()
        return 1
    end
    close_trim_logger()
    println(Core.stdout, "Simulation finished in ", round(time() - start; digits = 2), " s")
    return 0
end

end # module
