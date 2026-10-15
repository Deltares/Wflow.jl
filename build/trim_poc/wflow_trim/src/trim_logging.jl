# Logging for trimmed binaries.
#
# JuliaC disables the logging macros under `--trim` by overriding
# `current_logger_for_env(std_level, group, _module) = nothing`. Here we add more specific
# methods for a single concrete logger type, such that `@info` and friends dispatch
# statically to it. This is type piracy on Base, which is only acceptable since this module
# is the entry point of a trimmed binary and never loaded in a normal Julia session.

using Base.CoreLogging: CoreLogging, LogLevel, AbstractLogger, Debug, Info, Warn, Error
using Wflow: Dates

"Logger that writes to the terminal and to an optional log file."
struct TrimLogger <: AbstractLogger
    min_level::LogLevel
    silent::Bool
end

const LOGGER = Ref(TrimLogger(Info, false))
const LOG_FILE = Ref{IOStream}()

"Configure the logger from the TOML, the log file path is relative to `dir_output`."
function init_trim_logger(config::Wflow.Config)::Nothing
    (; loglevel, silent, path_log) = config.logging
    path = Wflow.output_path(config, path_log)
    mkpath(dirname(path))
    LOG_FILE[] = open(path, "w")
    LOGGER[] = TrimLogger(loglevel, silent)
    return nothing
end

function close_trim_logger()::Nothing
    isassigned(LOG_FILE) && close(LOG_FILE[])
    return nothing
end

CoreLogging.min_enabled_level(logger::TrimLogger) = logger.min_level
CoreLogging.shouldlog(::TrimLogger, level, _module, group, id) = true
CoreLogging.catch_exceptions(::TrimLogger) = true

function level_name(level::LogLevel)::String
    return if level == Debug
        "Debug"
    elseif level == Info
        "Info"
    elseif level == Warn
        "Warning"
    elseif level == Error
        "Error"
    else
        string("LogLevel(", level.level, ")")
    end
end

"""
Print values of commonly logged types; other values only by their type name. Like the
`ConsoleLogger`, strings are quoted for keyword arguments (`quote_strings = true`).
"""
function print_value(io::IO, @nospecialize(value), quote_strings::Bool = false)::Nothing
    if value isa String
        quote_strings ? show(io, value::String) : print(io, value::String)
    elseif value isa SubString{String}
        quote_strings ? show(io, value::SubString{String}) : print(io, value::SubString{String})
    elseif value isa Symbol
        print(io, value::Symbol)
    elseif value isa Bool
        print(io, value::Bool)
    elseif value isa Int
        print(io, value::Int)
    elseif value isa Float64
        print(io, value::Float64)
    elseif value isa Nothing
        print(io, "nothing")
    elseif value isa Dates.DateTime
        print(io, Wflow.format_datetime(value::Dates.DateTime))
    elseif value isa Dates.Second
        print(io, Dates.value(value::Dates.Second), " seconds")
    elseif value isa Exception
        print(io, error_message(value))
    elseif value isa Tuple{Exception, Any}
        # `exception = (e, catch_backtrace())`, the backtrace cannot be printed
        print(io, error_message(first(value::Tuple{Exception, Any})))
    else
        print(io, "<", nameof(typeof(value)), ">")
    end
    return nothing
end

"Write a log record in the `ConsoleLogger` layout."
function write_record(io::IO, level::LogLevel, @nospecialize(message), kwargs)::Nothing
    n_kwargs = length(kwargs)
    print(io, iszero(n_kwargs) ? "[ " : "┌ ", level_name(level), ": ")
    print_value(io, message)
    println(io)
    for (kwarg_idx, (key, value)) in enumerate(pairs(kwargs))
        print(io, kwarg_idx == n_kwargs ? "└   " : "│   ", key, " = ")
        print_value(io, value, true)
        println(io)
    end
    return nothing
end

function CoreLogging.handle_message(
        logger::TrimLogger, level::LogLevel, message, _module, group, id, file, line;
        kwargs...,
    )
    # Wflow.run logs the error that it rethrows, avoid printing it twice to the terminal
    if !logger.silent && id !== :wflow_run
        write_record(Core.stdout, level, message, kwargs)
    end
    if isassigned(LOG_FILE)
        write_record(LOG_FILE[], level, message, kwargs)
    end
    return nothing
end

"""
Only log messages from these modules. For other modules the logger is `nothing` at compile
time, such that their log statements stay dead code, which is required if the message
construction is not trim compatible, like in `ThreadingUtilities.__init__`.
"""
@inline is_logged_module(_module::Module) =
    _module === Wflow || _module === Wflow.NC || _module === @__MODULE__

# More specific than JuliaC's `current_logger_for_env(::LogLevel, group, _module)`. The
# logging macros pass a constant `Symbol` group and `Module`.
@inline function CoreLogging.current_logger_for_env(
        std_level::LogLevel, group::Symbol, _module::Module,
    )
    is_logged_module(_module) || return nothing
    logger = LOGGER[]
    return std_level >= logger.min_level ? logger : nothing
end

# The macros call `handle_message_nothrow`, which is `@nospecialize` and uses
# `@invokelatest`. JuliaC turns `invokelatest` into a direct call, and the methods below
# for the concrete logger make the dispatch static.
@noinline function CoreLogging.handle_message_nothrow(
        logger::TrimLogger, level, msg, _module, group, id, file, line; kwargs...,
    )
    try
        CoreLogging.handle_message(
            logger, level, msg, _module, group, id, file, line; kwargs...,
        )
    catch
        println(Core.stderr, "Error while writing a log message")
    end
    return nothing
end

# Called when evaluating the log message or its keyword arguments throws
function CoreLogging.logging_error(
        logger::TrimLogger, level, _module, group, id, filepath, line, err, real::Bool,
    )
    println(Core.stderr, "Exception while generating log record: ", error_message(err))
    return nothing
end
