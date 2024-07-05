"""
    parse_loglevel(input_level::AbstractString)::LogLevel
    parse_loglevel(input_level::Integer)::LogLevel

Parse a log level from either an integer or string.

# Examples
    parse_loglevel("info") -> Logging.Info
    parse_loglevel(0) -> LogLevel(0) (== Logging.Info)
"""
function parse_loglevel(input_level::AbstractString)::LogLevel
    level = lowercase(input_level)
    if level == "debug"
        return Logging.Debug
    elseif level == "info"
        return Logging.Info
    elseif level == "warn"
        return Logging.Warn
    elseif level == "error"
        return Logging.Error
    else
        error("loglevel $input_level not recognized")
    end
end

parse_loglevel(input_level::Integer) = LogLevel(input_level)

"Print a log message to a single line, for Delft-FEWS"
function format_message(io::IO, args)::Nothing
    kwargs = IOBuffer()
    for p in args.kwargs
        print(kwargs, string(" |", p[1], " = ", p[2], "|"))
    end
    println(
        io,
        now(),
        " | ",
        args._module,
        " | [",
        args.level,
        "] ",
        args.message,
        String(take!(kwargs)),
    )
    return nothing
end

"Initialize a logger, which is different if `fews_run` is set in the Config."
function init_logger(config::Config; silent = false)::Tuple{TeeLogger, IOStream}
    loglevel = parse_loglevel(get(config, "loglevel", "info"))
    path_log = output_path(config, get(config, "path_log", "log.txt"))
    mkpath(dirname(path_log))
    log_handle = open(path_log, "w")
    fews_run = get(config, "fews_run", false)::Bool

    file_logger = if fews_run
        # Format the log message to be printed on a single line.
        MinLevelLogger(FormatLogger(format_message, log_handle), loglevel)
    else
        # use ConsoleLogger instead of FileLogger to be able to print full stacktraces there
        # https://github.com/JuliaLogging/LoggingExtras.jl/issues/46#issuecomment-803716480
        ConsoleLogger(
            IOContext(log_handle, :compact => false, :limit => false, :color => false),
            loglevel;
            show_limited = false,
        )
    end

    terminal_logger = if silent
        NullLogger()
    else
        # avoid debug information overflowing the terminal, -1 is the level of ProgressLogging
        # avoid double stacktraces by filtering the @error catch in Wflow.run
        EarlyFilteredLogger(
            log -> log.id !== :wflow_run,
            MinLevelLogger(TerminalLogger(), LogLevel(-1)),
        )
    end

    logger = TeeLogger(
        # avoid progresslogger and NCDatasets debug messages in the file
        EarlyFilteredLogger(
            log -> log.group !== :ProgressLogging && log._module != NCDatasets,
            file_logger,
        ),
        terminal_logger,
    )
    with_logger(logger) do
        @debug "Logger initialized." loglevel silent fews_run path_log
    end
    return logger, log_handle
end
