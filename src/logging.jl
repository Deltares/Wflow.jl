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
        "timestamp = ",
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
function init_logger(config::Config)::Tuple{TeeLogger,IOStream}
    loglevel = parse_loglevel(get(config, "loglevel", "info"))
    path_log = output_path(config, get(config, "path_log", "log.txt"))
    log_handle = open(path_log, "w")
    fews_run = get(config, "fews_run", false)::Bool

    file_logger = if fews_run
        # Format the log message to be printed on a single line.
        FormatLogger(format_message, log_handle)
    else
        f_logger = FileLogger(log_handle)
        # Add a timestamp as a keyword argument to a log message.
        # Added as an internal option for debugging performance issues, which is
        # normally off since it creates a lot of noise in the log file and makes it
        # non reproducable.
        add_timestamps = false
        if add_timestamps
            f_logger = TransformerLogger(f_logger) do orig_log
                new_kwargs = (orig_log.kwargs..., timestamp = now())
                merge(orig_log, (kwargs = new_kwargs,))
            end
        end
        f_logger
    end

    logger = TeeLogger(
        # avoid progresslogger and NCDatasets debug messages in the file
        EarlyFilteredLogger(
            log -> log.group !== :ProgressLogging && log._module != NCDatasets,
            MinLevelLogger(file_logger, loglevel),
        ),
        # avoid debug information overflowing the terminal, -1 is the level of ProgressLogging
        MinLevelLogger(TerminalLogger(), LogLevel(-1)),
    )
    # Delft-FEWS already gets a timestamp from format_message, otherwise add it as a kwarg.
    return logger, log_handle
end
