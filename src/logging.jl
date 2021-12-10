"Parse log level"
function parse_loglevel(input_level::AbstractString)
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

# always add add timestamp to FileLogger
# (level, message, _module, group, id, file, line, kwargs)
timestamp_logger(logger) =
    TransformerLogger(logger) do orig_log
        new_kwargs = (orig_log.kwargs..., timestamp = now())
        new_log = merge(orig_log, (kwargs = new_kwargs,))
        return new_log
    end

"Initialize logger for Julia/wflow_cli"
function init_logger(log_handle, loglevel)
    logger = TeeLogger(
        # avoid progresslogger and NCDatasets debug messages in the file
        EarlyFilteredLogger(
            log -> log.group !== :ProgressLogging && log._module != NCDatasets,
            MinLevelLogger(FileLogger(log_handle, always_flush = true), loglevel),
        ),
        # avoid debug information overflowing the terminal, -1 is the level of ProgressLogging
        MinLevelLogger(TerminalLogger(), LogLevel(-1)),
    )

    tslogger = timestamp_logger(logger)
    return tslogger
end

"Initialize logger for Delft-FEWS run"
function init_logger_fews(log_handle, loglevel)
    # create single line logging for FEWS
    function format_message(io, args)
        println(
            io,
            "timestamp = $(now()) | ",
            args._module,
            " | ",
            "[",
            args.level,
            "] ",
            args.message,
        )
    end

    logger = TeeLogger(
        # avoid progresslogger and NCDatasets debug messages in the file
        EarlyFilteredLogger(
            log -> log.group !== :ProgressLogging && log._module != NCDatasets,
            MinLevelLogger(
                FormatLogger(format_message, log_handle, always_flush = true),
                loglevel,
            ),
        ),
        # avoid debug information overflowing the terminal, -1 is the level of ProgressLogging
        MinLevelLogger(TerminalLogger(), LogLevel(-1)),
    )

    return logger
end

"Initialize logger either running from Delft-FEWS or Julia/wflow_cli"
function init_logger(config::Config)
    loglevel = parse_loglevel(get(config, "loglevel", "info"))
    log_handle = open(output_path(config, "log.txt"), "w")
    fews_run = get(config, "fews_run", false)::Bool

    if fews_run
        return init_logger_fews(log_handle, loglevel), log_handle
    else
        return init_logger(log_handle, loglevel), log_handle
    end
end
