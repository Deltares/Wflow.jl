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
    (; loglevel) = config.logging
    path_log = output_path(config, config.logging.path_log)
    mkpath(dirname(path_log))
    log_handle = open(path_log, "w")

    fews_run = config.fews_run__flag
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
