using TerminalLoggers
using LoggingExtras
using ProgressLogging
using Dates
using TOML
using Test
using Wflow

tomlpath = raw"c:\Users\visser_mn\.julia\dev\Wflow\test\sbm_simple.toml"

# ProgressLogging shows in VSCode
# @info shows in terminal
Wflow.run(tomlpath);

# ProgressLogging shows in terminal
# @info shows in terminal
with_logger(TerminalLogger()) do
    Wflow.run(tomlpath)
end

config = Wflow.Config(tomlpath)

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

@test parse_loglevel("InfO") == Logging.Info
@test parse_loglevel(0) == Logging.Info

# minimum loglevel to include in the logfile, Debug by default
loglevel = parse_loglevel(get(config, "loglevel", "debug"))

log_handle = open("log.txt", "w")
logger = TeeLogger(
    # EarlyFilteredLogger(x->x._module != NCDatasets, MinLevelLogger(FileLogger("log-t$n.txt"), Logging.Debug)),
    # avoid progresslogger in the file, instead we put more useful debug statements
    # at the start of every timestep
    EarlyFilteredLogger(log->log.group !== :ProgressLogging,
        MinLevelLogger(
            FileLogger(log_handle, always_flush=true),
            loglevel)
    ),
    # avoid debug information overflowing the terminal, -1 is the level of ProgressLogging
    MinLevelLogger(TerminalLogger(), LogLevel(-1)),
)

# always add add timestamp to FileLogger
# (level, message, _module, group, id, file, line, kwargs)
timestamp_logger(logger) = TransformerLogger(logger) do orig_log
    new_kwargs = (orig_log.kwargs..., timestamp=now())
    new_log = merge(orig_log, (kwargs=new_kwargs,))
    return new_log
end

tslogger = timestamp_logger(logger)
with_logger(tslogger) do
    @debug "Verbose debugging information.  Invisible by default" pi
    @info "An informational message" 1.2
    @warn "Something was odd.  You should pay attention" rand()
    @error "A non fatal error occurred" randn()
    @progress for _ = 1:50; sleep(0.1); end
    log(-1)
    close(log_handle)
end


# TODO check how errors/stacktraces are logged, should go in the log file
# to catch stacktraces, put everything in a try catch, like wflow_cli
# https://github.com/JuliaLogging/LoggingExtras.jl/issues/46

# TODO create single line logging for FEWS
# https://publicwiki.deltares.nl/display/FEWSDOC/05+General+Adapter+Module#id-05GeneralAdapterModule-logFile

with_logger(logger) do
    Wflow.run(tomlpath)
end

Wflow.run(tomlpath)
# Wflow.run(tomlpath, logging=false)

with_logger(TerminalLogger()) do
    @info 1+1
end

@progress for _ = 1:50; sleep(0.1); end


s = open("test1.txt","w")
println(s, now())
flush(s)

lg = TeeLogger(
    FileLogger("test2.txt", append=false, always_flush=true),
    TerminalLogger()
)
with_logger(lg) do
    @info "to file" now()
end
