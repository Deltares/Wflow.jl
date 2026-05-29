using TestItemRunner
using LoggingExtras
using Base.Threads

@info "testing Wflow with" nthreads() VERSION

const pattern = length(ARGS) > 0 ? ARGS[1] : ""

with_logger(NullLogger()) do
    @run_package_tests filter = ti -> occursin(pattern, ti.name)
end
