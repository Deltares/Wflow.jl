using TestItemRunner
using LoggingExtras
using Base.Threads

@info "testing Wflow with" nthreads() VERSION

with_logger(NullLogger()) do
    @run_package_tests
end
