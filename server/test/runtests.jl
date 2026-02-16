using TestItemRunner
using Logging: with_logger, NullLogger

with_logger(NullLogger()) do
    @run_package_tests
end
