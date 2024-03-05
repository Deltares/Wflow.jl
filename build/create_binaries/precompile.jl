# this script is run during app creation to precompile the code executed in this script that
# way the app will have less latency. This script assumes the test data (download_test_data.jl)
# has been run before.

using Wflow

testdir = abspath(dirname(pathof(Wflow)), "..", "test")
Wflow.run(joinpath(testdir, "sbm_config.toml"))
Wflow.run(joinpath(testdir, "sbm_gwf_config.toml"))
Wflow.run(joinpath(testdir, "hbv_config.toml"))
Wflow.run(joinpath(testdir, "sediment_config.toml"))
