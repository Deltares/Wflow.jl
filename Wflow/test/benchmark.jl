# These benchmarks need to be revised, and run not as part of the tests, but separately.
# https://juliaci.github.io/PkgBenchmark.jl/stable/

#=

"Prints a benchmark results just like btime"
function print_benchmark(trialmin)
    trialtime = BenchmarkTools.time(trialmin)
    trialallocs = BenchmarkTools.allocs(trialmin)
    println(
        "  ",
        BenchmarkTools.prettytime(trialtime),
        " (",
        trialallocs,
        " allocation",
        trialallocs == 1 ? "" : "s",
        ": ",
        BenchmarkTools.prettymemory(BenchmarkTools.memory(trialmin)),
        ")",
    )
end

# test/run_hbv.jl
benchmark = @benchmark Wflow.update(model)
trialmin = BenchmarkTools.minimum(benchmark)
println("HBV Model update")
print_benchmark(trialmin)

# test/run_sbm_gwf.jl
benchmark = @benchmark Wflow.run(tomlpath)
trialmin = BenchmarkTools.minimum(benchmark)
println("SBM GWF Model update (run)")
print_benchmark(trialmin)

# test/run_sbm.jl
benchmark = @benchmark Wflow.update(model)
trialmin = BenchmarkTools.minimum(benchmark)
println("SBM Model update")
print_benchmark(trialmin)

# test/run_sediment.jl
benchmark = @benchmark Wflow.run(tomlpath)
trialmin = BenchmarkTools.minimum(benchmark)
println("Sediment Model update (run)")
print_benchmark(trialmin)
=#
