using Revise
using Wflow
tomlpath = joinpath(@__DIR__, "test", "sediment_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sediment_model(config)
#@unpack network = model

model = Wflow.update(model)

#flush(model.writer.csv_io)

Wflow.close_files(model)

# benchmark = @benchmark Wflow.update(model)

# "Prints a benchmark results just like btime"
# function print_benchmark(trialmin)
#     trialtime = BenchmarkTools.time(trialmin)
#     trialallocs = BenchmarkTools.allocs(trialmin)
#     println(
#         "  ",
#         BenchmarkTools.prettytime(trialtime),
#         " (",
#         trialallocs,
#         " allocation",
#         trialallocs == 1 ? "" : "s",
#         ": ",
#         BenchmarkTools.prettymemory(BenchmarkTools.memory(trialmin)),
#         ")",
#     )
# end

# trialmin = BenchmarkTools.minimum(benchmark)

# println("Sediment Model update")
# print_benchmark(trialmin)
# # @profview Wflow.update(model)

# Wflow.close_files(model, delete_output = true)
