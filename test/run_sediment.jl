using Wflow
tomlpath = joinpath(@__DIR__, "test", "sediment_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sediment_model(config)

model = Wflow.update(model)
#flush(model.writer.csv_io)

Wflow.close_files(model)
