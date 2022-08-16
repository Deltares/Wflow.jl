# Step 5: Analyzing the model output

After running the model example from the previous step 4, the model results can be found in
`data/output_moselle_simple.csv`.

If required, it is also possible to output NetCDF files as output, by modifying the TOML
file. An example is shown below:

```toml
# Spatial output
[output]
path = "data/output.nc"

[output.lateral.river]
q_av = "q_river"

[output.lateral.land]
q = "q_land"
h = "h_land"

# Scalar output (mapped to the specified map)
[netcdf]
path = "data/output_scalar.nc"

[[netcdf.variable]]
name = "Q"
map = "gauges"
parameter = "lateral.river.q_av"

[[netcdf.variable]]
name = "prec"
map = "subcatchment"
parameter = "vertical.precipitation"
reducer = "mean"
```

Using your own preferred programming language, the model output files can be easily read and
visualized.