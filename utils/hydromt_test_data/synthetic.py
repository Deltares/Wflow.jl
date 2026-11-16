"""Synthetic 7 x 5 cell test case for the `sbm_gwf` model type (no HydroMT model).

Writes `staticmaps-sbm-groundwater.nc` and `forcing-sbm-groundwater-part{1,2}.nc`. The
forcing values are read from `synthetic_sbm_gwf_forcing.csv`.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

X = 124488.0 + 50.0 * np.arange(5)
Y = 463139.0 + 50.0 * np.arange(7)

LDD = [
    [8, 8, 7, 8, 8],
    [6, 8, 9, 8, 4],
    [8, 8, 9, 8, 4],
    [5, 4, 4, 4, 4],
    [3, 2, 2, 2, 2],
    [2, 2, 6, 2, 4],
    [6, 2, 2, 2, 4],
]
RIVER = [
    [0, 0, 0, 0, 0],
    [0, 1, 0, 1, 0],
    [0, 1, 0, 1, 0],
    [1, 1, 1, 1, 1],
    [0, 1, 0, 1, 0],
    [0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0],
]
DEM = [
    [1.99, 1.937, 1.96, 1.948, 1.88],
    [1.95, 1.7, 1.95, 1.8, 1.99],
    [1.97, 1.675, 1.98, 1.775, 1.92],
    [1.6, 1.65, 1.7, 1.75, 1.8],
    [1.93, 1.675, 1.9, 1.775, 1.94],
    [1.86, 1.7, 1.99, 1.8, 1.9],
    [1.98, 1.97, 1.99, 1.9343, 1.94],
]
# Main channel (row 3) and tributaries
MAIN_CHANNEL_ROW = 3
INACTIVE_CELLS = [(0, 0), (6, 3), (6, 4)]
DRAIN_CELL = (3, 2)
CONSTANT_HEAD_COLUMN = 4

CONSTANTS = {
    "wflow_riverlength": 50.0,
    "soilthickness": 2000.0,
    "kv": 250.0,
    "f": 0.001,
    "M": 500.0,
    "thetaS": 0.45,
    "thetaR": 0.1,
    "rootingdepth": 300.0,
    "cmax": 1.0,
    "EoverR": 0.11,
    "canopygapfraction": 0.4,
    "PathFrac": 0.0,
    "InfiltCapSoil": 100.0,
    "Slope": 0.0005,
    "RiverSlope": 0.001,
    "N": 0.15,
    "N_river": 0.035,
    "k": 2.5,
}


def staticmaps() -> xr.Dataset:
    shape = (Y.size, X.size)
    river = np.array(RIVER, dtype=np.float64)
    is_river = river == 1
    dem = np.array(DEM, dtype=np.float32)
    river_bottom = dem - np.float32(0.5)
    main_channel = np.zeros(shape, dtype=bool)
    main_channel[MAIN_CHANNEL_ROW, :] = True
    width = np.where(is_river, 2.0, 0.0).astype(np.float32)
    width[MAIN_CHANNEL_ROW, ::2] = 4.0
    conductance = np.where(is_river, np.where(main_channel, 320.0, 160.0), 0.0).astype(np.float32)
    subcatch = np.ones(shape, dtype=np.int32)
    for cell in INACTIVE_CELLS:
        subcatch[cell] = -9999
    drain = np.zeros(shape)
    drain[DRAIN_CELL] = 1.0
    elev_drain = np.full(shape, np.nan, dtype=np.float32)
    elev_drain[DRAIN_CELL] = dem[DRAIN_CELL] - np.float32(0.1)
    cond_drain = np.full(shape, np.nan, dtype=np.float32)
    cond_drain[DRAIN_CELL] = 1.0
    constant_head = np.full(shape, np.nan, dtype=np.float32)
    constant_head[:, CONSTANT_HEAD_COLUMN] = 1.8

    maps = {
        "wflow_subcatch": subcatch,
        "wflow_ldd": np.array(LDD, dtype=np.int8),
        "wflow_river": river,
        "river_bottom": river_bottom,
        "drain": drain,
        "elev_drain": elev_drain,
        "cond_drain": cond_drain,
        "wflow_dem": dem,
        "wflow_riverwidth": width,
        **{name: np.full(shape, value, dtype=np.float32) for name, value in CONSTANTS.items()},
        "infiltration_conductance": conductance,
        "exfiltration_conductance": conductance.copy(),
        "specific_yield": np.full(shape, 0.35, dtype=np.float32),
        "constant_head": constant_head,
        "bankfull_elevation": np.where(is_river, dem.astype(np.float64), 0.0),
        "bankfull_depth": np.where(is_river, river_bottom, dem).astype(np.float64),
    }
    ds = xr.Dataset({name: (("y", "x"), data) for name, data in maps.items()}, coords={"x": X, "y": Y})
    ds["soillayers"] = ("layer", np.array([100.0, 300.0, 800.0, np.nan]))
    return ds


def encoding(ds: xr.Dataset) -> dict:
    enc = {}
    for name, da in ds.data_vars.items():
        if name == "wflow_subcatch":
            enc[name] = {"_FillValue": np.int32(-9999)}
        elif name == "wflow_ldd":
            enc[name] = {"_FillValue": np.int8(-15)}
        elif name in ("specific_yield", "bankfull_elevation", "bankfull_depth"):
            enc[name] = {"_FillValue": np.nan}
        else:
            enc[name] = {"_FillValue": da.dtype.type(-9999)}
    return enc


def forcing(forcing_csv: Path) -> dict[int, xr.Dataset]:
    table = pd.read_csv(forcing_csv)
    parts = {}
    for part, df in table.groupby("file"):
        times = sorted(df["time"].unique())
        data_vars = {}
        for var in ("P", "PET"):
            values = df[df["variable"] == var].sort_values("time").iloc[:, 3:].to_numpy()
            data_vars[var] = (("time", "y", "x"), values.reshape(len(times), Y.size, X.size))
        ds = xr.Dataset(data_vars, coords={"time": np.asarray(times, dtype=np.int64), "y": Y, "x": X})
        ds["time"].attrs.update(units="days since 2000-06-01", calendar="proleptic_gregorian")
        parts[int(part)] = ds
    return parts


def write(input_dir: Path, forcing_csv: Path):
    ds = staticmaps()
    ds.to_netcdf(input_dir / "staticmaps-sbm-groundwater.nc", encoding=encoding(ds))
    for part, ds_forcing in forcing(forcing_csv).items():
        enc = {var: {"_FillValue": np.nan} for var in ("P", "PET")}
        ds_forcing.to_netcdf(input_dir / f"forcing-sbm-groundwater-part{part}.nc", encoding=enc)
