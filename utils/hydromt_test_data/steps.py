"""Pipeline steps that complement the HydroMT model methods in the workflow files.

Steps in a workflow file prefixed with ``pipeline.`` are resolved to the functions in this
module. They receive the HydroMT model as first argument, followed by the keyword
arguments from the workflow file. All other steps are HydroMT(-Wflow) model methods.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from hydromt_wflow.utils import get_config, get_wflow_var_fullname, set_config

logger = logging.getLogger("hydromt_test_data")

INPUT_SECTIONS = ("forcing", "static", "cyclic")
# HydroMT-Wflow standard names that differ from the ones in this Wflow version
STANDARD_NAME_ALIASES = {"river__slope": "river_bottom_surface__slope"}


def map_name(model, std_name: str) -> str:
    """Return the netCDF variable name linked to a Wflow standard name in the config."""
    key = get_wflow_var_fullname(std_name, model.config.data)
    if key == std_name:
        raise KeyError(f"`{std_name}` is not set in the model config")
    value = model.config.get_value(key)
    return value["netcdf_variable_name"] if isinstance(value, dict) else value


def _grid(model, std_name: str) -> xr.DataArray:
    return model.staticmaps.data[map_name(model, std_name)]


def _active_cells(model) -> xr.DataArray:
    basins = _grid(model, "subbasin_location__count")
    return (basins != basins.raster.nodata) & (basins > 0)


def _set_static(model, da: xr.DataArray, name: str, std_name=None, section="static"):
    """Add a map to the staticmaps and (optionally) link it to a Wflow standard name."""
    model.staticmaps.set(da, name=name)
    if std_name is not None:
        key = "input" if section == "input" else f"input.{section}"
        model.config.set(f"{key}.{std_name}", name)


def _nearest_cell(da: xr.DataArray, x: float, y: float, candidates=None):
    """Return the (row, col) index of the cell nearest to (x, y), optionally restricted."""
    xs = da[da.raster.x_dim].values
    ys = da[da.raster.y_dim].values
    xx, yy = np.meshgrid(xs, ys)
    distance = (xx - x) ** 2 + (yy - y) ** 2
    if candidates is not None:
        distance = np.where(candidates, distance, np.inf)
    row, col = np.unravel_index(np.argmin(distance), distance.shape)
    if not np.isfinite(distance[row, col]):
        raise ValueError(f"No candidate cell found near ({x}, {y})")
    return row, col


def _snap_candidates(model, snap_to):
    if snap_to is None:
        return _active_cells(model).values
    if snap_to == "river":
        return (_grid(model, "river_location__mask") == 1).values
    if snap_to == "reservoir_outlet":
        if get_wflow_var_fullname("reservoir_location__count", model.config.data) == "reservoir_location__count":
            raise ValueError(
                "The model has no reservoirs to snap to, check the output of the "
                "setup_reservoirs_* steps (e.g. min_area or the reservoirs within the basin)"
            )
        outlets = _grid(model, "reservoir_location__count")
        return ((outlets != outlets.raster.nodata) & (outlets > 0)).values
    raise ValueError(f"Unknown snap_to option `{snap_to}`")


def stack_layers(model, std_name: str, sources: list[str]):
    """Stack 2D maps into a layered map, e.g. ksat at depth into layered ksat."""
    layers = model.staticmaps.data["layer"].values if "layer" in model.staticmaps.data.dims else None
    da = xr.concat([model.staticmaps.data[source] for source in sources], dim="layer")
    da = da.assign_coords(layer=layers if layers is not None and len(layers) == len(sources) else np.arange(1, len(sources) + 1))
    da.raster.set_nodata(np.nan)
    _set_static(model, da.astype(np.float32), std_name, std_name)


def point_map(
    model,
    std_name: str,
    points: list[list[float]],
    section: str = "input",
    snap_to: str | None = None,
    nodata: int = -9999,
):
    """Add an integer map with values at given [x, y, value] points and nodata elsewhere."""
    template = _grid(model, "subbasin_location__count")
    data = np.full(template.shape, nodata, dtype=np.int32)
    candidates = _snap_candidates(model, snap_to)
    for x, y, value in points:
        row, col = _nearest_cell(template, x, y, candidates)
        data[row, col] = value
    da = template.copy(data=data).astype(np.int32)
    da.raster.set_nodata(nodata)
    _set_static(model, da, std_name, std_name, section)


def point_timeseries(
    model,
    std_name: str,
    x: float,
    y: float,
    values: list[float] | None = None,
    segments: list[list[float]] | None = None,
    length: int | None = None,
    snap_to: str | None = None,
    background: str = "zero",
):
    """Add a cyclic map with a time series at one (snapped) point.

    The time series at the point is given by ``values`` or by ``segments``, a list of
    ``[start_index, value]`` pairs that hold until the next start index (``length`` sets
    the total number of time steps). Other active cells are zero (``background="zero"``)
    or nodata (``background="nodata"``).
    """
    if values is None:
        if segments is None or length is None:
            raise ValueError("Provide `values` or `segments` and `length`")
        series = np.zeros(length, dtype=np.float32)
        starts = [int(start) for start, _ in segments] + [length]
        for (start, value), end in zip(segments, starts[1:]):
            series[int(start) : end] = value
    else:
        series = np.asarray(values, dtype=np.float32)
    mask = _active_cells(model)
    row, col = _nearest_cell(mask, x, y, _snap_candidates(model, snap_to))
    fill = np.float32(0.0) if background == "zero" else np.float32(np.nan)
    data = np.where(mask.values, fill, np.float32(np.nan))
    data = np.repeat(data[None, :, :], series.size, axis=0)
    data[:, row, col] = series
    # reuse the coordinates of an existing cyclic dimension of the same length
    time = np.arange(1, series.size + 1)
    staticmaps = model.staticmaps.data
    for dim, size in staticmaps.sizes.items():
        if dim.startswith("time") and size == series.size and dim in staticmaps.coords:
            time = staticmaps[dim].values
            break
    da = xr.DataArray(data, coords={"time": time, **mask.coords}, dims=("time", *mask.dims))
    da.raster.set_nodata(np.nan)
    _set_static(model, da, std_name, std_name, section="cyclic")


def lst_parameters(model, mapping_fn: str, landuse_map: str = "meta_landuse"):
    """Add land surface temperature parameters by a lookup on the model land use map."""
    table = pd.read_csv(mapping_fn, index_col=0)
    landuse = model.staticmaps.data[landuse_map]
    mask = _active_cells(model).values
    missing = sorted(set(np.unique(landuse.values[mask])) - set(table.index))
    if missing:
        logger.warning(f"Land use classes {missing} missing in {mapping_fn}, set to NaN")
    std_names = {
        "z0m": "atmosphere_bottom_air_flow__roughness_length",
        "z0h": "atmosphere_bottom_air_heat_flow__roughness_length",
        "d0": "atmosphere_bottom_air_flow__zero_plane_displacement_length",
        "skin_layer_height": "atmosphere_bottom_air__above_ground_height",
    }
    for column, std_name in std_names.items():
        data = np.full(landuse.shape, np.nan, dtype=np.float32)
        for landuse_class, value in table[column].items():
            data[landuse.values == landuse_class] = value
        da = landuse.copy(data=data)
        da.raster.set_nodata(np.nan)
        _set_static(model, da, std_name, std_name)


STEFAN_BOLTZMANN = 5.670374419e-8  # [W m-2 K-4]
DAILY_RIGHT_LABELED = {"time": "1D", "closed": "right", "label": "right"}


def _era5_daily_catalog(model, start, end, hourly_source: str, ssr_source: str | None, albedo: float):
    """Daily wind speed and net radiation from hourly ERA5 in the HydroMT data catalog.

    Net radiation follows FAO-56 (Allen et al., 1998), Rn = Rns - Rnl:

    * Rns, net shortwave radiation: (1 - ``albedo``) times the incoming shortwave
      radiation (ssrd), or the daily mean of hourly ERA5 ssr from ``ssr_source``. The
      hourly ssr files on the Deltares P drive are unchunked global arrays, so reading a
      regional time series from them is slow.
    * Rnl, net longwave radiation (eq. 39): from daily maximum and minimum air temperature,
      actual vapour pressure (from the dew point temperature) and the ratio of incoming
      shortwave radiation to clear-sky radiation (0.75 + 2e-5 z) * Ra, with Ra the top of
      atmosphere incident solar radiation and z the model land elevation. ERA5 net thermal
      radiation (str) is not available on the Deltares P drive.
    """
    import time

    time_range = ((start - pd.Timedelta(hours=23)).isoformat(), end.isoformat())
    kwargs = {"geom": model.region, "buffer": 2, "time_range": time_range, "single_var_as_array": False}
    catalog = model.data_catalog
    hourly = catalog.get_rasterdataset(
        hourly_source, variables=["wind10_u", "wind10_v", "temp", "temp_dew", "kin", "kout"], **kwargs
    )
    temp_dew = hourly["temp_dew"]
    vapour_pressure = 0.6108 * np.exp(17.27 * temp_dew / (temp_dew + 237.3))  # [kPa]
    resample = DAILY_RIGHT_LABELED
    daily = xr.Dataset(
        {
            "wind": np.hypot(hourly["wind10_u"], hourly["wind10_v"]).resample(**resample).mean(),
            "temp_max": hourly["temp"].resample(**resample).max(),
            "temp_min": hourly["temp"].resample(**resample).min(),
            "vapour_pressure": vapour_pressure.resample(**resample).mean(),
            "kin": hourly["kin"].resample(**resample).mean(),
            "kout": hourly["kout"].resample(**resample).mean(),
        }
    ).sel(time=slice(start, end))
    tic = time.perf_counter()
    logger.info(f"Reading {hourly.sizes['time']} hours of ERA5 from {hourly_source}")
    daily = daily.load()
    logger.info(f"Read {hourly_source} in {time.perf_counter() - tic:.0f} s")
    if ssr_source is None:
        daily["net_shortwave"] = (1.0 - albedo) * daily["kin"]
        shortwave_name = f"(1 - {albedo}) ssrd"
    else:
        tic = time.perf_counter()
        logger.info(f"Reading {hourly.sizes['time']} hours of ERA5 ssr from {ssr_source}")
        ssr = catalog.get_rasterdataset(ssr_source, variables=["ssr"], **kwargs)["ssr"]
        daily["net_shortwave"] = ssr.resample(**resample).mean().sel(time=slice(start, end)).load()
        logger.info(f"Read {ssr_source} in {time.perf_counter() - tic:.0f} s")
        shortwave_name = "ssr"
    daily = daily.raster.reproject_like(model.staticmaps.data, method="nearest_index")

    elevation = _grid(model, "land_surface__elevation").fillna(0.0)
    clear_sky = (0.75 + 2e-5 * elevation) * daily["kout"]
    relative_shortwave = (daily["kin"] / clear_sky.where(clear_sky > 0)).clip(0.3, 1.0).fillna(0.3)
    temp_max_k, temp_min_k = daily["temp_max"] + 273.15, daily["temp_min"] + 273.15
    net_longwave = (
        STEFAN_BOLTZMANN
        * 0.5
        * (temp_max_k**4 + temp_min_k**4)
        * (0.34 - 0.14 * np.sqrt(daily["vapour_pressure"]))
        * (1.35 * relative_shortwave - 0.35)
    )
    radiation = daily["net_shortwave"] - net_longwave
    radiation.attrs.update(long_name=f"Net radiation ({shortwave_name} - FAO-56 net longwave radiation)")
    return xr.Dataset({"wind": daily["wind"], "net_radiation": radiation})


def era5_lst_forcing(
    model,
    hourly_source: str = "era5_hourly_zarr",
    ssr_source: str | None = None,
    albedo: float = 0.23,
):
    """Add daily wind speed and net radiation forcing from hourly ERA5 (model period).

    Daily values are right-labeled, consistent with the other Wflow forcing: the value at
    day D covers the hourly ERA5 values in (D - 1 day, D]. Wind speed is the daily mean
    of the hourly 10 m wind speed [m s-1]. Net radiation is the daily mean FAO-56 net
    radiation [W m-2], see ``_era5_daily_catalog``.
    """
    start = pd.Timestamp(model.config.get_value("time.starttime"))
    end = pd.Timestamp(model.config.get_value("time.endtime"))
    daily = _era5_daily_catalog(model, start, end, hourly_source, ssr_source, albedo)
    daily["wind"].attrs.update(unit="m s-1")
    daily["net_radiation"].attrs.update(unit="W m-2")
    model.forcing.set(daily.astype(np.float32))
    std_names = {
        "wind": "land_surface_air_flow__speed",
        "net_radiation": "land_surface_total_net_radiation__energy_flux",
    }
    for name, std_name in std_names.items():
        model.config.set(f"input.forcing.{std_name}", name)


def groundwater_maps(
    model,
    ksat_horizontal_ratio: float = 100.0,
    riverbed_thickness_fraction: float = 0.5,
    max_soil_thickness: float = 2000.0,
    aquifer_thickness: float = 5000.0,
    constant_head_depth: float = 0.5,
):
    """Add the maps used by the `sbm_gwf` model type, derived from the sbm maps.

    * river_bottom__elevation: land elevation minus bankfull river depth [m]
    * subsurface_surface_water__horizontal_saturated_hydraulic_conductivity:
      vertical ksat at the surface [mm d-1 -> m d-1] times ``ksat_horizontal_ratio``
    * subsurface__horizontal_saturated_hydraulic_conductivity_scale_parameter:
      ksat scale parameter f [mm-1 -> m-1]
    * river_water__infiltration_conductance and river_water__exfiltration_conductance:
      river area times ksat [m d-1] divided by the riverbed thickness, which is
      ``riverbed_thickness_fraction`` times the bankfull river depth [m2 d-1]
    * model_constant_boundary_condition__hydraulic_head: land elevation minus
      ``constant_head_depth`` at the basin outlet(s) [m]
    * SoilThickness_gw: soil thickness, extended to ``aquifer_thickness`` where the soil
      thickness equals ``max_soil_thickness`` [mm]
    * specific_yield: saturated minus residual water content [-]
    """
    elevation = _grid(model, "land_surface__elevation")
    river = _grid(model, "river_location__mask") == 1
    river_depth = _grid(model, "river_bank_water__depth")
    river_width = _grid(model, "river__width")
    river_length = _grid(model, "river__length")
    ksat = _grid(model, "soil_surface_water__vertical_saturated_hydraulic_conductivity")
    f = _grid(model, "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter")
    soil_thickness = _grid(model, "soil__thickness")
    theta_s = _grid(model, "soil_water__saturated_volume_fraction")
    theta_r = _grid(model, "soil_water__residual_volume_fraction")
    ldd = _grid(model, "basin__local_drain_direction")
    mask = _active_cells(model)

    def masked(da, where):
        out = da.where(where).astype(np.float32)
        out.raster.set_nodata(np.nan)
        return out

    ksat_m = ksat / 1000.0
    conductance = river_width * river_length * ksat_m / (riverbed_thickness_fraction * river_depth)
    maps = {
        "zb_river": ("river_bottom__elevation", masked(elevation - river_depth, river)),
        "kh_surface": (
            "subsurface_surface_water__horizontal_saturated_hydraulic_conductivity",
            masked(ksat_m * ksat_horizontal_ratio, mask),
        ),
        "gwf_f": (
            "subsurface__horizontal_saturated_hydraulic_conductivity_scale_parameter",
            masked(f * 1000.0, mask),
        ),
        "riverbed_cond": ("river_water__infiltration_conductance", masked(conductance, river)),
        "constant_head": (
            "model_constant_boundary_condition__hydraulic_head",
            masked(elevation - constant_head_depth, mask & (ldd == 5)),
        ),
        "SoilThickness_gw": (
            None,
            masked(xr.where(soil_thickness >= max_soil_thickness, aquifer_thickness, soil_thickness), mask),
        ),
        "specific_yield": (None, masked(theta_s - theta_r, mask)),
    }
    for name, (std_name, da) in maps.items():
        _set_static(model, da, name, std_name)
    model.config.set("input.static.river_water__exfiltration_conductance", "riverbed_cond")


# Export helpers


def normalize_config(config: dict) -> dict:
    """Rename config entries of HydroMT-Wflow standard names that differ in this Wflow version."""
    for old, new in STANDARD_NAME_ALIASES.items():
        key = get_wflow_var_fullname(old, config)
        if key != old:
            set_config(config, f"{key.rsplit('.', 1)[0]}.{new}", _pop(config, key))
    return config


def apply_naming(config: dict, staticmaps: xr.Dataset, forcing: xr.Dataset, naming: dict):
    """Rename HydroMT output to the netCDF variable names used by the Wflow tests.

    ``naming`` holds a ``maps`` table (old name -> new name, for maps without a Wflow
    standard name) and per config section (``input``, ``static``, ``cyclic``, ``forcing``)
    a table of Wflow standard name -> netCDF variable name. Variables linked in the config
    are copied to the new name and the config is updated. Variables that are no longer
    referenced afterwards are dropped. Entries for which no data is found are reported,
    but kept in the config (e.g. forcing produced by a Wflow run).
    """
    staticmaps = staticmaps.rename({k: v for k, v in naming.get("maps", {}).items() if k in staticmaps})
    replaced = set()
    missing = []
    for section in ("input",) + INPUT_SECTIONS:
        data = forcing if section == "forcing" else staticmaps
        for std_name, new_name in naming.get(section, {}).items():
            key = get_wflow_var_fullname(std_name, config)
            current = None
            if key != std_name:
                if isinstance(get_config(config, key), dict) and "value" in get_config(config, key):
                    continue  # constant value (setup_constant_pars), not a map
                current = _pop(config, key)
                if isinstance(current, dict):
                    current = current.get("netcdf_variable_name")
            if current is not None and current != new_name and current in data:
                data[new_name] = data[current]
                replaced.add(current)
            target = "input" if section == "input" else f"input.{section}"
            set_config(config, f"{target}.{std_name}", new_name)
    # check after all renames, as several standard names can share a netCDF variable
    for section in ("input",) + INPUT_SECTIONS:
        data = forcing if section == "forcing" else staticmaps
        if section == "forcing" and len(data.data_vars) == 0:
            continue
        for std_name, new_name in naming.get(section, {}).items():
            value = get_config(config, get_wflow_var_fullname(std_name, config))
            if new_name not in data and not (isinstance(value, dict) and "value" in value):
                missing.append(f"{section}.{std_name} = {new_name}")
    referenced = set(_input_names(config))
    staticmaps = staticmaps.drop_vars([v for v in replaced if v in staticmaps and v not in referenced])
    forcing = forcing.drop_vars([v for v in replaced if v in forcing and v not in referenced])
    return config, staticmaps, forcing, missing


def _pop(config: dict, key: str):
    """Remove a dotted key from a nested config and return its value."""
    *parents, last = key.split(".")
    table = config
    for part in parents:
        table = table[part]
    return table.pop(last)


def _input_names(config: dict):
    for key, value in config.get("input", {}).items():
        if isinstance(value, dict) and key in INPUT_SECTIONS:
            for entry in value.values():
                yield entry.get("netcdf_variable_name") if isinstance(entry, dict) else entry
        elif isinstance(value, str) and not key.startswith("path_"):
            yield value


# Post-processing of Wflow run output


def duplicate_variables(path: Path, duplicates: dict[str, str]):
    """Add copies of variables to a netCDF file (new name: existing name)."""
    with xr.open_dataset(path) as ds:
        ds = ds.load()
    for new_name, name in duplicates.items():
        ds[new_name] = ds[name]
    ds.to_netcdf(path)


def lattometres(lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Length of a degree longitude and latitude [m] at latitude ``lat``, as in Wflow."""
    rad = np.deg2rad(lat)
    latlen = 111132.92 - 559.82 * np.cos(2 * rad) + 1.175 * np.cos(4 * rad) - 0.0023 * np.cos(6 * rad)
    longlen = 111412.84 * np.cos(rad) - 93.5 * np.cos(3 * rad) + 0.118 * np.cos(5 * rad)
    return longlen, latlen


def overland_flow_geometry(ldd: np.ndarray, lat: np.ndarray, cell_size: float, river: np.ndarray, river_width: np.ndarray, river_length: np.ndarray):
    """Overland flow width and length [m] as computed by Wflow (``LandParameters``).

    ``ldd`` and the river maps are 2D arrays (lat, lon), ``lat`` the latitude of the rows
    and ``cell_size`` the cell size [degree]. The flow length depends on the flow direction.
    The surface flow width is the cell width perpendicular to the flow direction, or for
    river cells the land area (cell area minus river area) divided by the flow length.
    """
    longlen, latlen = lattometres(lat)
    x_length = np.broadcast_to((longlen * cell_size)[:, None], ldd.shape)
    y_length = np.broadcast_to((latlen * cell_size)[:, None], ldd.shape)
    diagonal = np.hypot(x_length, y_length)
    north_south, east_west = np.isin(ldd, (2, 8)), np.isin(ldd, (4, 6))
    flow_length = np.where(north_south, y_length, np.where(east_west, x_length, diagonal))
    flow_width = np.where(north_south, x_length, np.where(east_west, y_length, x_length * y_length / diagonal))
    area = x_length * y_length
    river_fraction = np.where(river, np.minimum(np.nan_to_num(river_length * river_width) / area, 1.0), 0.0)
    surface_flow_width = np.where(river, (1.0 - river_fraction) * area / flow_length, flow_width)
    return surface_flow_width, flow_length


def rhine_kinematic_wave_maps(output_grid: Path, staticmaps: Path, config: dict, target: Path, names: list[str]):
    """Write the overland flow kinematic wave test data (staticmaps-rhine.nc).

    ``names`` are variables in the last time step of the Wflow grid output. The overland
    flow width (Bw) and length (DCL) are not available as Wflow output and are computed
    from the static maps (``overland_flow_geometry``). The variables are written without
    coordinates and fill values, with dimensions (x, y) and y ordered from north to south,
    as read by `Wflow/test/routing_process.jl`.
    """
    inputs = config["input"]

    def static_name(std_name):
        value = inputs.get(std_name, inputs.get("static", {}).get(std_name))
        return value["netcdf_variable_name"] if isinstance(value, dict) else value

    def dims(da):
        return next(d for d in da.dims if d in ("y", "lat", "latitude")), next(d for d in da.dims if d in ("x", "lon", "longitude"))

    def north_to_south(values: np.ndarray, lat: np.ndarray) -> np.ndarray:
        values = values[::-1] if lat[0] < lat[-1] else values
        return values.T

    with xr.open_dataset(staticmaps, mask_and_scale=False) as ds_static, xr.open_dataset(output_grid) as ds_out:
        ldd_da = ds_static[static_name("basin__local_drain_direction")]
        y_dim, x_dim = dims(ldd_da)
        lat, lon = ds_static[y_dim].values, ds_static[x_dim].values
        ldd = ldd_da.values.astype(np.uint8)
        river = ds_static[static_name("river_location__mask")].values
        river = np.isfinite(river) & (river == 1)
        river_width = ds_static[static_name("river__width")].where(river).values.astype(np.float64)
        river_length = ds_static[static_name("river__length")].where(river).values.astype(np.float64)
        arrays = {}
        for name in names:
            da = ds_out[name].isel(time=-1)
            out_y, out_x = dims(da)
            # Wflow writes y ascending, HydroMT staticmaps can be north-up
            da = da.transpose(out_y, out_x).sortby(out_y, ascending=bool(lat[0] < lat[-1]))
            if not (np.allclose(da[out_y], lat) and np.allclose(da[out_x], lon)):
                raise ValueError("Wflow output and staticmaps grids differ")
            arrays[name] = da.values.astype(np.float64)
    width, length = overland_flow_geometry(ldd, lat, float(np.abs(np.diff(lon)).mean()), river, river_width, river_length)
    inactive = ~np.isin(ldd, np.arange(1, 10))
    arrays["Bw"], arrays["DCL"] = np.where(inactive, np.nan, width), np.where(inactive, np.nan, length)
    ldd = np.where(inactive, 255, ldd).astype(np.uint8)
    out = xr.Dataset(
        {name: (("x", "y"), north_to_south(values, lat)) for name, values in {"ldd": ldd, **arrays}.items()}
    )
    encoding = {name: {"_FillValue": None} for name in out.data_vars}
    out.to_netcdf(target, encoding=encoding)
