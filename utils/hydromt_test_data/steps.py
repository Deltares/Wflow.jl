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
from hydromt_wflow.utils import get_wflow_var_fullname, set_config

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


def basemaps(model, grid_bounds: list[float] | None = None, **kwargs):
    """Run ``setup_basemaps``, optionally on a model grid with a fixed extent.

    HydroMT aligns the model grid to multiples of the model resolution. Some test models
    were built with an older HydroMT version that did not, so their grid is shifted by a
    fraction of a cell. With ``grid_bounds`` ([west, south, east, north], on the
    hydrography grid) the hydrography is clipped to these bounds, extended by whole model
    cells to contain the basin, before upscaling. HydroMT clips the upscaled grid to the
    active cells, which is extended to ``grid_bounds`` afterwards (``_pad_staticmaps``).
    """
    if grid_bounds is None:
        model.setup_basemaps(**kwargs)
        return

    from hydromt_wflow import workflows

    parse_region = workflows.parse_region

    def parse_region_fixed_extent(data_catalog, hydrography_fn, resolution, **region_kwargs):
        geometries, xy, ds_aligned = parse_region(
            data_catalog=data_catalog, hydrography_fn=hydrography_fn, resolution=resolution, **region_kwargs
        )
        if "meta_basins_highres" not in geometries:
            raise ValueError("grid_bounds requires a model resolution coarser than the hydrography")
        # extend grid_bounds by whole model cells (keeps the grid) to contain the basin
        basin = ds_aligned.raster.clip_mask(ds_aligned["mask"])
        bw, bs, be, bn = basin.raster.bounds
        west, south, east, north = grid_bounds

        def n_cells(overshoot):
            return int(np.ceil(max(overshoot, 0.0) / resolution - 1e-6)) + 1

        extent = (
            west - n_cells(west - bw) * resolution,
            south - n_cells(south - bs) * resolution,
            east + n_cells(be - east) * resolution,
            north + n_cells(bn - north) * resolution,
        )
        ds = data_catalog.get_rasterdataset(hydrography_fn)
        margin = abs(ds.raster.res[0]) / 4
        ds = ds.raster.clip_bbox((extent[0] + margin, extent[1] + margin, extent[2] - margin, extent[3] - margin))
        ds.coords["mask"] = ds.raster.geometry_mask(geometries["meta_basins_highres"])
        if int(ds["mask"].sum()) != int(ds_aligned["mask"].sum()):
            raise ValueError(f"The hydrography extent {extent} does not contain the basin")
        return geometries, xy, ds

    workflows.parse_region = parse_region_fixed_extent
    try:
        model.setup_basemaps(**kwargs)
    finally:
        workflows.parse_region = parse_region
    _pad_staticmaps(model, grid_bounds)


def _pad_staticmaps(model, grid_bounds: list[float]):
    """Extend the staticmaps to ``grid_bounds`` with nodata.

    HydroMT clips the model grid to the active cells, which can remove border rows or
    columns of the original grid.
    """
    ds = model.staticmaps.data
    x_dim, y_dim = ds.raster.x_dim, ds.raster.y_dim
    res = abs(ds.raster.res[0])
    west, south, east, north = grid_bounds
    xs = west + res * (np.arange(round((east - west) / res)) + 0.5)
    ys = north - res * (np.arange(round((north - south) / res)) + 0.5)
    if ds.raster.res[1] > 0:
        ys = ys[::-1]
    w, s, e, n = ds.raster.bounds
    tolerance = res / 10
    if w < west - tolerance or s < south - tolerance or e > east + tolerance or n > north + tolerance:
        raise ValueError(
            f"The model grid {ds.raster.bounds} exceeds grid_bounds {grid_bounds}: the "
            "upscaled basin has active cells outside grid_bounds"
        )
    if ds.raster.shape == (ys.size, xs.size):
        return

    def fill_value(da):
        if da.dtype == bool:
            return False
        if da.raster.nodata is not None:
            return da.raster.nodata
        return np.nan if np.issubdtype(da.dtype, np.floating) else 0

    fill = {name: fill_value(da) for name, da in ds.data_vars.items()}
    padded = ds.reindex({x_dim: xs, y_dim: ys}, method="nearest", tolerance=tolerance, fill_value=fill)
    logger.info(f"Extending the model grid from {ds.raster.shape} to {padded.raster.shape} cells (grid_bounds)")
    model.staticmaps._data = padded


def gauges(model, gauges_fn: str, **kwargs):
    """Run ``setup_gauges`` where the lowest gauge ID wins when gauges share a cell.

    ``setup_gauges`` writes the gauges in the order of the data source, so the last gauge
    of a cell is kept. The existing test models kept the gauge with the lowest ID.
    """
    from hydromt.error import NoDataStrategy

    gdf = model.data_catalog.get_geodataframe(gauges_fn, geom=model.basins, handle_nodata=NoDataStrategy.IGNORE)
    if gdf is None:
        logger.info(f"No {gauges_fn} gauges found within the model domain")
        return
    kwargs.setdefault("basename", Path(gauges_fn).stem.replace("_", "-"))
    model.setup_gauges(gauges_fn=gdf.sort_index(ascending=False), **kwargs)


def constant_maps(model, **values):
    """Add maps with a constant value on all active cells.

    HydroMT-Wflow v1 ``setup_constant_pars`` writes constants as values in the TOML file,
    while the Wflow tests read these parameters from the staticmaps. The maps are named
    after the Wflow standard name and renamed on export (see ``naming/*.yml``).
    """
    mask = _active_cells(model)
    for std_name, value in values.items():
        da = xr.where(mask, np.float32(value), np.float32(np.nan)).astype(np.float32)
        da.raster.set_nodata(np.nan)
        _set_static(model, da, std_name, std_name)


def outlet_maps(model, **values):
    """Add maps with a constant value at the basin outlet(s) and nodata elsewhere."""
    ldd = _grid(model, "basin__local_drain_direction")
    outlets = _active_cells(model) & (ldd == 5)
    for std_name, value in values.items():
        da = xr.where(outlets, np.float32(value), np.float32(np.nan)).astype(np.float32)
        da.raster.set_nodata(np.nan)
        _set_static(model, da, std_name, std_name)


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


ERA5_CLOUD_SOURCES = {
    # DestinE Earth Data Hub, chunks of 1080 h x 60 x 60 (45 days x 15 deg x 15 deg),
    # authentication with an API key in ~/.netrc (https://earthdatahub.destine.eu)
    "earthdatahub": {
        "url": "https://data.earthdatahub.destine.eu/era5/era5-single-levels-atmosphere-v0.zarr",
        "storage_options": {"client_kwargs": {"trust_env": True}},
        "time": "valid_time",
        "variables": {"u10": "u10", "v10": "v10", "ssr": "ssr", "str": "str"},
    },
    # Google ARCO-ERA5, chunks of 1 h x the whole globe (721 x 1440), no authentication
    "arco": {
        "url": "gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3",
        "storage_options": {"token": "anon"},
        "time": "time",
        "variables": {
            "u10": "10m_u_component_of_wind",
            "v10": "10m_v_component_of_wind",
            "ssr": "surface_net_solar_radiation",
            "str": "surface_net_thermal_radiation",
        },
    },
}
STEFAN_BOLTZMANN = 5.670374419e-8  # [W m-2 K-4]
DAILY_RIGHT_LABELED = {"time": "1D", "closed": "right", "label": "right"}


def _era5_daily_cloud(model, start, end, source: str, url: str | None, buffer: float, num_workers: int):
    """Daily wind speed and net radiation (ssr + str) from a cloud ERA5 Zarr store."""
    import dask

    settings = ERA5_CLOUD_SOURCES[source]
    url = url or settings["url"]
    names = settings["variables"]
    xmin, ymin, xmax, ymax = model.staticmaps.data.raster.bounds
    try:
        ds = xr.open_zarr(url, storage_options=settings["storage_options"], chunks={})
    except Exception as e:
        hint = (
            " Add your Earth Data Hub API key to ~/.netrc, see "
            "https://earthdatahub.destine.eu/getting-started#configuring-netrc"
            if source == "earthdatahub"
            else ""
        )
        raise RuntimeError(f"Could not open {url}.{hint}") from e
    ds = ds[list(names.values())].rename({v: k for k, v in names.items()})
    ds = ds.rename({settings["time"]: "time"}) if settings["time"] != "time" else ds
    # ERA5 longitudes are in [0, 360), select before converting to [-180, 180)
    lon_min, lon_max = (xmin - buffer) % 360, (xmax + buffer) % 360
    if lon_min <= lon_max:
        ds = ds.sel(longitude=slice(lon_min, lon_max))
    else:
        ds = xr.concat(
            [ds.sel(longitude=slice(lon_min, 360)), ds.sel(longitude=slice(0, lon_max))],
            dim="longitude",
        )
    ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180)
    descending = bool(ds.latitude[0] > ds.latitude[-1])
    lat_slice = slice(ymax + buffer, ymin - buffer) if descending else slice(ymin - buffer, ymax + buffer)
    ds = ds.sel(time=slice(start - pd.Timedelta(hours=23), end), latitude=lat_slice)
    wind = np.hypot(ds["u10"], ds["v10"]).resample(**DAILY_RIGHT_LABELED).mean()
    radiation = (ds["ssr"] + ds["str"]).resample(**DAILY_RIGHT_LABELED).sum() / 86400.0
    daily = xr.Dataset({"wind": wind, "net_radiation": radiation}).sel(time=slice(start, end))
    n_chunks = sum(int(np.prod(ds[v].data.numblocks)) for v in names)
    logger.info(f"Downloading {ds.sizes['time']} hours of ERA5 from {url} ({n_chunks} chunks)")
    with dask.config.set(scheduler="threads", num_workers=num_workers):
        daily = daily.load()
    daily.raster.set_crs(4326)
    daily = daily.raster.reproject_like(model.staticmaps.data, method="nearest_index")
    daily["net_radiation"].attrs.update(long_name="Net radiation (ssr + str)")
    return daily


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
    start: str | None = None,
    end: str | None = None,
    source: str = "catalog",
    hourly_source: str = "era5_hourly_zarr",
    ssr_source: str | None = None,
    albedo: float = 0.23,
    url: str | None = None,
    filename: str = "forcing_lst.nc",
    buffer: float = 0.5,
    num_workers: int = 16,
):
    """Add daily wind speed and net radiation forcing from hourly ERA5.

    Daily values are right-labeled, consistent with the other Wflow forcing: the value at
    day D covers the hourly ERA5 values in (D - 1 day, D]. Wind speed is the daily mean
    of the hourly 10 m wind speed [m s-1]. Net radiation is the daily mean net radiation
    [W m-2]. By default the model period is used (``start``, ``end``).

    Sources:

    * ``catalog`` (default): hourly ERA5 from the HydroMT data catalog (``hourly_source``
      for wind, temperature, dew point and shortwave radiation) with FAO-56 net radiation,
      see ``_era5_daily_catalog``. Net shortwave radiation uses ``albedo``, or ERA5 ssr
      from ``ssr_source`` (e.g. ``era5_hourly``, slow on the Deltares P drive). The
      Deltares P drive has no ERA5 net thermal radiation.
    * ``earthdatahub``: ERA5 ssr + str from the DestinE Earth Data Hub, chunked for
      regional time series, requires an API key in ~/.netrc.
    * ``arco``: ERA5 ssr + str from Google ARCO-ERA5 (no key). The store is chunked per
      hour for the whole globe, so each hour costs four global chunks of ~4 MB regardless
      of the model size (~70 GB for a year): use a short period with this source.

    The result is written to ``<model root>/<filename>`` instead of the model forcing,
    because HydroMT-Wflow does not write forcing with missing values (outside the period).
    It is merged into the model forcing on export (``extra_forcing`` in the workflow
    ``export`` section).
    """
    start = pd.Timestamp(start or model.config.get_value("time.starttime"))
    end = pd.Timestamp(end or model.config.get_value("time.endtime"))
    if source == "catalog":
        daily = _era5_daily_catalog(model, start, end, hourly_source, ssr_source, albedo)
    else:
        daily = _era5_daily_cloud(model, start, end, source, url, buffer, num_workers)
    daily["wind"].attrs.update(unit="m s-1")
    daily["net_radiation"].attrs.update(unit="W m-2")
    daily = daily.astype(np.float32)
    path = Path(model.root.path, filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    daily.to_netcdf(path)
    std_names = {
        "wind": "land_surface_air_flow__speed",
        "net_radiation": "land_surface_total_net_radiation__energy_flux",
    }
    for name, std_name in std_names.items():
        model.config.set(f"input.forcing.{std_name}", name)


def merge_forcing(forcing: xr.Dataset, path: Path) -> xr.Dataset:
    """Add the variables of a forcing file on the same grid, NaN outside its period."""
    with xr.open_dataset(path) as ds:
        extra = ds.load()
    if extra.raster.crs is None:
        extra.raster.set_crs(4326)
    renames = {extra.raster.x_dim: forcing.raster.x_dim, extra.raster.y_dim: forcing.raster.y_dim}
    extra = extra.rename({k: v for k, v in renames.items() if k != v})
    for dim in (forcing.raster.x_dim, forcing.raster.y_dim):
        if extra.sizes[dim] != forcing.sizes[dim] or not np.allclose(extra[dim], forcing[dim]):
            raise ValueError(f"{path} is not on the model forcing grid ({dim})")
        extra[dim] = forcing[dim]
    extra = extra.drop_vars("spatial_ref", errors="ignore").reindex(time=forcing["time"])
    return forcing.assign({name: extra[name] for name in extra.data_vars})


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
            if new_name not in data:
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
