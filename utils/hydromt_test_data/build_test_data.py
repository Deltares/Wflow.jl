"""Rebuild the Wflow test models (Wflow/test/*.toml and Wflow/test/data/input) with HydroMT.

Usage (from the repository root):

    pixi run -e hydromt build-test-data [STAGE ...] [--output DIR] [--only NAME ...]

Stages (default: all, in this order):

    build    HydroMT model builds/updates (`hydromt/*.yml`, `models.yml`) into `<output>/models`
    export   rename to the netCDF names of the tests and write `<output>/test/data/input`
    configs  write the test TOML files (`configs/*.yml`) and the Wflow run TOML files
    run      Wflow runs that produce the states, the sediment forcing and Rhine test data
    compare  compare the generated files with `Wflow/test`

See README.md for details.
"""

import argparse
import copy
import logging
import os
import shlex
import shutil
import subprocess
import sys
import tomllib
import urllib.request
from pathlib import Path

import cftime
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from hydromt.readers import read_workflow_yaml
from hydromt.writers import write_nc, write_toml
from hydromt_wflow import WflowSbmModel, WflowSedimentModel
from hydromt_wflow.utils import get_config, get_wflow_var_fullname, set_config

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import steps  # noqa: E402
import synthetic  # noqa: E402
from testdata import TestDataComponent  # noqa: E402

REPO = HERE.parents[1]
TEST_DIR = REPO / "Wflow" / "test"
STAGES = ["build", "export", "configs", "run", "compare"]
MODELS = {"wflow_sbm": WflowSbmModel, "wflow_sediment": WflowSedimentModel}
# Build order matters: moselle_sediment updates moselle_sbm
WORKFLOWS = ["moselle_sbm", "moselle_sediment", "piave_sbm", "rhine_sbm"]
# Run order matters: the sediment spin-up uses the forcing of the sediment forcing run
RUNS = [
    "moselle_spinup",
    "moselle_sediment_forcing",
    "moselle_sediment_spinup",
    "piave_spinup",
    "piave_gwf_spinup",
    "sbm_gwf_spinup",
    "rhine_kinematic_wave",
]
HYDROMT_WFLOW_TAG = "v1.1.0"
DEMAND_URL = (
    "https://raw.githubusercontent.com/Deltares/hydromt_wflow/"
    f"{HYDROMT_WFLOW_TAG}/examples/data/demand"
)
DEMAND_FILES = [
    "data_catalog.yml",
    "glcnmo.tif",
    "irrigated_area.tif",
    "lisflood_gwbodies.tif",
    "lisflood_gwfrac.tif",
    "lisflood_ncfrac.tif",
    "lisflood_waterregions.tif",
    "pcr_globwb.nc",
]
DEFAULT_JULIA = [
    "pixi",
    "run",
    "--environment",
    "default",
    "julia",
    "--startup-file=no",
    f"--project={REPO / 'Wflow'}",
]

logger = logging.getLogger("hydromt_test_data")


class Paths:
    def __init__(self, output: Path):
        self.output = output
        self.models = output / "models"
        self.staging = output / "staging"
        self.test = output / "test"
        self.input = self.test / "data" / "input"
        self.runs = output / "runs"
        self.cache = output / "cache"

    def placeholders(self) -> dict[str, str]:
        return {"{output}": self.output.as_posix()}


def load_yaml(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f) or {}


def resolve(obj, placeholders: dict[str, str]):
    """Replace placeholders in all strings of a (nested) object."""
    if isinstance(obj, str):
        for key, value in placeholders.items():
            obj = obj.replace(key, value)
        return obj
    if isinstance(obj, dict):
        return {k: resolve(v, placeholders) for k, v in obj.items()}
    if isinstance(obj, list):
        return [resolve(v, placeholders) for v in obj]
    return obj


def load_model_spec(name: str, paths: Paths) -> dict:
    """Pipeline settings of a model (models.yml): workflow, base model, data, export."""
    return resolve(load_yaml(HERE / "models.yml")[name], paths.placeholders())


# Stage: build


def download_demand_catalog(paths: Paths):
    """Download the water demand example data of HydroMT-Wflow (Piave region)."""
    folder = paths.cache / "hydromt_wflow_demand"
    folder.mkdir(parents=True, exist_ok=True)
    for filename in DEMAND_FILES:
        target = folder / filename
        if not target.exists():
            logger.info(f"Downloading {DEMAND_URL}/{filename}")
            urllib.request.urlretrieve(f"{DEMAND_URL}/{filename}", target)
    catalog = load_yaml(folder / "data_catalog.yml")
    # avoid a name clash with `pcr_globwb` in the Deltares data catalog
    if "pcr_globwb" in catalog:
        catalog["hydromt_wflow_example_pcr_globwb"] = catalog.pop("pcr_globwb")
        with open(folder / "data_catalog.yml", "w") as f:
            yaml.safe_dump(catalog, f, sort_keys=False)
    return (folder / "data_catalog.yml").as_posix()


def build(names: list[str], paths: Paths, data_libs: list[str]):
    """Build or update the models with HydroMT (``Model.build`` / ``Model.update``).

    The workflow files in ``hydromt/`` are plain HydroMT workflows. The custom steps
    (``testdata.*``) are provided by the ``testdata`` component that is added to the model.
    """
    for name in names:
        spec = load_model_spec(name, paths)
        modeltype, model_init, steps_ = read_workflow_yaml(HERE / "hydromt" / spec["workflow"])
        libs = data_libs + model_init.pop("data_libs", [])
        if "demand" in spec.get("data_libs", []):
            libs.append(download_demand_catalog(paths))
        root = paths.models / name
        logger.info(f"Building {name} in {root}")
        if "base" in spec:
            model = MODELS[modeltype](root=str(paths.models / spec["base"]), mode="r", data_libs=libs, **model_init)
            model.add_component("testdata", TestDataComponent(model))
            model.update(model_out=root, steps=steps_, forceful_overwrite=True)
        else:
            model = MODELS[modeltype](root=str(root), mode="w+", data_libs=libs, **model_init)
            model.add_component("testdata", TestDataComponent(model))
            model.build(steps=steps_)


# Stage: export


def write_staticmaps(ds: xr.Dataset, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    write_nc(ds, path, compress=True, gdal_compliant=True, rename_dims=True, force_overwrite=True)


def write_forcing(ds: xr.Dataset, path: Path):
    """Write forcing like HydroMT-Wflow (float32, 2 decimals, one time step per chunk)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    ds = ds.drop_vars(["mask", "idx_out"], errors="ignore").round(2)
    chunksizes = (1, ds.raster.ycoords.size, ds.raster.xcoords.size)
    encoding = {v: {"dtype": "float32", "chunksizes": chunksizes} for v in ds.data_vars}
    encoding["time"] = {"units": "days since 1900-01-01T00:00:00"}
    write_nc(
        ds,
        path,
        compress=True,
        gdal_compliant=True,
        rename_dims=True,
        force_overwrite=True,
        to_netcdf_kwargs={"encoding": encoding},
    )


def export_model(name: str, paths: Paths):
    settings = load_model_spec(name, paths).get("export", {})
    modeltype, model_init, _ = read_workflow_yaml(HERE / "hydromt" / load_model_spec(name, paths)["workflow"])
    model_init.pop("data_libs", None)
    model = MODELS[modeltype](root=str(paths.models / name), mode="r", **model_init)
    model.config.read()
    model.staticmaps.read()
    forcing = xr.Dataset()
    if settings.get("forcing"):
        model.forcing.read()
        forcing = model.forcing.data
    config = steps.normalize_config(copy.deepcopy(model.config.data))
    staticmaps = model.staticmaps.data
    if "naming" in settings:
        naming = load_yaml(HERE / "naming" / f"{settings['naming']}.yml")
        config, staticmaps, forcing, missing = steps.apply_naming(config, staticmaps, forcing, naming)
        for entry in missing:
            logger.warning(f"{name}: no data found for {entry}")
    if "staticmaps" in settings:
        logger.info(f"Writing {settings['staticmaps']}")
        write_staticmaps(staticmaps, paths.input / settings["staticmaps"])
        set_config(config, "input.path_static", settings["staticmaps"])
    for item in settings.get("forcing", []):
        target = paths.input / item["path"]
        logger.info(f"Writing {target}")
        write_forcing(forcing.sel(time=slice(str(item["start"]), str(item["end"]))), target)
    paths.staging.mkdir(parents=True, exist_ok=True)
    write_toml(paths.staging / f"{name}.toml", config)
    model.close()


def write_noleap_forcing(paths: Paths):
    """Write a 10 x 10 cell, 5 day subset of the Moselle forcing with a noleap calendar."""
    with xr.open_dataset(paths.input / "forcing-moselle.nc") as ds:
        y_dim, x_dim = ds.raster.y_dim, ds.raster.x_dim
        ix = int(np.abs(ds[x_dim].values - 6.251667).argmin())
        iy = int(np.abs(ds[y_dim].values - 49.5925).argmin())
        ds = ds[["precip", "pet", "temp"]].isel(
            time=slice(0, 5), **{x_dim: slice(ix, ix + 10), y_dim: slice(iy, iy + 10)}
        )
        ds = ds.rename({x_dim: "x", y_dim: "y"}).load()
    times = [cftime.DatetimeNoLeap(2000, 1, day) for day in range(2, 7)]
    ds = ds.assign_coords(time=times)
    encoding = {v: {"_FillValue": np.float32(np.nan), "dtype": "float32"} for v in ds.data_vars}
    # the _FillValue of the time dimension is on purpose (Wflow/test/io.jl)
    encoding["time"] = {"units": "days since 2000-01-02 00:00:00", "calendar": "noleap", "dtype": "float64", "_FillValue": np.nan}
    ds.to_netcdf(paths.input / "forcing-calendar-noleap.nc", encoding=encoding)


def write_reservoir_curves(paths: Paths):
    """Write the reservoir storage and rating curve test tables."""
    source = HERE / "reservoir_curves"
    for filename in ("reservoir_sh_1.csv", "reservoir_sh_2.csv"):
        shutil.copyfile(source / filename, paths.input / filename)
    curve = pd.read_csv(source / "reservoir_hq_2_curve.csv", dtype=str)
    days = [str(day) for day in range(1, 366)]
    lines = ["H," + ",".join(days)]
    lines += [h + "," + ",".join([q] * len(days)) for h, q in zip(curve["H"], curve["Q"])]
    (paths.input / "reservoir_hq_2.csv").write_text("\n".join(lines) + "\n", newline="\n")


def export(names: list[str], paths: Paths):
    paths.input.mkdir(parents=True, exist_ok=True)
    for name in names:
        export_model(name, paths)
    if (paths.input / "forcing-moselle.nc").exists():
        write_noleap_forcing(paths)
    synthetic.write(paths.input, HERE / "synthetic_sbm_gwf_forcing.csv")
    write_reservoir_curves(paths)


# Stage: configs


def load_spec(kind: str, name: str, paths: Paths) -> dict:
    """Load a TOML specification, merged with the specification it extends."""
    spec = resolve(load_yaml(HERE / kind / f"{name}.yml"), paths.placeholders())
    if "extends" not in spec:
        return spec
    base = load_spec("configs", spec.pop("extends"), paths)
    merged = {**base, **{k: v for k, v in spec.items() if k not in ("config", "remove")}}
    merged["config"] = {**base.get("config", {}), **spec.get("config", {})}
    merged["remove"] = base.get("remove", []) + spec.get("remove", [])
    return merged


def add_inputs(config: dict, inputs: dict[str, list[str]], staging: dict, name: str):
    for section, std_names in inputs.items():
        target = "input" if section == "input" else f"input.{section}"
        for std_name in std_names:
            key = get_wflow_var_fullname(std_name, staging)
            if key == std_name:
                raise KeyError(f"{name}: `{std_name}` not found in the exported model config")
            set_config(config, f"{target}.{std_name}", copy.deepcopy(get_config(staging, key)))


def make_config(kind: str, name: str, paths: Paths) -> tuple[dict, dict]:
    """Create a Wflow TOML configuration from a specification.

    The `inputs` are taken from the config of the exported HydroMT model (`model`), all
    other entries from `config` (dotted keys), followed by the removal of `remove` keys.
    The inputs are inserted after the `input.path_*` entries to follow the TOML layout.
    """
    spec = load_spec(kind, name, paths)
    staging = {}
    if spec.get("model"):
        with open(paths.staging / f"{spec['model']}.toml", "rb") as f:
            staging = tomllib.load(f)
    inputs = spec.get("inputs")
    config = copy.deepcopy(staging) if inputs == "all" else {}
    items = list(spec.get("config", {}).items())
    path_keys = [i for i, (key, _) in enumerate(items) if key.startswith("input.path_")]
    insert_at = path_keys[-1] if path_keys else len(items) - 1
    for i, (key, value) in enumerate(items):
        set_config(config, key, copy.deepcopy(value))
        if i == insert_at and isinstance(inputs, dict):
            add_inputs(config, inputs, staging, name)
    if not items and isinstance(inputs, dict):
        add_inputs(config, inputs, staging, name)
    for key in spec.get("remove", []):
        *parents, last = key.split(".")
        table = get_config(config, ".".join(parents)) if parents else config
        if isinstance(table, dict):
            table.pop(last, None)
    return config, spec


def write_configs(paths: Paths, runs: list[str]):
    paths.test.mkdir(parents=True, exist_ok=True)
    for spec_path in sorted((HERE / "configs").glob("*.yml")):
        config, _ = make_config("configs", spec_path.stem, paths)
        write_toml(paths.test / f"{spec_path.stem}.toml", config)
    paths.runs.mkdir(parents=True, exist_ok=True)
    for name in runs:
        config, _ = make_config("runs", name, paths)
        write_toml(paths.runs / f"{name}.toml", config)


# Stage: run


def run_wflow(names: list[str], paths: Paths, julia: list[str]):
    for name in names:
        toml_path = paths.runs / f"{name}.toml"
        logger.info(f"Running Wflow with {toml_path}")
        subprocess.run([*julia, "-e", "using Wflow; Wflow.run(ARGS[1])", str(toml_path)], check=True, cwd=REPO)
        config, spec = make_config("runs", name, paths)
        dir_output = Path(config["dir_output"])
        for product in spec.get("products", []):
            target = paths.input / product["to"]
            shutil.copyfile(dir_output / product["from"], target)
            if "duplicate" in product:
                steps.duplicate_variables(target, product["duplicate"])
        if "postprocess" in spec:
            settings = dict(spec["postprocess"])
            function = getattr(steps, settings.pop("function"))
            function(
                output_grid=dir_output / settings["output_grid"],
                staticmaps=Path(config["dir_input"]) / config["input"]["path_static"],
                config=config,
                target=paths.input / settings["target"],
                names=settings["names"],
            )


# Stage: compare


def flatten(d: dict, prefix: str = "") -> dict:
    out = {}
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            out.update(flatten(v, key))
        else:
            out[key] = v
    return out


def compare(paths: Paths):
    """Report differences between the generated and the existing test files."""
    print("\n# TOML files")
    for original in sorted(TEST_DIR.glob("*.toml")):
        if original.name == "Project.toml":
            continue
        generated = paths.test / original.name
        if not generated.exists():
            print(f"- {original.name}: not generated")
            continue
        with open(original, "rb") as f1, open(generated, "rb") as f2:
            a, b = flatten(tomllib.load(f1)), flatten(tomllib.load(f2))
        diffs = [f"only in original: {k}" for k in a.keys() - b.keys()]
        diffs += [f"only in generated: {k}" for k in b.keys() - a.keys()]
        diffs += [f"{k}: {a[k]!r} != {b[k]!r}" for k in a.keys() & b.keys() if a[k] != b[k]]
        print(f"- {original.name}: {'identical' if not diffs else ''}")
        for diff in sorted(diffs):
            print(f"    {diff}")

    print("\n# Input data")
    input_dir = TEST_DIR / "data" / "input"
    if not input_dir.exists():
        print("Wflow/test/data/input not found, run `pixi run download-test-data` first")
        return
    for original in sorted(input_dir.iterdir()):
        generated = paths.input / original.name
        if not generated.exists():
            print(f"- {original.name}: not generated")
            continue
        if original.suffix == ".csv":
            a, b = (p.read_bytes().replace(b"\r\n", b"\n").rstrip(b"\n") for p in (original, generated))
            print(f"- {original.name}: {'identical' if a == b else 'different'}")
            continue
        with xr.open_dataset(original, decode_times=False) as a, xr.open_dataset(generated, decode_times=False) as b:
            print(f"- {original.name}: dims {dict(a.sizes)} -> {dict(b.sizes)}")
            only_a, only_b = sorted(set(a.data_vars) - set(b.data_vars)), sorted(set(b.data_vars) - set(a.data_vars))
            if only_a:
                print(f"    missing variables: {only_a}")
            if only_b:
                print(f"    additional variables: {only_b}")
            for var in sorted(set(a.data_vars) & set(b.data_vars)):
                va, vb = a[var], b[var]
                if not (np.issubdtype(va.dtype, np.number) and np.issubdtype(vb.dtype, np.number)):
                    continue
                stats = [float(v.mean()) if v.size else np.nan for v in (va, vb)]
                if va.shape == vb.shape and np.allclose(va.values, vb.values, equal_nan=True):
                    print(f"    {var}: identical")
                else:
                    print(f"    {var}: shape {va.shape} -> {vb.shape}, mean {stats[0]:.6g} -> {stats[1]:.6g}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stages", nargs="*", metavar="STAGE", help=f"one or more of {', '.join(STAGES)} or all (default)")
    parser.add_argument("--output", type=Path, default=HERE / "output", help="output directory")
    parser.add_argument("--only", nargs="+", help="limit the build/export/run stages to these names")
    parser.add_argument(
        "--data-libs",
        nargs="+",
        default=["deltares_data"],
        help="HydroMT data catalogs (default: deltares_data)",
    )
    parser.add_argument(
        "--julia",
        help="command to start Julia with the Wflow project "
        f"(default: {' '.join(DEFAULT_JULIA)})",
    )
    args = parser.parse_args()
    invalid = [stage for stage in args.stages if stage not in STAGES + ["all"]]
    if invalid:
        parser.error(f"invalid stage(s) {invalid}, choose from {', '.join(STAGES)} or all")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    paths = Paths(args.output.resolve())
    stages = STAGES if not args.stages or "all" in args.stages else args.stages
    workflows = [n for n in WORKFLOWS if not args.only or n in args.only]
    runs = [n for n in RUNS if not args.only or n in args.only]
    julia = shlex.split(args.julia, posix=os.name != "nt") if args.julia else DEFAULT_JULIA
    if "build" in stages:
        build(workflows, paths, args.data_libs)
    if "export" in stages:
        export(workflows, paths)
    if "configs" in stages:
        write_configs(paths, RUNS)
    if "run" in stages:
        run_wflow(runs, paths, julia)
    if "compare" in stages:
        compare(paths)


if __name__ == "__main__":
    main()
