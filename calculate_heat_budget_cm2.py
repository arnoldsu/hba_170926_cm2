#!/usr/bin/env python3
"""Vialard et al. (2001)-inspired mixed-layer heat-budget decomposition."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import xarray as xr

RHO0 = 1035.0
CP = 3989.24495
EARTH_RADIUS = 6_371_000.0
MAX_MLD = 50.0
READ_DEPTH = 60.0

FILE_VARIABLE = {}
MODEL_TOKEN = "bj594_piControl"
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = PROJECT_ROOT / "data"

def open_var(root: Path, stem: str, chunks: int):
    file_stem, variable = FILE_VARIABLE.get(stem, (stem, stem))
    hits = sorted(root.glob(f"{file_stem}_Omon_*{MODEL_TOKEN}*.nc"))
    if not hits:
        raise FileNotFoundError(f"No native CM2 file for {stem} under {root}")
    path = hits[0]
    ds = xr.open_dataset(path)
    chunk_map = {d: (-1 if d != "time" else chunks) for d in ds.dims}
    ds = ds.chunk(chunk_map)
    da = ds[variable]
    rename = {}
    for old, new in (("xt_ocean", "lon"), ("xu_ocean", "lon"), ("yt_ocean", "lat"), ("yu_ocean", "lat")):
        if old in da.dims:
            rename[old] = new
    da = da.rename(rename)
    if "lon" in da.coords:
        da = da.assign_coords(lon=da.lon % 360).sortby("lon").sel(lon=slice(120, 280))
    if "lat" in da.coords:
        da = da.sel(lat=slice(-25, 25))
    for zname in ("st_ocean", "sw_ocean"):
        if zname in da.coords:
            da = da.sel({zname: slice(0, READ_DEPTH)})
    # Guard against duplicate timestamps before aligning all CM2 inputs.
    if "time" in da.dims:
        duplicate = da.get_index("time").duplicated(keep="first")
        if duplicate.any():
            print(f"{stem}: dropping {int(duplicate.sum())} duplicate time records")
            da = da.isel(time=np.flatnonzero(~duplicate))
    return da

def last_where(field, condition, zdim):
    """Deepest valid value satisfying condition, without eager 4-D indexing."""
    out = xr.full_like(field.isel({zdim: 0}, drop=True), np.nan)
    for k in range(field.sizes[zdim]):
        out = xr.where(condition.isel({zdim: k}), field.isel({zdim: k}), out)
    return out

def first_where(field, condition, zdim):
    """Shallowest valid value satisfying condition."""
    out = xr.full_like(field.isel({zdim: 0}, drop=True), np.nan)
    for k in range(field.sizes[zdim] - 1, -1, -1):
        out = xr.where(condition.isel({zdim: k}), field.isel({zdim: k}), out)
    return out

def main():
    def monthly_bar_prime(field):
        """Calendar-month climatology mapped to time, and monthly anomaly."""
        climatology = field.groupby("time.month").mean("time", skipna=True)
        bar = climatology.sel(month=field.time.dt.month)
        if "month" in bar.coords:
            bar = bar.drop_vars("month")
        bar = bar.assign_coords(time=field.time)
        return bar, field - bar

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    p.add_argument("--output", default="ACCESS_CM2_HBA_VD2001_native_25S25N_Wm2.nc")
    p.add_argument("--nmonths", type=int, default=None, help="number of months; default uses all common months")
    p.add_argument("--chunk-time", type=int, default=24)
    p.add_argument("--kv", type=float, default=1e-5, help="constant E2 diffusivity, m2 s-1")
    args = p.parse_args()

    get = lambda name: open_var(args.data_dir, name, args.chunk_time)
    temp, mld, u, v, wt = get("temp"), get("mld"), get("u"), get("v"), get("wt")
    sw, coupler, pme, frazil = get("sw_heat"), get("sfc_hflux_coupler"), get("sfc_hflux_pme"), get("frazil_3d")
    target = {"lat": temp.lat, "lon": temp.lon}
    u = u.interp(target)
    v = v.interp(target)
    aligned = xr.align(temp, mld, u, v, wt, sw, coupler, pme, frazil, join="inner")
    if not aligned[0].sizes.get("time", 0):
        raise ValueError("No common monthly timestamps across required inputs")
    if args.nmonths:
        aligned = tuple(x.isel(time=slice(0, args.nmonths)) for x in aligned)
    temp, mld, u, v, wt, sw, coupler, pme, frazil = aligned
    # ACCESS-CM2 files are labelled deg_C but contain Kelvin-valued temperatures.
    temp = temp - 273.15
    mld = mld.clip(min=0, max=MAX_MLD)
    print(f"common period: {temp.time.values[0]} to {temp.time.values[-1]} ({temp.sizes['time']} months)")
    zdim = "st_ocean"
    z = temp[zdim]

    # Infer cell interfaces from tracer-cell centres. ACCESS upper levels are
    # 5,15,... m, giving 0,10,... m interfaces.
    edges = np.empty(z.size + 1)
    zv = z.values.astype(float)
    edges[0] = max(0.0, zv[0] - 0.5 * (zv[1] - zv[0]))
    edges[1:-1] = 0.5 * (zv[:-1] + zv[1:])
    edges[-1] = zv[-1] + 0.5 * (zv[-1] - zv[-2])
    top = xr.DataArray(edges[:-1], dims=zdim, coords={zdim: z})
    bottom = xr.DataArray(edges[1:], dims=zdim, coords={zdim: z})
    # Integrate only the fraction of each native layer above the capped MLD.
    layer_base = xr.where(mld < bottom, mld, bottom)
    thick = xr.where(mld.notnull(), (layer_base - top).clip(min=0), np.nan)
    inside = thick > 0
    thick = thick.where(inside)
    h = thick.sum(zdim, skipna=True).where(inside.any(zdim)).rename("H_EXTN")
    tm = ((temp * thick).sum(zdim, skipna=True) / h).rename("T_MIX")

    # Real-time derivatives in SI seconds.
    sec = ((temp.time - temp.time[0]) / np.timedelta64(1, "s")).astype("float64")
    tm_sec = tm.assign_coords(time=sec)
    h_sec = h.assign_coords(time=sec)
    dtmp_dt = tm_sec.differentiate("time").assign_coords(time=temp.time)
    dh_dt = h_sec.differentiate("time").assign_coords(time=temp.time)
    tend = (RHO0 * CP * h * dtmp_dt).rename("TEND")

    # Surface forcing retained by the mixed layer.
    sw_below = sw.sum(zdim, skipna=True) - sw.where(inside).sum(zdim, skipna=True)
    shf = (coupler + pme + frazil.isel({zdim: 0}) - sw_below).rename("SHF")

    # Horizontal gradients on a regular lon/lat analysis grid. Differentiate
    # per degree, then convert degrees to metres locally.
    deg_to_rad = np.pi / 180.0
    dx_per_degree = EARTH_RADIUS * deg_to_rad * np.cos(np.deg2rad(temp.lat))
    dy_per_degree = EARTH_RADIUS * deg_to_rad
    dtdx = temp.differentiate("lon") / dx_per_degree
    dtdy = temp.differentiate("lat") / dy_per_degree
    ax = (-RHO0 * CP * (u * dtdx * thick).sum(zdim, skipna=True)).rename("AX")
    ay = (-RHO0 * CP * (v * dtdy * thick).sum(zdim, skipna=True)).rename("AY")

    # Temperature immediately below MLD and vertical velocity at its base.
    tb = first_where(temp, z >= mld, zdim)
    wdim = "sw_ocean" if "sw_ocean" in wt.dims else zdim
    wz = wt[wdim]
    wb = last_where(wt, wz <= mld, wdim)
    delta = tb - tm
    az = (RHO0 * CP * wb * delta).rename("AZ")
    e1 = (RHO0 * CP * xr.where(dh_dt > 0, dh_dt, 0) * delta).rename("E1")

    # Constant-Kv approximation at the first tracer point below MLD.
    dtdz = temp.differentiate(zdim)
    grad_b = first_where(dtdz, z >= mld, zdim)
    e2 = (RHO0 * CP * args.kv * grad_b).rename("E2")

    # Vialard et al. (2001)-inspired mean/anomaly decomposition. With Omon
    # input, bar is calendar-month climatology and prime is monthly anomaly.
    temp_b, temp_a = monthly_bar_prime(temp)
    u_b, u_a = monthly_bar_prime(u)
    v_b, v_a = monthly_bar_prime(v)
    wt_b, wt_a = monthly_bar_prime(wt)
    tparts = {"B": temp_b, "A": temp_a}
    uparts = {"B": u_b, "A": u_a}
    vparts = {"B": v_b, "A": v_a}
    wtparts = {"B": wt_b, "A": wt_a}
    dtdx_parts = {k: value.differentiate("lon") / dx_per_degree for k, value in tparts.items()}
    dtdy_parts = {k: value.differentiate("lat") / dy_per_degree for k, value in tparts.items()}
    tm_parts = {k: (value * thick).sum(zdim, skipna=True) / h for k, value in tparts.items()}
    tb_parts = {k: first_where(value, z >= mld, zdim) for k, value in tparts.items()}
    delta_parts = {k: tb_parts[k] - tm_parts[k] for k in tparts}
    wb_parts = {k: last_where(value, wz <= mld, wdim) for k, value in wtparts.items()}
    decomp = {}
    # First letter: temperature (A=Tprime, B=Tbar); second: velocity.
    for tkey, vkey in (("B", "B"), ("A", "B"), ("B", "A"), ("A", "A")):
        suffix = tkey + vkey
        decomp[f"AX_{suffix}"] = (-RHO0 * CP * (uparts[vkey] * dtdx_parts[tkey] * thick).sum(zdim, skipna=True)).rename(f"AX_{suffix}")
        decomp[f"AY_{suffix}"] = (-RHO0 * CP * (vparts[vkey] * dtdy_parts[tkey] * thick).sum(zdim, skipna=True)).rename(f"AY_{suffix}")
        decomp[f"AZ_{suffix}"] = (RHO0 * CP * wb_parts[vkey] * delta_parts[tkey]).rename(f"AZ_{suffix}")
        decomp[f"ADV_{suffix}"] = (decomp[f"AX_{suffix}"] + decomp[f"AY_{suffix}"] + decomp[f"AZ_{suffix}"]).rename(f"ADV_{suffix}")

    adv = (ax + ay + az).rename("ADV")
    adv_decomp_sum = sum(decomp[f"ADV_{s}"] for s in ("BB", "AB", "BA", "AA")).rename("ADV_DECOMP_SUM")
    adv_decomp_error = (adv - adv_decomp_sum).rename("ADV_DECOMP_ERROR")

    rhs = (shf + ax + ay + az + e1 + e2).rename("RHS")
    residual = (tend - rhs).rename("RESIDUAL")
    out = xr.Dataset({"T_MIX": tm, "H_EXTN": h, "TEND": tend, "SHF": shf,
                      "AX": ax, "AY": ay, "AZ": az, "ADV": adv,
                      "E1": e1, "E2": e2, "RHS": rhs, "RESIDUAL": residual,
                      "ADV_DECOMP_SUM": adv_decomp_sum,
                      "ADV_DECOMP_ERROR": adv_decomp_error, **decomp})
    out.T_MIX.attrs.update(units="degC", long_name="mixed-layer mean temperature")
    out.H_EXTN.attrs.update(units="m", long_name="discrete mixed-layer thickness")
    for name in ("TEND", "SHF", "AX", "AY", "AZ", "ADV", "E1", "E2", "RHS", "RESIDUAL",
                 "ADV_DECOMP_SUM", "ADV_DECOMP_ERROR", *decomp):
        out[name].attrs["units"] = "W m-2"
    labels = {"BB": "-ubar dot grad(Tbar)", "AB": "-ubar dot grad(Tprime)",
              "BA": "-uprime dot grad(Tbar)", "AA": "-uprime dot grad(Tprime)"}
    for suffix, label in labels.items():
        for term in ("AX", "AY", "AZ", "ADV"):
            out[f"{term}_{suffix}"].attrs.update(
                decomposition=label,
                convention="first letter temperature; second letter velocity")
    out.attrs.update(
        equation="RESIDUAL = TEND - (SHF + AX + AY + AZ + E1 + E2)",
        advection_decomposition="ADV = ADV_BB + ADV_AB + ADV_BA + ADV_AA",
        decomposition_mean="calendar-month climatology over selected calculation period",
        decomposition_warning="Omon input cannot resolve within-month eddy covariance in Vialard et al. (2001)",
        rho0=RHO0, heat_capacity=CP, constant_Kv=args.kv,
        source_grid="native ACCESS-CM2 grid; u/v interpolated to tracer cells",
        region="120-280E, 25S-25N (nearest native cell centres)",
        data_directory=str(args.data_dir.resolve()),
        maximum_mixed_layer_depth_m=MAX_MLD,
        read_depth_m=READ_DEPTH,
        warning="E1 and constant-Kv E2 are reconstructed candidates, not verified one-to-one online diagnostics",
    )
    encoding = {name: {"zlib": True, "complevel": 1} for name in out.data_vars}
    write = out.to_netcdf(args.output, encoding=encoding, compute=False)
    write.compute(scheduler="single-threaded")
    print(f"wrote {args.output}")

if __name__ == "__main__":
    main()
