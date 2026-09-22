#!/usr/bin/env python3
"""Validate (and optionally link) ACCESS-CM2 inputs for the heat budget."""
from __future__ import annotations
import argparse
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODEL_TOKEN = "bj594_piControl"
from netCDF4 import Dataset

REQUIRED = {
    "temp": "potential temperature", "mld": "mixed-layer depth", "u": "zonal velocity",
    "v": "meridional velocity", "wt": "vertical velocity", "sw_heat": "3-D shortwave heating",
    "sfc_hflux_coupler": "coupled surface heat flux", "sfc_hflux_pme": "PME heat flux",
    "frazil_3d": "frazil heating", "temp_vdiffuse_sbc": "surface-boundary heat diagnostic",
}
OPTIONAL = {
    "temp_tendency": "model tendency reference",
    "temp_xflux_adv": "native x advective heat transport",
    "temp_yflux_adv": "native y advective heat transport",
    "temp_zflux_adv": "native z advective heat transport",
    "temp_vdiffuse_diff_cbt": "vertical mixing candidate for E2",
    "temp_nonlocal_KPP": "non-local KPP candidate for E2",
    "temp_vdiffuse_k33": "K33 vertical mixing candidate for E2",
    "temp_sigma_diff": "sigma diffusion candidate for E2",
}

def find(root: Path, stem: str):
    hits = sorted(root.glob(f"{stem}_Omon_*{MODEL_TOKEN}*.nc"))
    return hits[0] if hits else None

def inspect(path: Path, variable: str):
    with Dataset(path) as nc:
        if variable not in nc.variables:
            return False, "variable absent"
        v = nc.variables[variable]
        return True, f"shape={v.shape}, units={getattr(v, 'units', '?')}"

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, default=PROJECT_ROOT / "data")
    p.add_argument("--link-dir", type=Path, help="create symlinks here; omit for validation only")
    args = p.parse_args()
    missing = []
    for group, needed in (("REQUIRED", REQUIRED), ("OPTIONAL", OPTIONAL)):
        print(f"\n{group}")
        for stem, purpose in needed.items():
            path = find(args.source, stem)
            if path is None:
                print(f"  MISSING {stem}: {purpose}")
                if group == "REQUIRED": missing.append(stem)
                continue
            ok, meta = inspect(path, stem)
            print(f"  {'OK' if ok else 'BAD'} {stem}: {path.name}; {meta}")
            if group == "REQUIRED" and not ok: missing.append(stem)
            if args.link_dir and ok:
                args.link_dir.mkdir(parents=True, exist_ok=True)
                link = args.link_dir / path.name
                if not link.exists(): link.symlink_to(path.resolve())
    if missing:
        raise SystemExit("Required inputs missing/bad: " + ", ".join(missing))

if __name__ == "__main__":
    main()
