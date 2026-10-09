"""Adapter for gridded / HRU netCDF model output (requires ``xarray``).

Configuration (YAML)::

    adapter: netcdf
    engine: h5netcdf             # or netcdf4 / scipy (whatever your xarray can open)
    time_dim: time
    cell_dims: [lat, lon]        # or [hru], [gridcell], ...
    mask_var: RUNOFF_SURF        # cells where this canonical variable is all-NaN are skipped
    variables:                   # canonical name: variable name, or a dict
      S_RAIN: RAIN
      S_SNOW: QSNOMELT
      RUNOFF_SURF: QOVER
      W1_LIQ: {var: SOILLIQ, sum_over: levsoi}
      W1_ICE: {var: SOILICE, sum_over: levsoi}
      Q1_BASE: QDRAI
      E1_TR: {var: QVEGT, scale: 86400}
    scale: 86400                 # optional global multiplier (kg m-2 s-1 -> mm/day)

A dict entry may contain ``var``, ``isel`` (``{dim: index or [indices]}``),
``sum_over`` (dim name or list) and ``scale``.  Per-layer canonical columns
are built by selecting the appropriate ``isel`` indices of a layered variable.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterator, List, Tuple

import numpy as np
import pandas as pd

from .base import Adapter


class NetCDFAdapter(Adapter):
    name = "netcdf"
    patterns = ("*.nc", "*.nc4")

    def __init__(self, config=None) -> None:
        super().__init__(config)
        c = self.config
        if "variables" not in c:
            raise ValueError("netcdf adapter needs a 'variables' section")
        self.variables: Dict[str, Any] = {k.upper(): v for k, v in c["variables"].items()}
        self.time_dim = c.get("time_dim", "time")
        self.cell_dims: List[str] = list(c.get("cell_dims", ["lat", "lon"]))
        self.engine = c.get("engine")
        self.mask_var = str(c.get("mask_var", "RUNOFF_SURF")).upper()
        self.scale = float(c.get("scale", 1.0))

    # ---------------------------------------------------------------- open
    def _open(self, source):
        try:
            import xarray as xr
        except ImportError as exc:  # pragma: no cover
            raise ImportError("The netcdf adapter needs xarray: pip install 'snowtrack[netcdf]'") from exc
        src = Path(source)
        kwargs = {"engine": self.engine} if self.engine else {}
        if src.is_dir():
            files = sorted(p for pat in self.patterns for p in src.glob(pat))
            if not files:
                raise FileNotFoundError(f"no netCDF files in {src}")
            return xr.open_mfdataset(files, combine="by_coords", **kwargs) if len(files) > 1 else xr.open_dataset(files[0], **kwargs)
        return xr.open_dataset(src, **kwargs)

    def _variable(self, ds, spec):
        if isinstance(spec, str):
            da = ds[spec]
            return da * self.scale
        var = ds[spec["var"]]
        if "isel" in spec:
            sel = {k: (v if isinstance(v, int) else list(v)) for k, v in spec["isel"].items()}
            var = var.isel(**sel)
        if "sum_over" in spec:
            dims = spec["sum_over"]
            var = var.sum(dim=dims if isinstance(dims, str) else list(dims), skipna=False)
        return var * float(spec.get("scale", 1.0)) * self.scale

    # ---------------------------------------------------------- iteration
    def iter_items(self, source) -> Iterator[Tuple[str, Any]]:
        ds = self._open(source)
        arrays = {canon: self._variable(ds, spec) for canon, spec in self.variables.items()}
        time = pd.DatetimeIndex(pd.to_datetime(ds[self.time_dim].values))
        mask_da = arrays.get(self.mask_var)
        sizes = [ds.sizes[d] for d in self.cell_dims]
        for idx in np.ndindex(*sizes):
            sel = {d: int(i) for d, i in zip(self.cell_dims, idx)}
            if mask_da is not None:
                col = mask_da.isel(**sel).values
                if np.all(np.isnan(col)):
                    continue
            coords = [float(ds[d].values[i]) if d in ds.coords else int(i) for d, i in zip(self.cell_dims, idx)]
            name = "cell_" + "_".join(f"{c:.4f}" if isinstance(c, float) else str(c) for c in coords)
            table = {"DATE": time}
            for canon, da in arrays.items():
                table[canon] = np.asarray(da.isel(**sel).values, dtype=float)
            yield name, pd.DataFrame(table)

    def load(self, payload: Any) -> pd.DataFrame:
        if isinstance(payload, pd.DataFrame):
            return payload
        raise TypeError("netcdf adapter payloads are DataFrames produced by iter_items")
