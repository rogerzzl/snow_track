"""Generic column-mapping adapter for any tabular per-cell model output.

Configuration (YAML)::

    adapter: mapping
    read:                      # passed to pandas.read_csv (all optional)
      sep: ","                 # "\\s+" for whitespace
      header: 0                # null + names: [...] for header-less files
      comment: "#"
      skiprows: 0
    time:
      date: date               # a single date column ...
      # year: yr; month: mo; day: dy     # ... or separate columns
    columns:                   # canonical name: raw column name OR pandas expression
      S_RAIN: rain
      S_SNOW: melt
      RUNOFF_SURF: qsurf
      INFIL: infil             # optional
      W1_LIQ: "sm_upper + sm_lower"
      Q1_BASE: qbase
      E1_TR: et_soil
    scale:                     # optional multipliers after mapping (unit conversion)
      S_RAIN: 1.0
    keep: [Q_TOTAL]            # optional raw columns passed through to the output

Expressions are evaluated with ``DataFrame.eval``; wrap odd column names in
backticks (``` `soil moisture` ```).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import pandas as pd

from .base import Adapter


class MappingAdapter(Adapter):
    name = "mapping"

    def __init__(self, config=None) -> None:
        super().__init__(config)
        if "columns" not in self.config:
            raise ValueError("mapping adapter needs a 'columns' section (canonical name -> raw column/expression)")
        self.columns: Dict[str, str] = {k.upper(): str(v) for k, v in self.config["columns"].items()}
        self.read_kwargs = dict(self.config.get("read", {}))
        self.time = {k.lower(): v for k, v in self.config.get("time", {"date": "DATE"}).items()}
        self.scale = {k.upper(): float(v) for k, v in self.config.get("scale", {}).items()}
        self.keep = list(self.config.get("keep", []))

    def read_raw(self, path: Path) -> pd.DataFrame:
        kwargs = dict(self.read_kwargs)
        sep = kwargs.pop("sep", ",")
        if len(sep) > 1 or sep.startswith("\\"):
            kwargs.setdefault("engine", "python")
        kwargs.setdefault("comment", "#")
        return pd.read_csv(path, sep=sep, compression="infer", **kwargs)

    def load(self, payload: Any) -> pd.DataFrame:
        raw = self.read_raw(Path(payload)) if not isinstance(payload, pd.DataFrame) else payload
        return self.to_canonical(raw)

    def to_canonical(self, raw: pd.DataFrame) -> pd.DataFrame:
        out = pd.DataFrame(index=raw.index)
        if "date" in self.time:
            out["DATE"] = pd.to_datetime(raw[self.time["date"]])
        else:
            for key in ("year", "month", "day"):
                if key not in self.time:
                    raise ValueError("time section needs 'date' or 'year'/'month'/'day'")
                out[key.upper()] = raw[self.time[key]]
            if "hour" in self.time:
                out["HOUR"] = raw[self.time["hour"]]
        for canon, expr in self.columns.items():
            if expr in raw.columns:
                values = raw[expr]
            else:
                values = raw.eval(expr)
            out[canon] = pd.to_numeric(values, errors="coerce") * self.scale.get(canon, 1.0)
        for k in self.keep:
            if k in raw.columns:
                out[k.upper()] = raw[k]
        return out.reset_index(drop=True)
