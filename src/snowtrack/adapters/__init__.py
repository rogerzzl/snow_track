"""Adapters turn model output into the canonical per-cell schema.

An adapter yields ``(cell_name, payload)`` items for a source (a directory of
per-cell files, a single file, or a netCDF dataset) and knows how to turn a
payload into a canonical ``pandas.DataFrame`` (:meth:`Adapter.load`).

Built-in adapters
-----------------
``canonical``    files already in the canonical schema (CSV / whitespace, optionally gzipped)
``vic-classic``  VIC 4.x classic-driver ASCII ``fluxes_*`` files (as used in Zhu et al.)
``mapping``      any tabular per-cell output, with a YAML column mapping (HBV, SWAT, GR4J, ...)
``netcdf``       gridded or HRU netCDF output (CLM, Noah-MP, VIC image driver, ...)

Register your own with :func:`register_adapter`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Type

import yaml

from .base import Adapter
from .canonical import CanonicalAdapter
from .mapping import MappingAdapter
from .netcdf import NetCDFAdapter
from .vic_classic import VicClassicAdapter

ADAPTERS: Dict[str, Type[Adapter]] = {
    CanonicalAdapter.name: CanonicalAdapter,
    VicClassicAdapter.name: VicClassicAdapter,
    MappingAdapter.name: MappingAdapter,
    NetCDFAdapter.name: NetCDFAdapter,
}


def register_adapter(cls: Type[Adapter]) -> Type[Adapter]:
    ADAPTERS[cls.name] = cls
    return cls


def load_config(path: Optional[str]) -> dict:
    """Read a YAML (or JSON) adapter configuration file; ``None`` gives ``{}``."""
    if path is None:
        return {}
    text = Path(path).read_text(encoding="utf-8")
    data = yaml.safe_load(text) or {}
    if not isinstance(data, dict):
        raise ValueError(f"{path}: configuration must be a mapping at top level")
    return data


def get_adapter(name: Optional[str] = None, config: Optional[dict] = None) -> Adapter:
    """Instantiate an adapter by name (or by the ``adapter`` key of ``config``)."""
    config = dict(config or {})
    name = name or config.get("adapter") or "canonical"
    name = name.replace("_", "-").lower()
    if name not in ADAPTERS:
        raise KeyError(f"Unknown adapter {name!r}. Available: {sorted(ADAPTERS)}")
    return ADAPTERS[name](config)


__all__ = [
    "Adapter", "CanonicalAdapter", "VicClassicAdapter", "MappingAdapter", "NetCDFAdapter",
    "ADAPTERS", "register_adapter", "get_adapter", "load_config",
]
