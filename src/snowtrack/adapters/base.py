"""Adapter base class and shared file helpers."""

from __future__ import annotations

import glob
import os
from pathlib import Path
from typing import Any, Iterator, List, Optional, Sequence, Tuple

import pandas as pd

DEFAULT_PATTERNS: Tuple[str, ...] = ("*.csv", "*.csv.gz", "*.txt", "*.txt.gz", "*.dat", "*.dat.gz")
SKIP_PREFIXES: Tuple[str, ...] = ("tracked_", "annual_", "monthly_", "overall_", "domain_", "routed_")


def list_files(source: os.PathLike | str, patterns: Sequence[str] = DEFAULT_PATTERNS,
               skip_prefixes: Sequence[str] = SKIP_PREFIXES) -> List[Path]:
    """Files in ``source`` matching any pattern (a single file is returned as is)."""
    src = Path(source)
    if src.is_file():
        return [src]
    if not src.is_dir():
        raise FileNotFoundError(f"Input not found: {src}")
    found = set()
    for pat in patterns:
        found.update(Path(p) for p in glob.glob(str(src / pat)))
    files = sorted(p for p in found if p.is_file() and not p.name.startswith(tuple(skip_prefixes)))
    return files


def cell_name_from_path(path: os.PathLike | str) -> str:
    """``fluxes_31.0500_97.1500`` / ``grid_001.csv.gz`` -> ``fluxes_31.0500_97.1500`` / ``grid_001``."""
    name = Path(path).name
    for suffix in (".gz", ".csv", ".txt", ".dat"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


class Adapter:
    """Base class.  Subclasses set ``name`` and implement :meth:`load`."""

    name = "base"
    patterns: Tuple[str, ...] = DEFAULT_PATTERNS

    def __init__(self, config: Optional[dict] = None) -> None:
        self.config = dict(config or {})
        pats = self.config.get("pattern")
        if pats:
            self.patterns = tuple([pats] if isinstance(pats, str) else pats)

    # ---- iteration over cells --------------------------------------------
    def iter_items(self, source: os.PathLike | str) -> Iterator[Tuple[str, Any]]:
        """Yield ``(cell_name, payload)``; the default payload is a file path."""
        for path in list_files(source, self.patterns):
            yield cell_name_from_path(path), path

    def load(self, payload: Any) -> pd.DataFrame:
        """Turn a payload into a canonical DataFrame."""
        raise NotImplementedError

    # ---- convenience -------------------------------------------------------
    def read_cell(self, path: os.PathLike | str) -> pd.DataFrame:
        return self.load(Path(path))

    def describe(self) -> str:
        return f"{self.name} adapter"
