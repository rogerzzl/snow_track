"""Adapter for files that are already in the canonical schema."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .base import Adapter


def read_table(path: Path, sep: str | None = None, comment: str = "#", **kwargs) -> pd.DataFrame:
    """Read CSV or whitespace-delimited text (optionally ``.gz``) with a header row."""
    if sep is None:
        sep = r"\s+|,"
        engine = "python"
    else:
        engine = "python" if len(sep) > 1 else None
    return pd.read_csv(path, sep=sep, engine=engine, comment=comment, compression="infer", **kwargs)


class CanonicalAdapter(Adapter):
    name = "canonical"

    def load(self, payload: Any) -> pd.DataFrame:
        path = Path(payload)
        df = read_table(path, sep=self.config.get("sep"), comment=self.config.get("comment", "#"))
        df.columns = [str(c).strip().upper() for c in df.columns]
        return df
