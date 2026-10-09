"""snowtrack: model-agnostic source tracking of snowmelt through surface and subsurface runoff pathways.

The package partitions bulk storages and fluxes produced by any hydrological or
land-surface model into contributions from individual water sources (rain,
seasonal snowmelt, glacier melt, ...) and reports, for every time step, how much
of each source leaves a model cell through the surface pathway (direct runoff)
and through the subsurface pathway (drainage/baseflow).

Typical use::

    from snowtrack import CellData, SourceTracker, TrackerConfig
    cell = CellData.from_dataframe(df)            # df in the canonical schema
    out = SourceTracker(TrackerConfig()).run(cell)

or from the command line::

    snowtrack run --input-dir my_cells --output-dir tracked
"""

from .schema import CellData, UNKNOWN, SchemaError, canonical_column_summary
from .tracker import SourceTracker, TrackerConfig

__version__ = "2.0.0"
__all__ = [
    "CellData",
    "SchemaError",
    "SourceTracker",
    "TrackerConfig",
    "UNKNOWN",
    "canonical_column_summary",
    "__version__",
]
