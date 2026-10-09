"""Route tracked cell components to a basin outlet.

Default scheme: **Lohmann channel routing** (Lohmann et al., 1996, 1998), the
routing model used with VIC.  The runoff of every cell is convolved with the
impulse response of the linearised Saint-Venant equation (diffusion wave)
over the flow distance from the cell to the gauge,

    h(x, t) = x / (2 sqrt(pi D t^3)) * exp(-(x - C t)^2 / (4 D t)),

with wave velocity ``C`` and diffusivity ``D``.  Surface components
(``RUNOFF_<SRC>_SURFACE``, ``RUNOFF_SURF``) are routed with the surface
velocity, subsurface components (``RUNOFF_<SRC>_SUBSURFACE``, ``BASEFLOW``)
with the (slower) subsurface velocity; totals are the sums of the routed
components.  ``h`` is the inverse-Gaussian density, so the mass arriving in
every time step is evaluated exactly from its cumulative distribution.

**The velocities (and the diffusivity) are calibration parameters.**  They
are calibrated against gauged discharge together with the hydrological model
and must be given explicitly; there are no defaults.

Flow distance: by default the great-circle distance from the cell centre
(parsed from the file name ``..._<lat>_<lon>``) to the outlet; a cells table
(``cell_id, distance_m[, weight|area]``) can supply flow-path distances from
a river network instead, or, with ``lag_surface, lag_base`` columns (whole
steps), prescribe the translation of every cell directly (``scheme="lag"``).
Output columns carry the suffix ``_OUTLET``.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

_LATLON_RE = re.compile(r"(-?\d+(?:\.\d+)?)_(-?\d+(?:\.\d+)?)")
SCHEMES = ("lohmann", "lag")

LEGACY_RENAMES = {
    "RUNOFF_SNOW": "RUNOFF_SNOW_TOTAL",
    "RUNOFF_SNOW_BASE": "RUNOFF_SNOW_SUBSURFACE",
}


@dataclass
class RoutingConfig:
    """Routing parameters.  ``surface_velocity`` and ``base_velocity`` (m/s) are
    calibration parameters and must be set; ``diffusivity`` (m^2/s) is the
    diffusion-wave diffusivity of the Lohmann scheme (also calibratable)."""

    surface_velocity: Optional[float] = None
    base_velocity: Optional[float] = None
    diffusivity: float = 800.0
    scheme: str = "lohmann"
    weight: str = "none"                       # or "area"
    grid_resolution_deg: Optional[float] = None
    max_kernel_steps: int = 3650
    kernel_mass_tolerance: float = 1.0e-4

    def __post_init__(self) -> None:
        if self.scheme not in SCHEMES:
            raise ValueError(f"scheme must be one of {SCHEMES}")

    def require_velocities(self) -> None:
        for name in ("surface_velocity", "base_velocity"):
            v = getattr(self, name)
            if v is None or not (v > 0):
                raise ValueError(
                    f"{name} is a calibration parameter and must be given explicitly (m/s > 0); "
                    "calibrate it against gauged discharge together with the hydrological model."
                )
        if not self.diffusivity > 0:
            raise ValueError("diffusivity must be > 0 (m^2/s)")


# ------------------------------------------------------------------ geometry
def parse_latlon(name: str) -> Optional[Tuple[float, float]]:
    m = _LATLON_RE.search(str(name))
    if not m:
        return None
    return float(m.group(1)), float(m.group(2))


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    lon1, lat1, lon2, lat2 = map(math.radians, (lon1, lat1, lon2, lat2))
    a = math.sin((lat2 - lat1) / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin((lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * math.asin(math.sqrt(a))


def grid_area_km2(lat: float, res_deg: float) -> float:
    R = 6371.0
    lat1 = math.radians(lat - res_deg / 2)
    lat2 = math.radians(lat + res_deg / 2)
    return R * R * abs(math.sin(lat2) - math.sin(lat1)) * math.radians(res_deg)


# ------------------------------------------------------------------- kernels
def _log_norm_cdf(z: float) -> float:
    """log Phi(z), stable for very negative z."""
    if z > -30.0:
        return math.log(0.5 * math.erfc(-z / math.sqrt(2.0)))
    return -0.5 * z * z - math.log(-z) - 0.5 * math.log(2.0 * math.pi) + math.log(1.0 - 1.0 / z ** 2 + 3.0 / z ** 4)


def inverse_gaussian_cdf(t: float, mu: float, lam: float) -> float:
    """CDF of the inverse-Gaussian (first-passage-time) distribution, mean ``mu``, shape ``lam``."""
    if t <= 0.0:
        return 0.0
    a = math.sqrt(lam / t)
    term1 = 0.5 * math.erfc(-(a * (t / mu - 1.0)) / math.sqrt(2.0))
    term2 = math.exp(2.0 * lam / mu + _log_norm_cdf(-a * (t / mu + 1.0)))
    return min(term1 + term2, 1.0)


def lohmann_kernel(distance_m: float, velocity_m_s: float, diffusivity_m2_s: float,
                   step_seconds: float = 86400.0, max_steps: int = 3650, mass_tol: float = 1.0e-4) -> np.ndarray:
    """Fraction of a cell's runoff of one step arriving at the outlet 0, 1, 2, ... steps later.

    Diffusion-wave impulse response of Lohmann et al. (1996) for flow distance
    ``distance_m``: an inverse-Gaussian distribution with mean ``x / C`` and
    shape ``x^2 / (2 D)``, integrated exactly over every time step.
    """
    if distance_m <= 0.0:
        return np.array([1.0])
    mu = distance_m / velocity_m_s
    lam = distance_m ** 2 / (2.0 * diffusivity_m2_s)
    weights: List[float] = []
    prev = 0.0
    for k in range(1, max_steps + 1):
        cur = inverse_gaussian_cdf(k * step_seconds, mu, lam)
        weights.append(cur - prev)
        prev = cur
        if 1.0 - cur <= mass_tol and k * step_seconds > mu:
            break
    w = np.asarray(weights, dtype=float)
    w[w < 0.0] = 0.0
    s = w.sum()
    return w / s if s > 0 else np.array([1.0])


def lag_steps(distance_km: float, velocity_m_s: float, step_seconds: float = 86400.0) -> int:
    days = distance_km * 1000.0 / (velocity_m_s * step_seconds)
    return 0 if days <= 0 else int(round(days))


def lag_kernel(steps: int) -> np.ndarray:
    k = np.zeros(max(int(steps), 0) + 1)
    k[-1] = 1.0
    return k


def convolve(series: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    x = np.where(np.isfinite(series), series, 0.0)
    if len(kernel) == 1 and kernel[0] == 1.0:
        return x.copy()
    return np.convolve(x, kernel)[: len(x)]


# ------------------------------------------------------------------- routing
def read_tracked(path) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["DATE"], compression="infer")
    df = df.rename(columns=LEGACY_RENAMES)
    return df.set_index("DATE").sort_index()


def _component_columns(columns: Iterable[str]) -> Tuple[List[str], List[str]]:
    surface, subsurface = [], []
    for c in columns:
        if c.endswith("_SURFACE") and c.startswith("RUNOFF_"):
            surface.append(c)
        elif c.endswith("_SUBSURFACE") and c.startswith("RUNOFF_"):
            subsurface.append(c)
    if "RUNOFF_SURF" in columns:
        surface.append("RUNOFF_SURF")
    if "BASEFLOW" in columns:
        subsurface.append("BASEFLOW")
    return surface, subsurface


def _strip(name: str) -> str:
    for suf in (".gz", ".csv"):
        if name.endswith(suf):
            name = name[: -len(suf)]
    if name.startswith("tracked_"):
        name = name[len("tracked_"):]
    return name


def cell_kernels(cell_id: str, outlet: Optional[Tuple[float, float]], cfg: RoutingConfig,
                 table: Optional[pd.DataFrame], step_seconds: float) -> Tuple[np.ndarray, np.ndarray, float, Optional[Tuple[float, float]]]:
    """Surface kernel, subsurface kernel, weight and lat/lon of one cell."""
    ll = parse_latlon(cell_id)
    row = table.loc[cell_id] if (table is not None and cell_id in table.index) else None
    weight = 1.0
    if row is not None:
        if "weight" in table.columns and pd.notna(row["weight"]):
            weight = float(row["weight"])
        elif "area" in table.columns and pd.notna(row["area"]):
            weight = float(row["area"])
    elif cfg.weight == "area" and cfg.grid_resolution_deg and ll is not None:
        weight = grid_area_km2(ll[0], cfg.grid_resolution_deg)

    if row is not None and {"lag_surface", "lag_base"} <= set(table.columns) and pd.notna(row["lag_surface"]):
        return lag_kernel(int(row["lag_surface"])), lag_kernel(int(row["lag_base"])), weight, ll

    if row is not None and "distance_m" in table.columns and pd.notna(row["distance_m"]):
        distance_m = float(row["distance_m"])
    else:
        if outlet is None or ll is None:
            raise ValueError(f"{cell_id}: need an outlet and lat/lon in the file name, or a cells table with distance_m / lags")
        distance_m = haversine_km(outlet[0], outlet[1], ll[0], ll[1]) * 1000.0

    cfg.require_velocities()
    if cfg.scheme == "lag":
        ks = lag_kernel(lag_steps(distance_m / 1000.0, cfg.surface_velocity, step_seconds))
        kb = lag_kernel(lag_steps(distance_m / 1000.0, cfg.base_velocity, step_seconds))
    else:
        ks = lohmann_kernel(distance_m, cfg.surface_velocity, cfg.diffusivity, step_seconds, cfg.max_kernel_steps, cfg.kernel_mass_tolerance)
        kb = lohmann_kernel(distance_m, cfg.base_velocity, cfg.diffusivity, step_seconds, cfg.max_kernel_steps, cfg.kernel_mass_tolerance)
    return ks, kb, weight, ll


def route(
    cells: Union[str, Path, Sequence[Tuple[str, Union[str, Path, pd.DataFrame]]]],
    outlet: Optional[Tuple[float, float]] = None,
    config: Optional[RoutingConfig] = None,
    cells_table: Optional[pd.DataFrame] = None,
    pattern: str = "tracked_*",
    step_seconds: float = 86400.0,
    log=print,
) -> pd.DataFrame:
    """Route all cells to the outlet; returns a DataFrame indexed by DATE (columns ``*_OUTLET``)."""
    cfg = config or RoutingConfig()
    if isinstance(cells, (str, Path)):
        paths = sorted(p for p in Path(cells).glob(pattern) if p.is_file())
        items: List[Tuple[str, Union[Path, pd.DataFrame]]] = [(p.name, p) for p in paths]
    else:
        items = list(cells)
    if not items:
        raise FileNotFoundError("no tracked cell files found")

    table = None
    if cells_table is not None:
        table = cells_table.copy()
        table.columns = [c.lower() for c in table.columns]
        table["cell_id"] = table["cell_id"].astype(str)
        table = table.set_index("cell_id")

    acc: Dict[str, pd.Series] = {}
    n_used = 0
    weights_sum = 0.0
    for name, payload in items:
        cell_id = _strip(name)
        try:
            ks, kb, w, _ = cell_kernels(cell_id, outlet, cfg, table, step_seconds)
        except ValueError as exc:
            if "calibration parameter" in str(exc):
                raise
            log(f"  skipping {name}: {exc}")
            continue
        if cfg.weight == "none" and (table is None or cell_id not in table.index):
            w = 1.0
        df = payload if isinstance(payload, pd.DataFrame) else read_tracked(payload)
        surf_cols, sub_cols = _component_columns(df.columns)
        for cols, kernel in ((surf_cols, ks), (sub_cols, kb)):
            for c in cols:
                s = pd.Series(convolve(df[c].to_numpy(float), kernel) * w, index=df.index)
                acc[c] = s if c not in acc else acc[c].add(s, fill_value=0.0)
        n_used += 1
        weights_sum += w
    if n_used == 0:
        raise ValueError("no cell could be routed")

    out = pd.DataFrame(acc).sort_index()
    if cfg.weight == "area" and weights_sum > 0:
        out = out / weights_sum            # basin-average depth
    sources = sorted({c[len("RUNOFF_"):-len("_SURFACE")] for c in out.columns if c.startswith("RUNOFF_") and c.endswith("_SURFACE")})
    for s in sources:
        out[f"RUNOFF_{s}_TOTAL"] = out.get(f"RUNOFF_{s}_SURFACE", 0.0) + out.get(f"RUNOFF_{s}_SUBSURFACE", 0.0)
    if "RUNOFF_SURF" in out.columns or "BASEFLOW" in out.columns:
        out["TOTAL_RUNOFF"] = out.get("RUNOFF_SURF", 0.0) + out.get("BASEFLOW", 0.0)
    out.columns = [f"{c}_OUTLET" for c in out.columns]
    out.index.name = "DATE"
    out.attrs["n_cells"] = n_used
    return out


def write_outlet(df: pd.DataFrame, path, outlet: Optional[Tuple[float, float]] = None, header_note: str = "") -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8", newline="") as f:
        if outlet is not None:
            f.write(f"# Routed runoff components at outlet ({outlet[0]:.4f}, {outlet[1]:.4f}); cells={df.attrs.get('n_cells', '?')}\n")
        if header_note:
            f.write(f"# {header_note}\n")
        f.write("# Generated by snowtrack route\n")
        df.reset_index().to_csv(f, index=False, float_format="%.6f")
    return p
