"""Multi-layer, dual-phase source tracker.

The tracker follows the scheme described in ``docs/method.md`` (Supplementary
Text S3 of Zhu et al.): source-resolved liquid and ice storages are carried
for every soil layer and every source; all bulk fluxes delivered by the host
model are partitioned among sources in proportion to the composition of the
store they leave, and the tracked storages are reconciled to the model's bulk
storages at the end of every step so that the decomposition is exactly
mass-consistent with the host model.

The number of layers (N >= 1), the presence of an ice phase and the set of
sources are all taken from the input (:class:`snowtrack.schema.CellData`).
With N = 1 and no ice columns the scheme reduces to a single well-mixed soil
store, which is the configuration that can be driven by most model outputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from .schema import UNKNOWN, CellData

COMPOSITIONS = ("end", "start", "legacy")


@dataclass(frozen=True)
class TrackerConfig:
    """Numerical and bookkeeping options.

    composition
        Which composition carries the withdrawals (percolation, baseflow,
        transpiration) that leave a layer during a step.

        ``"end"``    (default) inflow is mixed into the store first, then all
                     withdrawals carry the mixed composition (implicit,
                     well-mixed reservoir; inflow of a step can leave the
                     same step).
        ``"start"``  withdrawals carry the composition at the start of the
                     step (explicit; inflow of a step is only available from
                     the next step on).
        ``"legacy"`` withdrawals carry the start-of-step composition but the
                     *reported* baseflow carries the end-of-step composition
                     after reconciliation.  This is the convention of the code
                     that produced the archived results of Zhu et al. and is
                     what ``configs/vic_classic_zhu2026.yaml`` selects.
    epsilon
        Storage threshold below which a store is treated as empty (mm).
    strict
        Raise instead of clipping when a tracked storage would go more
        negative than ``negative_tolerance`` (a sign that the exported fluxes
        are inconsistent with the exported storages).
    strict_surface_closure
        Raise when ``sum(S_*) != INFIL + RUNOFF_SURF`` beyond
        ``surface_tolerance``; otherwise the unexplained part of infiltration
        or direct runoff is attributed to the UNKNOWN source.
    reconcile
        Rescale tracked storages to the bulk storages every step (recommended;
        switch off only for synthetic, perfectly closed inputs).
    spinup_years
        Years dropped from the *output* (the tracker still runs through them so
        that the UNKNOWN initial water is flushed).
    store_storage_columns
        Write the source-resolved storages (``W<n>_<SOURCE>_LIQ/ICE``) to the
        output.  Switch off for compact outputs.
    """

    composition: str = "end"
    epsilon: float = 1.0e-8
    strict: bool = False
    negative_tolerance: float = 1.0e-6
    strict_surface_closure: bool = False
    surface_tolerance: float = 1.0e-6
    reconcile: bool = True
    spinup_years: int = 0
    store_storage_columns: bool = True

    def __post_init__(self) -> None:
        if self.composition not in COMPOSITIONS:
            raise ValueError(f"composition must be one of {COMPOSITIONS}, got {self.composition!r}")


def _composition(store: np.ndarray, eps: float) -> np.ndarray:
    total = store.sum()
    if total < eps:
        return np.zeros_like(store)
    return store / total


class SourceTracker:
    """Run the source-tracking scheme for one cell (:meth:`run`)."""

    def __init__(self, config: TrackerConfig | None = None) -> None:
        self.config = config or TrackerConfig()

    # ------------------------------------------------------------------ main
    def run(self, cell: CellData) -> pd.DataFrame:
        cfg = self.config
        eps = cfg.epsilon
        T, N = cell.w_liq.shape
        sources: Tuple[str, ...] = tuple(cell.sources) + (UNKNOWN,)
        K = len(sources)          # active sources + UNKNOWN
        KA = K - 1                # active sources
        iu = K - 1                # index of UNKNOWN

        liq = np.zeros((N, K))
        ice = np.zeros((N, K))
        liq[:, iu] = np.maximum(cell.w_liq[0], 0.0)
        ice[:, iu] = np.maximum(cell.w_ice[0], 0.0)

        # output buffers
        surf = np.zeros((T, K))
        base = np.zeros((T, K))
        fw = np.zeros((T, K))                 # share of total column liquid storage
        liq_out = np.zeros((T, N, K)) if cfg.store_storage_columns else None
        ice_out = np.zeros((T, N, K)) if (cfg.store_storage_columns and cell.has_ice) else None
        resid_liq = np.zeros((T, N))
        resid_ice = np.zeros((T, N))
        clipped = np.zeros(T)
        injected = np.zeros(T)

        unit_unknown = np.zeros(K)
        unit_unknown[iu] = 1.0

        legacy = cfg.composition == "legacy"
        mix_first = cfg.composition == "end"

        for t in range(T):
            if t == 0:
                # no previous state: storages at t0 belong to UNKNOWN and so does any runoff of the first step
                surf[t, iu] = cell.runoff_surf[t]
                base[t, iu] = cell.q_base[t].sum()
                self._record_storage(t, liq, ice, fw, liq_out, ice_out, eps)
                continue

            # ---------------------------------------------------- surface
            s_tot = cell.supply_total[t]
            infil = cell.infil[t]
            r_surf = cell.runoff_surf[t]
            frac = np.zeros(K)
            if s_tot > eps:
                frac[:KA] = np.maximum(cell.supply[t], 0.0) / s_tot
            if cfg.strict_surface_closure and abs(s_tot - infil - r_surf) > cfg.surface_tolerance:
                raise ValueError(
                    f"{cell.name}: surface closure violated on {cell.dates[t].date()}: "
                    f"sum(S)={s_tot:.6f}, INFIL+RUNOFF_SURF={infil + r_surf:.6f}"
                )
            # attributable parts; the remainder is water the supply cannot explain
            attrib_surf = min(r_surf, max(s_tot - min(infil, s_tot), 0.0)) if s_tot > eps else 0.0
            attrib_infil = min(infil, s_tot) if s_tot > eps else 0.0
            surf[t] = attrib_surf * frac
            surf[t, iu] += r_surf - attrib_surf
            inflow = attrib_infil * frac
            inflow[iu] += infil - attrib_infil

            # -------------------------------------------- transport + ET
            base_step = np.zeros(K)
            base_layer_end: List[np.ndarray] = []
            for n in range(N):
                store = liq[n]
                q_perc = cell.q_perc[t, n]
                q_base = cell.q_base[t, n]
                e_tr = cell.e_tr[t, n]
                withdraw = q_perc + q_base + e_tr
                if mix_first:
                    mixed = store + inflow
                    c = _composition(mixed, eps)
                    if withdraw > 0 and c.sum() == 0.0:
                        c = unit_unknown
                    new = mixed - withdraw * c
                else:
                    c = _composition(store, eps)
                    if withdraw > 0 and c.sum() == 0.0:
                        c = unit_unknown
                    new = store + inflow - withdraw * c
                clipped[t] += self._clip_negative(new, cfg, cell, t, "liquid", n)
                liq[n] = new
                base_k = q_base * c
                base_step += base_k
                base_layer_end.append((n, q_base))
                inflow = q_perc * c          # percolation feeds the next layer

            # ----------------------------------------------- phase change
            if cell.has_ice:
                for n in range(N):
                    c_liq = _composition(liq[n], eps)
                    c_ice = _composition(ice[n], eps)
                    fr = cell.phi_freeze[t, n] if liq[n].sum() >= eps else 0.0
                    me = cell.phi_melt[t, n] if ice[n].sum() >= eps else 0.0
                    freeze = fr * c_liq
                    melt = me * c_ice
                    liq[n] = liq[n] - freeze + melt
                    ice[n] = ice[n] + freeze - melt
                    clipped[t] += self._clip_negative(liq[n], cfg, cell, t, "liquid", n)
                    clipped[t] += self._clip_negative(ice[n], cfg, cell, t, "ice", n)

            # ------------------------------------------------- reconcile
            if cfg.reconcile:
                for n in range(N):
                    resid_liq[t, n] = cell.w_liq[t, n] - liq[n].sum()
                    injected[t] += self._reconcile(liq[n], cell.w_liq[t, n], eps, iu, legacy)
                    if cell.has_ice:
                        resid_ice[t, n] = cell.w_ice[t, n] - ice[n].sum()
                        injected[t] += self._reconcile(ice[n], cell.w_ice[t, n], eps, iu, legacy)

            if legacy:
                # reported baseflow carries the end-of-step composition
                base_step = np.zeros(K)
                for n, q_base in base_layer_end:
                    c_end = _composition(liq[n], eps)
                    if q_base > 0 and c_end.sum() == 0.0:
                        c_end = unit_unknown
                    base_step += q_base * c_end
            base[t] = base_step

            self._record_storage(t, liq, ice, fw, liq_out, ice_out, eps)

        return self._assemble(cell, sources, surf, base, fw, liq_out, ice_out, resid_liq, resid_ice, clipped, injected)

    # --------------------------------------------------------------- helpers
    @staticmethod
    def _record_storage(t, liq, ice, fw, liq_out, ice_out, eps) -> None:
        column = liq.sum(axis=0)
        fw[t] = _composition(column, eps)
        if liq_out is not None:
            liq_out[t] = liq
        if ice_out is not None:
            ice_out[t] = ice

    @staticmethod
    def _clip_negative(values: np.ndarray, cfg: TrackerConfig, cell: CellData, t: int, phase: str, n: int) -> float:
        neg = values < 0.0
        if not neg.any():
            return 0.0
        worst = float(values[neg].min())
        if cfg.strict and worst < -cfg.negative_tolerance:
            raise ValueError(
                f"{cell.name}: source-resolved {phase} storage of layer {n + 1} became negative "
                f"({worst:.3e} mm) on {cell.dates[t].date()}. The exported fluxes withdraw more water "
                "than the layer holds; check the adapter mapping or run without --strict."
            )
        lost = float(-values[neg].sum())
        values[neg] = 0.0
        return lost

    @staticmethod
    def _reconcile(store: np.ndarray, target: float, eps: float, iu: int, legacy: bool) -> float:
        """Rescale ``store`` (in place) so that it sums to the bulk ``target``.

        Returns the amount of water newly assigned to UNKNOWN (0 unless the
        tracked store was empty while the bulk store is not).
        """
        target = max(float(target), 0.0)
        current = float(store.sum())
        if target < eps:
            store[:] = 0.0
            return 0.0
        if legacy:
            # exact transcription of the archived implementation:
            # fractions = clip(store / W_bulk, 0, 1), renormalised to 1
            f = np.clip(store / target, 0.0, 1.0)
            s = f.sum()
            if s > 0.0:
                f /= s
            store[:] = f * target
            return 0.0
        if current < eps:
            store[:] = 0.0
            store[iu] = target
            return target
        store *= target / current
        return 0.0

    def _assemble(self, cell, sources, surf, base, fw, liq_out, ice_out, resid_liq, resid_ice, clipped, injected) -> pd.DataFrame:
        cfg = self.config
        T, N = cell.w_liq.shape
        dates = cell.dates
        total_runoff = cell.runoff_surf + cell.q_base.sum(axis=1)
        total = surf + base

        data: Dict[str, np.ndarray] = {
            "DATE": dates,
            "YEAR": dates.year,
            "MONTH": dates.month,
            "DAY": dates.day,
            "TOTAL_RUNOFF": total_runoff,
            "RUNOFF_SURF": cell.runoff_surf,
            "BASEFLOW": cell.q_base.sum(axis=1),
        }
        safe_total = np.where(np.abs(total_runoff) < 1.0e-12, np.nan, total_runoff)
        for k, s in enumerate(sources):
            data[f"RUNOFF_{s}_TOTAL"] = total[:, k]
            data[f"RUNOFF_{s}_SURFACE"] = surf[:, k]
            data[f"RUNOFF_{s}_SUBSURFACE"] = base[:, k]
        for k, s in enumerate(sources):
            data[f"F_{s}_TOTAL"] = total[:, k] / safe_total
            data[f"F_{s}_SURFACE"] = surf[:, k] / safe_total
            data[f"F_{s}_SUBSURFACE"] = base[:, k] / safe_total
        for k, s in enumerate(sources):
            data[f"FW_{s}"] = fw[:, k]
        if liq_out is not None:
            for n in range(N):
                for k, s in enumerate(sources):
                    data[f"W{n + 1}_{s}_LIQ"] = liq_out[:, n, k]
                    if ice_out is not None:
                        data[f"W{n + 1}_{s}_ICE"] = ice_out[:, n, k]
        for n in range(N):
            data[f"BUDGET_RESIDUAL_LIQ{n + 1}"] = resid_liq[:, n]
            if cell.has_ice:
                data[f"BUDGET_RESIDUAL_ICE{n + 1}"] = resid_ice[:, n]
        data["CLIPPED_MASS"] = clipped
        data["UNKNOWN_INJECTED"] = injected
        data["RUNOFF_SOURCE_CLOSURE_ERROR"] = total.sum(axis=1) - total_runoff
        if liq_out is not None:
            for n in range(N):
                data[f"W{n + 1}_LIQ_CLOSURE_ERROR"] = liq_out[:, n, :].sum(axis=1) - cell.w_liq[:, n]
                if ice_out is not None:
                    data[f"W{n + 1}_ICE_CLOSURE_ERROR"] = ice_out[:, n, :].sum(axis=1) - cell.w_ice[:, n]

        out = pd.DataFrame(data)
        if cell.extra is not None:
            for c in cell.extra.columns:
                if c not in out.columns:
                    out[c] = cell.extra[c].to_numpy()
        if cfg.spinup_years > 0:
            first = pd.Timestamp(int(dates[0].year) + cfg.spinup_years, 1, 1)
            out = out.loc[out["DATE"] >= first].reset_index(drop=True)
        return out
