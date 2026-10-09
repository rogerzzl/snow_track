"""Write the output of a toy lumped conceptual model (HBV-like) for two catchments.

The model is *not* part of snowtrack; it stands in for whatever hydrological
model you run.  Its output has its own column names and a soil box plus two
response boxes.  ``mapping.yaml`` tells snowtrack how to read it.

    python make_toy_model_output.py          # writes model_output/<catchment>.csv
    snowtrack validate -i model_output -c mapping.yaml
    snowtrack run      -i model_output -o tracked -c mapping.yaml
"""

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def run_toy_model(seed: int, years: int = 6, start: int = 2005) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    dates = pd.date_range(f"{start}-01-01", f"{start + years - 1}-12-31", freq="D")
    doy = dates.dayofyear.to_numpy()
    temp = -6 + 18 * np.sin(2 * np.pi * (doy - 105) / 365.25) + rng.normal(0, 2, len(dates))
    precip = rng.gamma(0.5, 5.0, len(dates)) * (0.5 + np.clip(np.sin(2 * np.pi * (doy - 100) / 365.25), 0, 1))
    snowpack = sm = 80.0; suz = 10.0; slz = 60.0
    fc, beta, lp, k0, k1, k2, perc_max, uzl = 200.0, 2.0, 0.7, 0.2, 0.08, 0.02, 2.5, 20.0
    rows = []
    for i, d in enumerate(dates):
        snow = precip[i] if temp[i] < 0.5 else 0.0
        rain = precip[i] - snow
        snowpack += snow
        melt = min(snowpack, max(temp[i] - 0.5, 0) * 3.5)
        snowpack -= melt
        inp = rain + melt
        qdirect = inp * 0.05                                  # a little infiltration-excess runoff
        recharge = (inp - qdirect) * (sm / fc) ** beta
        sm += inp - qdirect - recharge
        pet = max(0.2 + 0.18 * temp[i], 0.0)
        aet = min(pet * min(sm / (lp * fc), 1.0), sm)
        sm -= aet
        perc = min(perc_max, suz + recharge)
        suz += recharge - perc
        q0 = k0 * max(suz - uzl, 0.0); q1 = k1 * suz
        suz -= q0 + q1
        slz += perc
        q2 = k2 * slz; slz -= q2
        rows.append(dict(date=d, rain=rain, snowmelt=melt, qdirect=qdirect, sm=sm, suz=suz, slz=slz,
                         q_upper=q0 + q1, q_lower=q2, perc=recharge, aet=aet, swe=snowpack))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    out = HERE / "model_output"
    out.mkdir(exist_ok=True)
    for name, seed in (("catchment_A", 1), ("catchment_B", 2)):
        run_toy_model(seed).to_csv(out / f"{name}.csv", index=False, float_format="%.6f")
        print("wrote", out / f"{name}.csv")
