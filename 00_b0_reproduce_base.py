"""B0: reproduce the base event study (flood-recovery-sindh src/did_day13.py) from the saved matched panel.

Panel: flood-recovery-sindh/data/processed/did_panel.parquet (cell x season, matched cells only).
Spec (as base): y_it = alpha_i + lambda_t + sum_{k != -2} beta_k treated_i 1[event_time = k];
PanelOLS, entity + time effects, SEs clustered by district. Outcomes NDVI (primary) and EVI anomalies.
Acceptance: coefficients within 0.001 of results/tables/day13_eventstudy_{ndvi,evi}.csv; SEs within 5 %.
The base applied pair-integrity and winsorisation before saving the panel; both are re-checked here, not
re-applied (re-applying would change the data).
Writes results/tables/b0_eventstudy_{ndvi,evi}.csv and b0_reproduction_check.csv.
"""
import numpy as np
import pandas as pd
from linearmodels.panel import PanelOLS

from _common import FLOOD, TABLES, log_env, set_seeds

set_seeds(); log_env(__file__)
P = pd.read_parquet(FLOOD / "data/processed/did_panel.parquet")
print(f"panel rows {len(P):,}, cells {P.cell_id.nunique():,}, pairs {P.matched_pair_id.nunique():,}")

# integrity: one treated + one control per pair in the reference season; 11 seasons per cell
ref = P[P.event_time == -2].groupby("matched_pair_id").treated_t3.agg(["count", "sum"])
bad = int(((ref["count"] != 2) | (ref["sum"] != 1)).sum())
spc = P.groupby("cell_id").size()
print(f"pairs failing integrity: {bad}; seasons per cell min {spc.min()} max {spc.max()}")
print("NDVI anomaly quantiles 0/0.001/0.999/1:", P.ndvi_anom.quantile([0, .001, .999, 1]).round(4).tolist())

rows = []
for y in ("ndvi_anom", "evi_anom"):
    df = P.dropna(subset=[y, "event_time"]).copy()
    df["event_time"] = df["event_time"].astype(int)
    ks = sorted(k for k in df.event_time.unique() if k != -2)
    nm = lambda k: f"et_m{abs(k)}" if k < 0 else f"et_p{k}"
    for k in ks:
        df[nm(k)] = ((df.treated_t3 == 1) & (df.event_time == k)).astype(int)
    df["cell_id_idx"] = df["cell_id"].astype(str)
    df["time_idx"] = df["event_time"]
    d = df.set_index(["cell_id_idx", "time_idx"])
    f = f"{y} ~ 1 + " + " + ".join(nm(k) for k in ks) + " + EntityEffects + TimeEffects"
    res = PanelOLS.from_formula(f, data=d).fit(cov_type="clustered", clusters=d["district"])
    ci = res.conf_int()
    out = pd.DataFrame([{"event_time": k, "coef": res.params[nm(k)], "se": res.std_errors[nm(k)],
                         "ci_low": ci.loc[nm(k), "lower"], "ci_high": ci.loc[nm(k), "upper"],
                         "pvalue": res.pvalues[nm(k)]} for k in ks]
                       + [{"event_time": -2, "coef": 0.0, "se": 0.0, "ci_low": 0.0, "ci_high": 0.0, "pvalue": np.nan}])
    out = out.sort_values("event_time").reset_index(drop=True)
    tag = y.split("_")[0]
    out.to_csv(TABLES / f"b0_eventstudy_{tag}.csv", index=False)
    base = pd.read_csv(FLOOD / f"results/tables/day13_eventstudy_{tag}.csv")
    m = out.merge(base, on="event_time", suffixes=("", "_base"))
    for _, r in m.iterrows():
        rows.append({"outcome": tag, "event_time": int(r.event_time), "coef": r.coef, "coef_base": r.coef_base,
                     "abs_diff": abs(r.coef - r.coef_base), "se": r.se, "se_base": r.se_base,
                     "se_rel_diff": abs(r.se / r.se_base - 1) if r.se_base else 0.0})
    print(f"{tag}: N = {res.nobs:,}")
C = pd.DataFrame(rows)
C["pass"] = (C.abs_diff <= 0.001) & (C.se_rel_diff <= 0.05)
C.to_csv(TABLES / "b0_reproduction_check.csv", index=False)
pd.set_option("display.width", 200)
print(C.round(5).to_string(index=False))
print("ACCEPTANCE:", "PASS" if C["pass"].all() else "FAIL")
