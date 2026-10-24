"""B0b (POST HOC, not pre-registered; decision log 2026-10-24): base event study with district x season effects.

Motivation: base matched pairs are spatially unconstrained (median pair distance 202 km, 7.3 % same district),
so treated and control cells can face different regional shocks (canal supply, 2025 rain). Replacing season
effects by district x season effects compares treated and control cells within the same district and season.
Same panel and outcome as B0; entity effects + district-season effects (PanelOLS other_effects); SEs clustered by
district. Writes results/tables/b0b_eventstudy_ndvi_district_season_fe.csv.
"""
import numpy as np
import pandas as pd
from linearmodels.panel import PanelOLS

from _common import FLOOD, TABLES, log_env, set_seeds

set_seeds(); log_env(__file__)
P = pd.read_parquet(FLOOD / "data/processed/did_panel.parquet")
P["event_time"] = P.event_time.astype(int)
ks = sorted(k for k in P.event_time.unique() if k != -2)
nm = lambda k: f"et_m{abs(k)}" if k < 0 else f"et_p{k}"
for k in ks:
    P[nm(k)] = ((P.treated_t3 == 1) & (P.event_time == k)).astype(int)
P["cell_id_idx"] = P.cell_id.astype(str); P["time_idx"] = P.event_time
d = P.set_index(["cell_id_idx", "time_idx"])
both = d.groupby("district").treated_t3.agg(["min", "max"])
print("districts with both treated and control cells:", int(((both["min"] == 0) & (both["max"] == 1)).sum()), "of", len(both))
ds = pd.Categorical(d["district"].astype(str) + "_" + d["season"].astype(str)).codes
rows = []
res = PanelOLS(d["ndvi_anom"], d[[nm(k) for k in ks]], entity_effects=True,
               other_effects=pd.DataFrame({"ds": ds}, index=d.index)).fit(cov_type="clustered", clusters=d["district"])
base = pd.read_csv(TABLES / "b0_eventstudy_ndvi.csv").set_index("event_time")
for k in ks:
    rows.append({"event_time": k, "coef": res.params[nm(k)], "se": res.std_errors[nm(k)], "pvalue": res.pvalues[nm(k)],
                 "coef_base_season_fe": base.loc[k, "coef"]})
out = pd.DataFrame(rows); out.to_csv(TABLES / "b0b_eventstudy_ndvi_district_season_fe.csv", index=False)
print(out.round(4).to_string(index=False))
