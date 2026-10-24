"""H8 (pre-registration amendment 1, frozen 2026-10-24 before estimation): flood-cut market access and the
persistent gap.

Effect modification only: cells cut off from every market town in the 2022 flood are also more flooded and
lower-lying. The model controls for flood duration (dose) x event time, normal market access and road distance,
but other confounding remains; no causal reading of h_k.

Model (matched panel, cell + season FE; district x season FE as a variant):
  y = a_i + l_t + sum_k [b_k D + h_k D Iso + e_k Iso + f_k D FDz + g_k D Mz + m_k Mz + r_k D Rz + s_k Rz] 1[k]
(FDz enters only through D x 1[k] x FDz: the 1[k] x FDz main effect is identified only by controls with partial
flood counts, and is included as well.)
Primary: h_P = mean(h_6, h_7, h_8); supported iff h_P < 0, WCB p < 0.05 and h_-1 p >= 0.10.
Writes results/tables/b8_flood_cut_access.csv.
"""
import numpy as np
import pandas as pd

from _common import DATA, FLOOD, TABLES, CID, log_env, set_seeds
from _es import EventStudy

set_seeds(); log_env(__file__)
P = pd.read_parquet(FLOOD / "data/processed/did_panel.parquet"); P["event_time"] = P.event_time.astype(int)
W = P.pivot(index=CID, columns="event_time", values="ndvi_anom"); cells = W.index.to_numpy(); ets = list(W.columns)
first = P.drop_duplicates(CID).set_index(CID).loc[cells]
D = first.treated_t3.to_numpy(float); dist = first.district.astype(str).to_numpy(); pair = first.matched_pair_id.to_numpy()
M = pd.read_parquet(FLOOD / "data/processed/master_dataset.parquet", columns=[CID, "flood_duration", "dist_to_road_m"]).set_index(CID).loc[cells]
cm = pd.read_parquet(DATA / "b2_cell_market.parquet").set_index(CID).loc[cells]
cf = pd.read_parquet(DATA / "b2_cell_market_flood.parquet").set_index(CID).loc[cells]
fl = pd.read_parquet(DATA / "b3_cell_flood_flags.parquet").set_index(CID).loc[cells]
z = lambda x: (x - x.mean()) / x.std()
Iso = cf.isolated_F.to_numpy(float)
FDz = z(M.flood_duration.to_numpy(float)); Mz = z(np.log1p(cm.M_hq_osm.to_numpy())); Rz = z(np.log1p(M.dist_to_road_m.to_numpy()))
dM = cf.dM_hq_osm.to_numpy(float); dMz = np.where(np.isfinite(dM), z(np.log1p(np.nan_to_num(dM, nan=0).clip(0))), 0.0)
print(f"Iso_F: treated {Iso[D == 1].mean():.3f}, controls {Iso[D == 0].mean():.3f}; "
      f"corr(Iso, flood_duration | treated) = {np.corrcoef(Iso[D == 1], M.flood_duration.to_numpy()[D == 1])[0, 1]:.3f}")
pr = pd.Series(fl.flag_2025.to_numpy(), index=pair).groupby(level=0).any()
noreflood = ~pd.Series(pair).map(pr).to_numpy()
KS = [k for k in ets if k != -2]
oneK = {k: (np.array(ets) == k).astype(float) for k in KS}


def extras(idx, mods, controls=True):
    ex = {}
    for name, v in mods.items():
        for k in KS:
            ex[f"h{name}_{k}"] = np.outer(D[idx] * v[idx], oneK[k]); ex[f"e{name}_{k}"] = np.outer(v[idx], oneK[k])
    if controls:
        for name, v in (("FD", FDz), ("M", Mz), ("R", Rz)):
            for k in KS:
                ex[f"c{name}_{k}"] = np.outer(D[idx] * v[idx], oneK[k]); ex[f"m{name}_{k}"] = np.outer(v[idx], oneK[k])
    return ex


rows = []


def run(label, mods, idx=None, fe="twfe", controls=True, boot=True):
    idx = np.ones(len(cells), bool) if idx is None else idx
    es = EventStudy(W.to_numpy()[idx], D[idx], ets, -2, extra=extras(idx, mods, controls),
                    fe_group=dist[idx] if fe == "dsfe" else None)
    name = list(mods)[0]
    for lab, ks in (("h_P (mean k6-8)", [6, 7, 8]), ("h_pre (k-1)", [-1]), ("h_0", [0]), ("h_3", [3])):
        R = np.zeros(len(es.names))
        for k in ks:
            R[es.names.index(f"h{name}_{k}")] = 1 / len(ks)
        t = es.test(R, dist[idx]); w = es.wcb(R, dist[idx]) if boot else {}
        rows.append({"spec": label, "fe": fe, "n_cells": int(idx.sum()), "combo": lab, "estimate": t["estimate"], "se": t["se"],
                     "p_district": t["p"], **{k: v for k, v in w.items() if k in ("p_wcb", "ci_wcb_lo", "ci_wcb_hi")}})
    R = np.zeros(len(es.names))
    for k in (6, 7, 8):
        R[es.names.index(f"b_{k}")] = 1 / 3
    t = es.test(R, dist[idx])
    rows.append({"spec": label, "fe": fe, "n_cells": int(idx.sum()), "combo": "b_P (gap, Iso=0, mods at mean)", "estimate": t["estimate"],
                 "se": t["se"], "p_district": t["p"]})
    print(pd.DataFrame(rows[-5:]).round(4).to_string(index=False), flush=True)


run("primary: Iso_F", {"Iso": Iso})
run("primary: Iso_F", {"Iso": Iso}, fe="dsfe")
run("Iso_F, excl. 2025 re-flooded pairs", {"Iso": Iso}, idx=noreflood)
run("Iso_F without FD/M/R controls", {"Iso": Iso}, controls=False)
run("secondary: dMz (added minutes)", {"dM": dMz}, idx=(Iso == 0))
T = pd.DataFrame(rows); T.to_csv(TABLES / "b8_flood_cut_access.csv", index=False)
