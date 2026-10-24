"""2025 re-flooding flag per cropland cell (prereg §1).

Raster: paper_a/data/raw/s1_flood/s1_flood_{2022,2025}_30m.tif (Paper A v1 mask; EPSG:32642, 30 m, bands
detections / flood_days / obs_days). The 300 m analysis cells are aligned to the 30 m pixel grid, so each cell is an
exact 10 x 10 pixel block. Rule (prereg): flagged if >= 50 % of observed pixels (obs_days > 0) have
detections >= 3. The same rule applied to 2022 is compared with the base treatment `treated_t3` as a check of the
aggregation (base flood_duration came from a different GEE export, so agreement need not be perfect).
Writes data/b3_cell_flood_flags.parquet and results/tables/b3_reflood_check.csv.
"""
import numpy as np
import pandas as pd
import rasterio

from _common import DATA, FLOOD, PAPER_A, TABLES, CID, log_env, set_seeds

set_seeds(); log_env(__file__)
C = pd.read_parquet(DATA / "b1_cells.parquet")
M = pd.read_parquet(FLOOD / "data/processed/master_dataset.parquet", columns=[CID, "treated_t3", "flood_duration", "matched_pair_id"])
out = C[[CID]].copy()
for y in (2022, 2025):
    with rasterio.open(PAPER_A / f"data/raw/s1_flood/s1_flood_{y}_30m.tif") as r:
        det = r.read(1); obs = r.read(3)
        left, top = r.transform.c, r.transform.f
    c0 = np.rint((C.x.to_numpy() - 150 - left) / 30).astype(int)       # cell minx -> pixel column
    r0 = np.rint((top - (C.y.to_numpy() + 150)) / 30).astype(int)      # cell maxy -> pixel row
    assert (np.abs((C.x.to_numpy() - 150 - left) / 30 - c0) < 1e-6).all(), "cell grid not aligned to raster"
    blk = np.add.outer(np.arange(10), np.zeros(10, int))
    rows = r0[:, None, None] + blk[None]; cols = c0[:, None, None] + blk.T[None]
    d, o = det[rows, cols], obs[rows, cols]
    nobs = (o > 0).sum((1, 2))
    share = np.where(nobs > 0, ((d >= 3) & (o > 0)).sum((1, 2)) / np.maximum(nobs, 1), np.nan)
    out[f"share_det3_{y}"] = share; out[f"flag_{y}"] = share >= 0.5; out[f"n_obs_px_{y}"] = nobs
X = out.merge(M, on=CID)
ct = pd.crosstab(X.treated_t3, X.flag_2022)
agree = (X.treated_t3.astype(bool) == X.flag_2022).mean()
m = X[X.matched_pair_id.notna()]
pairs_any = m.groupby("matched_pair_id").flag_2025.any()
rows = [{"check": "2022 rule vs treated_t3, agreement (all cells)", "value": agree},
        {"check": "2022 rule flagged share", "value": X.flag_2022.mean()},
        {"check": "treated_t3 share", "value": X.treated_t3.mean()},
        {"check": "2025 flagged share (all cells)", "value": X.flag_2025.mean()},
        {"check": "2025 flagged share, matched treated", "value": m[m.treated_t3 == 1].flag_2025.mean()},
        {"check": "2025 flagged share, matched controls", "value": m[m.treated_t3 == 0].flag_2025.mean()},
        {"check": "matched pairs with either member flagged 2025", "value": pairs_any.mean()},
        {"check": "cells with no observed pixels 2025", "value": (X.n_obs_px_2025 == 0).mean()}]
R = pd.DataFrame(rows); R.to_csv(TABLES / "b3_reflood_check.csv", index=False)
out.to_parquet(DATA / "b3_cell_flood_flags.parquet", index=False)
print(ct); print(R.round(4).to_string(index=False))
