"""B1: cell-village linkage and village-level recovery table.

Cells: flood-recovery-sindh analysis_grid.gpkg (555,309 cropland cells, 300 m, EPSG:32642), centroids.
Villages: paper_a/data/village_master.parquet (5,159, keys checked in paper_a/data/inspect_report.txt).
Links: (a) every cell centroid within 2 km of a village point (SWAI 2 km buffer rule; a cell can link to
several villages); (b) nearest village with its distance (cell-level assignment; prereg uses it within 5 km).
Village outcomes (2 km links):
  n_cells: cropland cells; n_matched / n_treated: matched cells / treated matched cells
  gap_k (k = -1..8): mean over treated matched cells of [treated ndvi_anom - its pair's control ndvi_anom]
    (flood-attributable, only where treated cells exist); gap_P = mean over k 6..8; gap_PK = mean over k 6, 8
  ttr_treated: mean ttr_months over treated matched cells (base metric; censoring share reported)
  L: descriptive persistent crop loss over ALL cropland cells = mean integrated NDVI (Kharif 2025, 2026) minus
    integrated NDVI Kharif 2021 (phenology.parquet); negative = loss. Not flood-attributable (prereg §4).
Ecological caveat: village values are means over nearby cells, not village farms; small-cell villages are
flagged (PB6: >= 10 cells primary, 5 / 25 sensitivity) and excluded transparently downstream.
Writes data/b1_cells.parquet, data/b1_cell_village_2km.parquet, data/b1_cell_nearest_village.parquet,
data/b1_village_recovery.parquet, results/tables/b1_coverage_by_district.csv.
"""
import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from _common import DATA, FLOOD, PAPER_A, TABLES, VID, CID, log_env, set_seeds

set_seeds(); log_env(__file__)
R2K, THR = 2000.0, (5, 10, 25)

cells_path = DATA / "b1_cells.parquet"
if cells_path.exists():
    C = pd.read_parquet(cells_path)
else:
    g = gpd.read_file(FLOOD / "data/processed/analysis_grid.gpkg", engine="pyogrio")
    c = g.geometry.centroid
    ll = gpd.GeoSeries(c, crs=32642).to_crs(4326)
    C = pd.DataFrame({CID: g[CID].to_numpy(), "district": g["district"].to_numpy(), "x": c.x.to_numpy(),
                      "y": c.y.to_numpy(), "lon": ll.x.to_numpy(), "lat": ll.y.to_numpy()})
    C.to_parquet(cells_path, index=False)
print(f"cells {len(C):,}")

V = pd.read_parquet(PAPER_A / "data/village_master.parquet")[[VID, "district", "lon", "lat"]]
vp = gpd.GeoSeries(gpd.points_from_xy(V.lon, V.lat), crs=4326).to_crs(32642)
V["x"], V["y"] = vp.x.to_numpy(), vp.y.to_numpy()
tv = cKDTree(np.c_[V.x, V.y])
cxy = np.c_[C.x, C.y]

d, j = tv.query(cxy)
NV = pd.DataFrame({CID: C[CID], VID: V[VID].to_numpy()[j], "dist_m": d})
NV.to_parquet(DATA / "b1_cell_nearest_village.parquet", index=False)
pairs = tv.query_ball_point(cxy, R2K)
li = np.repeat(np.arange(len(C)), [len(p) for p in pairs]); lj = np.concatenate([np.asarray(p, int) for p in pairs])
L2 = pd.DataFrame({CID: C[CID].to_numpy()[li], VID: V[VID].to_numpy()[lj],
                   "dist_m": np.hypot(cxy[li, 0] - V.x.to_numpy()[lj], cxy[li, 1] - V.y.to_numpy()[lj])})
L2.to_parquet(DATA / "b1_cell_village_2km.parquet", index=False)
print(f"cells within 2 km of a village: {L2[CID].nunique():,} ({L2[CID].nunique()/len(C):.1%}); "
      f"links {len(L2):,}; nearest village <= 5 km: {(NV.dist_m <= 5000).mean():.1%}")

# ---- cell-level outcomes ------------------------------------------------------------------
P = pd.read_parquet(FLOOD / "data/processed/did_panel.parquet",
                    columns=[CID, "event_time", "ndvi_anom", "treated_t3", "matched_pair_id"])
P["event_time"] = P.event_time.astype(int)
t = P[P.treated_t3 == 1].set_index(["matched_pair_id", "event_time"])
c = P[P.treated_t3 == 0].set_index(["matched_pair_id", "event_time"])["ndvi_anom"]
t = t.assign(gap=t["ndvi_anom"] - c.reindex(t.index).to_numpy()).reset_index()
G = t.pivot(index=CID, columns="event_time", values="gap")
G.columns = [f"gap_{k}" if k >= 0 else f"gap_m{-k}" for k in G.columns]
G["gap_P"] = G[["gap_6", "gap_7", "gap_8"]].mean(axis=1); G["gap_PK"] = G[["gap_6", "gap_8"]].mean(axis=1)
M = pd.read_parquet(FLOOD / "data/processed/master_dataset.parquet",
                    columns=[CID, "ttr_months", "ttr_censored", "treated_t3", "matched_pair_id"])
M["matched"] = M.matched_pair_id.notna()
ph = pd.read_parquet(FLOOD / "data/processed/phenology.parquet", columns=[CID, "season", "integrated_ndvi"],
                     filters=[("season", "in", ["Kharif_2021", "Kharif_2025", "Kharif_2026"])])
ph = ph.pivot(index=CID, columns="season", values="integrated_ndvi")
cellL = (ph[["Kharif_2025", "Kharif_2026"]].mean(axis=1) - ph["Kharif_2021"]).rename("L")
CO = M.set_index(CID).join(G).join(cellL)
CO.to_parquet(DATA / "b1_cell_outcomes.parquet")

# ---- village aggregation ------------------------------------------------------------------
X = L2.join(CO, on=CID)
tr = X[X.matched & (X.treated_t3 == 1)]
gcols = [c for c in G.columns]
VR = pd.DataFrame({"n_cells": X.groupby(VID).size(), "n_matched": X[X.matched].groupby(VID).size(),
                   "n_treated": tr.groupby(VID).size(), "L": X.groupby(VID)["L"].mean(),
                   "ttr_treated": tr.groupby(VID)["ttr_months"].mean(),
                   "ttr_censored_share": tr.groupby(VID)["ttr_censored"].mean()}).join(tr.groupby(VID)[gcols].mean())
VR = V[[VID, "district"]].set_index(VID).join(VR)
for k in ("n_cells", "n_matched", "n_treated"):
    VR[k] = VR[k].fillna(0).astype(int)
VR.reset_index().to_parquet(DATA / "b1_village_recovery.parquet", index=False)

cov = VR.groupby("district").agg(villages=("n_cells", "size"), **{
    f"cells_ge{n}": ("n_cells", lambda s, n=n: int((s >= n).sum())) for n in THR},
    **{f"treated_ge{n}": ("n_treated", lambda s, n=n: int((s >= n).sum())) for n in THR})
cov.loc["TOTAL"] = cov.sum()
cov.to_csv(TABLES / "b1_coverage_by_district.csv")
pd.set_option("display.width", 200)
print(cov.to_string())
print(VR[["n_cells", "n_treated", "L", "gap_P", "ttr_treated"]].describe().round(3).to_string())
