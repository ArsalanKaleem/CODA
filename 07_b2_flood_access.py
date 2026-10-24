"""B2 secondary (prereg §1): flood-conditioned market access M_F, 2022 event, Paper A base closure rule.

Exposure: as paper_a/05_a1_segment_exposure.py, points every 15 m along each road piece (endpoints included), band
flood_days of paper_a/data/raw/s1_flood/s1_flood_2022_30m.tif (v1 mask); points with obs_days == 0 are masked
(permanent water) and never count as inundated.
Rule unit: Paper A closes OSMnx junction-to-junction edges. The full graph (02b) has vertex-pair pieces, so
exposure is aggregated over junction-to-junction CHAINS (pieces joined through nodes of degree 2 whose two pieces
share road class and bridge flag). A chain is closed if its unmasked share with flood_days >= 3 is >= 0.6
(D3, T0.6); chains with no unmasked point stay open (as Paper A, nan -> 0); bridge chains are never closed.
M_F: min travel time to the primary town set (hq_osm) on the graph without closed pieces, from the snapped node;
if that node cannot reach a town, the best node within 1 km is used (Paper A decision 6); none -> isolated (inf).
This is an exposure scenario (radar water on the road), not observed closures.
Writes data/b2_flood_closed_segs.parquet, data/b2_{cell,village}_market_flood.parquet,
results/tables/b2_flood_access_summary.csv.
"""
import numpy as np
import pandas as pd
import rasterio
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

from _common import DATA, PAPER_A, TABLES, CID, VID, log_env, set_seeds

set_seeds(); log_env(__file__)
SP, DAYS, THR, R_ISO = 15.0, 3, 0.6, 1000.0
N = pd.read_parquet(DATA / "b2full_graph_nodes.parquet"); E = pd.read_parquet(DATA / "b2full_graph_edges.parquet")
nx_, ny_ = N.x.to_numpy(), N.y.to_numpy(); n = len(N)
S = E.drop_duplicates("seg")[["seg", "u", "v", "length_m", "fclass", "bridge"]].reset_index(drop=True)
print(f"undirected pieces {len(S):,}")

# ---- chains -----------------------------------------------------------------------------------------------
inc = pd.DataFrame({"node": np.r_[S.u, S.v], "piece": np.r_[S.index, S.index]})
deg = inc.groupby("node").size()
key = S.fclass.astype(str) + "|" + S.bridge.astype(str)
inc["key"] = key.to_numpy()[inc.piece]
nk = inc.groupby("node").key.nunique()
through = deg.index[(deg == 2) & (nk.reindex(deg.index) == 1)]
t = inc[inc.node.isin(through)].sort_values(["node", "piece"])
pa = t.piece.to_numpy()[0::2]; pb = t.piece.to_numpy()[1::2]
_, chain = connected_components(coo_matrix((np.ones(len(pa)), (pa, pb)), shape=(len(S), len(S))), directed=False)
S["chain"] = chain
print(f"chains {S.chain.nunique():,}; median pieces per chain {S.groupby('chain').size().median():.0f}")

# ---- exposure ---------------------------------------------------------------------------------------------
with rasterio.open(PAPER_A / "data/raw/s1_flood/s1_flood_2022_30m.tif") as r:
    fd = r.read(2); obs = r.read(3); tr = r.transform
x0, y0, x1, y1 = nx_[S.u], ny_[S.u], nx_[S.v], ny_[S.v]
L = np.hypot(x1 - x0, y1 - y0)
k = np.maximum(2, np.ceil(L / SP).astype(int) + 1)
pid = np.repeat(np.arange(len(S)), k)
f = (np.arange(k.sum()) - np.repeat(np.cumsum(k) - k, k)) / np.repeat(k - 1, k)
px = x0[pid] + (x1 - x0)[pid] * f; py = y0[pid] + (y1 - y0)[pid] * f
col = np.floor((px - tr.c) / tr.a).astype(int); row = np.floor((py - tr.f) / tr.e).astype(int)
ok = (row >= 0) & (row < fd.shape[0]) & (col >= 0) & (col < fd.shape[1])
v = np.zeros(len(px), np.int16); m = np.ones(len(px), bool)
v[ok] = fd[row[ok], col[ok]]; m[ok] = obs[row[ok], col[ok]] == 0
del fd, obs
ch = S.chain.to_numpy()[pid]
um = ~m
num = np.bincount(ch[um], weights=(v[um] >= DAYS).astype(float), minlength=chain.max() + 1)
den = np.bincount(ch[um], minlength=chain.max() + 1)
frac = np.where(den > 0, num / np.maximum(den, 1), 0.0)
S["chain_frac_fd_ge3"] = frac[S.chain]
S["closed"] = (S.chain_frac_fd_ge3 >= THR) & ~S.bridge.astype(bool)
S[["seg", "chain", "chain_frac_fd_ge3", "closed"]].to_parquet(DATA / "b2_flood_closed_segs.parquet", index=False)
print(f"sample points {len(px):,}; closed pieces {S.closed.mean():.2%}, closed km {S.length_m[S.closed].sum()/1000:,.0f} "
      f"of {S.length_m.sum()/1000:,.0f}")

# ---- routing ----------------------------------------------------------------------------------------------
closed_seg = set(S.seg[S.closed])
Eo = E[~E.seg.isin(closed_seg)]
Gr = csr_matrix((Eo.travel_s.to_numpy(), (Eo.v.to_numpy(), Eo.u.to_numpy())), shape=(n, n))
T = pd.read_parquet(DATA / "b2_towns.parquet")
src = np.unique(T.node[T.in_hq_osm].to_numpy())
tF = dijkstra(Gr, directed=True, indices=src, min_only=True) / 60
tree = cKDTree(np.c_[nx_, ny_])
tree_r = cKDTree(np.c_[nx_[np.isfinite(tF)], ny_[np.isfinite(tF)]]); reach_t = tF[np.isfinite(tF)]


def flood_access(base, xy):
    mf = tF[base.node.to_numpy()]
    bad = ~np.isfinite(mf)
    if bad.any():                                      # Paper A decision 6: best reachable node within 1 km
        lists = tree_r.query_ball_point(xy[bad], R_ISO)
        mf[bad] = [reach_t[l].min() if len(l) else np.inf for l in lists]
    out = base[[base.columns[0], "node", "snap_m", "M_hq_osm"]].copy()
    out["M_F_hq_osm"] = mf; out["isolated_F"] = ~np.isfinite(mf)
    out["dM_hq_osm"] = np.where(np.isfinite(mf), mf - out.M_hq_osm, np.nan)
    return out


import geopandas as gpd
C = pd.read_parquet(DATA / "b1_cells.parquet")
CM = pd.read_parquet(DATA / "b2_cell_market.parquet").set_index(CID).loc[C[CID]].reset_index()
CF = flood_access(CM, np.c_[C.x, C.y]); CF.to_parquet(DATA / "b2_cell_market_flood.parquet", index=False)
V = pd.read_parquet(PAPER_A / "data/village_master.parquet")[[VID, "lon", "lat"]]
vp = gpd.GeoSeries(gpd.points_from_xy(V.lon, V.lat), crs=4326).to_crs(32642)
VM = pd.read_parquet(DATA / "b2_village_market.parquet").set_index(VID).loc[V[VID]].reset_index()
VF = flood_access(VM, np.c_[vp.x, vp.y]); VF.to_parquet(DATA / "b2_village_market_flood.parquet", index=False)

M = pd.read_parquet(DATA.parent.parent / "flood-recovery-sindh/data/processed/master_dataset.parquet",
                    columns=[CID, "treated_t3", "matched_pair_id"])
X = CF.merge(M, on=CID)
rows = []
for name, d in (("all cells", X), ("matched treated", X[X.matched_pair_id.notna() & (X.treated_t3 == 1)]),
                ("matched controls", X[X.matched_pair_id.notna() & (X.treated_t3 == 0)]), ("villages", VF)):
    rows.append({"group": name, "n": len(d), "isolated_share": d.isolated_F.mean(),
                 "dM_median": np.nanmedian(d.dM_hq_osm), "dM_p90": np.nanpercentile(d.dM_hq_osm, 90),
                 "dM_gt10min_share": np.nanmean(d.dM_hq_osm > 10), "M_median": d.M_hq_osm.median(),
                 "M_F_median_reachable": np.nanmedian(np.where(d.isolated_F, np.nan, d.M_F_hq_osm))})
R = pd.DataFrame(rows); R.to_csv(TABLES / "b2_flood_access_summary.csv", index=False)
print(R.round(3).to_string(index=False))
