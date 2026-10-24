"""B2 step 2: market towns and travel time from every cropland cell and village.

Graph (--graph): `full` (PRIMARY, researcher decision 2026-10-24: all Geofabrik drive roads, 02b) or `major`
(sensitivity, major roads only, 02).
Towns (deviation from prereg, researcher decision 2026-10-24: GeoNames has no taluka HQ for Sindh):
  PRIMARY `hq_osm` = district HQ (GeoNames PPLA / PPLA2, admin1 05 = Sindh) + OSM place=town or city (Geofabrik
  2026-07-12 places layer) inside Sindh, as the taluka-level market-town approximation.
  Sensitivity: `hq` alone; `hq_pop{20,50,100}k` = HQ + GeoNames places with population >= threshold.
Market access M = min over towns of network travel time (minutes) from the cell/village TO the town (Dijkstra on
the reversed graph), snapped to the nearest graph node (snap untimed, distance stored; PB7).
Gravity G(h) = sum_j pop_j exp(-ln2 / h * t_ij) over the PRIMARY towns, pop_j = WorldPop within 5 km;
half-lives h = 30 / 60 / 120 min are PLACEHOLDERS (PB2 unsourced) and are not used in tests.
Writes data/b2_towns.parquet, data/b2_{cell,village}_market[_major].parquet.
"""
import argparse
import os

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask as rmask
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

from _common import DATA, RHAI, SWAI, PAPER_A, VID, CID, log_env, set_seeds

ap = argparse.ArgumentParser(); ap.add_argument("--graph", default="full", choices=["full", "major"])
args = ap.parse_args()
set_seeds(); log_env(__file__)
GDIR = os.environ.get("GEOFABRIK_DIR", r"C:\Users\Abbasi\Desktop\pakistan-260712-free.shp")
THRESH = (20000, 50000, 100000)
HALF_LIVES = (30, 60, 120)          # PLACEHOLDER (PB2)
pre = "b2full" if args.graph == "full" else "b2"
suf = "" if args.graph == "full" else "_major"

N = pd.read_parquet(DATA / f"{pre}_graph_nodes.parquet"); E = pd.read_parquet(DATA / f"{pre}_graph_edges.parquet")
n = len(N)
Gr = csr_matrix((E.travel_s.to_numpy(), (E.v.to_numpy(), E.u.to_numpy())), shape=(n, n))
tree = cKDTree(np.c_[N.x, N.y])
sindh = gpd.read_file(RHAI / "data/sindh_boundary.gpkg").to_crs(4326).union_all()

# ---- towns ---------------------------------------------------------------------------------------------
cols = ["geonameid", "name", "asciiname", "alt", "lat", "lon", "fclass", "fcode", "cc", "cc2", "a1", "a2", "a3",
        "a4", "pop", "elev", "dem", "tz", "mod"]
g = pd.read_csv(RHAI / "data/geonames_PK.txt", sep="\t", names=cols, dtype=str, low_memory=False,
                usecols=["geonameid", "asciiname", "lat", "lon", "fclass", "fcode", "a1", "pop"])
g = g[g.fclass == "P"].assign(pop=lambda d: d["pop"].astype(int))
g = gpd.GeoDataFrame(g, geometry=gpd.points_from_xy(g.lon.astype(float), g.lat.astype(float)), crs=4326)
g = g[g.within(sindh) & (g.a1 == "05")]          # a1 05 = Sindh (drops e.g. Hub, Balochistan)
g["hq"] = g.fcode.isin(["PPLA", "PPLA2", "PPLA3"])
g = g[g.hq | (g["pop"] >= min(THRESH))]
gn = pd.DataFrame({"town_id": "gn" + g.geonameid, "name": g.asciiname, "source": "geonames", "fcode": g.fcode,
                   "pop_gn": g["pop"], "hq": g.hq, "osm_town": False, "geometry": g.geometry})
o = gpd.read_file(os.path.join(GDIR, "gis_osm_places_free_1.shp"), engine="pyogrio")
o = o[o.fclass.isin(["town", "city"]) & o.within(sindh)]
om = pd.DataFrame({"town_id": "osm" + o.osm_id.astype(str), "name": o.name, "source": "osm", "fcode": o.fclass,
                   "pop_gn": 0, "hq": False, "osm_town": True, "geometry": o.geometry})
T = gpd.GeoDataFrame(pd.concat([gn, om], ignore_index=True), crs=4326).to_crs(32642)
with rasterio.open(SWAI / "sindh_population.tif") as r:
    wp = []
    for geom in T.buffer(5000).to_crs(r.crs):
        a, _ = rmask(r, [geom], crop=True, filled=True, nodata=0)
        a = a.astype(float); a[~np.isfinite(a) | (a < 0)] = 0
        wp.append(a.sum())
T["worldpop_5km"] = wp
td, ti = tree.query(np.c_[T.geometry.x, T.geometry.y])
T["node"], T["snap_m"] = ti, td
sets = {"hq_osm": (T.hq | T.osm_town).to_numpy(), "hq": T.hq.to_numpy()} | \
       {f"hq_pop{t // 1000}k": (T.hq | (T.pop_gn >= t)).to_numpy() for t in THRESH}
for k, m in sets.items():
    T[f"in_{k}"] = m
if args.graph == "full":
    T.drop(columns="geometry").assign(x=T.geometry.x, y=T.geometry.y).to_parquet(DATA / "b2_towns.parquet", index=False)
print({k: int(m.sum()) for k, m in sets.items()}, f"town snap max {td.max():.0f} m, >1 km: {(td > 1000).sum()}")

# ---- targets: cells and villages -------------------------------------------------------------------------
C = pd.read_parquet(DATA / "b1_cells.parquet")
V = pd.read_parquet(PAPER_A / "data/village_master.parquet")[[VID, "lon", "lat"]]
vp = gpd.GeoSeries(gpd.points_from_xy(V.lon, V.lat), crs=4326).to_crs(32642)
cs, ci = tree.query(np.c_[C.x, C.y]); vs, vi = tree.query(np.c_[vp.x, vp.y])
need = np.unique(np.r_[ci, vi]); pos = np.searchsorted(need, np.r_[ci, vi])
res = {}
for k, m in sets.items():
    d = dijkstra(Gr, directed=True, indices=np.unique(T.node.to_numpy()[m]), min_only=True)
    res[f"M_{k}"] = d[need] / 60
prim = np.flatnonzero(sets["hq_osm"])
Dm = np.empty((len(prim), len(need)), np.float32)
for r_, j in enumerate(prim):                       # per-town times for gravity and nearest-town id
    Dm[r_] = dijkstra(Gr, directed=True, indices=int(T.node.iat[j]))[need] / 60
pw = T.worldpop_5km.to_numpy()[prim]
for h in HALF_LIVES:
    res[f"G_h{h}_PLACEHOLDER"] = (pw[:, None] * np.exp(-np.log(2) / h * Dm.astype(float))).sum(0)
res["nearest_town"] = T.town_id.to_numpy()[prim][Dm.argmin(0)]
R = pd.DataFrame(res)
nC = len(C)
CM = pd.concat([C[[CID]].reset_index(drop=True), pd.DataFrame({"snap_m": cs, "node": ci}),
                R.iloc[pos[:nC]].reset_index(drop=True)], axis=1)
VM = pd.concat([V[[VID]].reset_index(drop=True), pd.DataFrame({"snap_m": vs, "node": vi}),
                R.iloc[pos[nC:]].reset_index(drop=True)], axis=1)
CM.to_parquet(DATA / f"b2_cell_market{suf}.parquet", index=False)
VM.to_parquet(DATA / f"b2_village_market{suf}.parquet", index=False)
pd.set_option("display.width", 220)
for name, df in (("cells", CM), ("villages", VM)):
    print(name, f"unreachable (inf) M_hq_osm: {int(np.isinf(df.M_hq_osm).sum())}")
    print(df[[c for c in df.columns if c.startswith(("M_", "snap"))]].replace(np.inf, np.nan)
          .describe(percentiles=[.5, .9, .99]).round(1).to_string())
