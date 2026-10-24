"""B2 step 1: Sindh-wide major-road graph from the Geofabrik extract in flood-recovery-sindh (read only).

Source: flood-recovery-sindh/data/raw/osm/osm_roads.gpkg (12,531 lines, Geofabrik Pakistan; classes motorway
to tertiary and their links; fields fclass, oneway F/T/B, bridge, layer).
Build (as paper_a/19_a2_dadu_geofabrik.py): lines noded with shapely.union_all (every crossing becomes a
junction; over- and underpasses become false junctions, which can only add connectivity), RHAI class speeds
(rhai_05_route_by_district.py HWY_SPEEDS; `_link` classes take the parent class speed; maxspeed tags ignored as
in Paper A step 19), one-way from `oneway` (F forward, T reverse, B both), largest STRONGLY connected component (a weak
component left 764 Ghotki cells snapped to one-way dead ends that cannot reach any town; decision log 2026-10-24).
Edge travel time = length / speed. Off-network leg from a cell or village to its snapped node is untimed
(RHAI and Paper A convention, PB7); the snap distance is stored so models can control for it.
Writes data/b2_graph_nodes.parquet (node, x, y EPSG:32642) and data/b2_graph_edges.parquet (u, v, length_m,
travel_s, fclass, bridge).
"""
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from _common import DATA, FLOOD, log_env, set_seeds

set_seeds(); log_env(__file__)
SPEED = {"motorway": 90, "trunk": 70, "primary": 55, "secondary": 45, "tertiary": 35, "unclassified": 25,
         "residential": 20, "track": 15, "path": 8, "service": 20}

roads = gpd.read_file(FLOOD / "data/raw/osm/osm_roads.gpkg", engine="pyogrio").explode(index_parts=False).to_crs(32642)
roads["highway"] = roads["fclass"].str.replace("_link", "", regex=False)
print(f"lines {len(roads):,}, length {roads.length.sum()/1000:,.0f} km; classes not in speed table: "
      f"{sorted(set(roads.highway) - set(SPEED))}")

noded = shapely.union_all(roads.geometry.to_numpy())
segs = gpd.GeoDataFrame(geometry=list(shapely.get_parts(noded)), crs=32642)
segs = segs[segs.length > 0.5].reset_index(drop=True)
mid = gpd.GeoDataFrame(geometry=segs.geometry.interpolate(0.5, normalized=True), crs=32642)
j = gpd.sjoin_nearest(mid, roads[["highway", "fclass", "oneway", "bridge", "geometry"]], how="left", max_distance=1.0)
j = j[~j.index.duplicated(keep="first")]
segs = segs.join(j[["highway", "fclass", "oneway", "bridge"]])
print(f"noded segments {len(segs):,}; without attributes {segs.highway.isna().sum()}")
segs["speed_kph"] = segs.highway.map(SPEED).fillna(SPEED["unclassified"])

# nodes = segment end points (rounded to 1 cm)
a = np.array([g.coords[0] for g in segs.geometry]); b = np.array([g.coords[-1] for g in segs.geometry])
keys = np.round(np.r_[a, b], 2)
uniq, inv = np.unique(keys, axis=0, return_inverse=True)
u, v = inv[:len(segs)], inv[len(segs):]
L = segs.length.to_numpy(); tt = L / (segs.speed_kph.to_numpy() * 1000 / 3600)
ow = segs.oneway.astype(str).to_numpy()
fw, rv = ow != "T", ow != "F"
E = pd.concat([pd.DataFrame({"u": u[fw], "v": v[fw], "length_m": L[fw], "travel_s": tt[fw],
                             "fclass": segs.fclass.to_numpy()[fw], "bridge": segs.bridge.to_numpy()[fw]}),
               pd.DataFrame({"u": v[rv], "v": u[rv], "length_m": L[rv], "travel_s": tt[rv],
                             "fclass": segs.fclass.to_numpy()[rv], "bridge": segs.bridge.to_numpy()[rv]})],
              ignore_index=True)
E = E[E.u != E.v]
n = len(uniq)
ncomp, lab = connected_components(coo_matrix((np.ones(len(E)), (E.u, E.v)), shape=(n, n)), directed=True,
                                  connection="strong")
big = np.bincount(lab).argmax()
keep = lab == big
print(f"nodes {n:,}, strong components {ncomp:,}; largest keeps {keep.mean():.1%} of nodes")
newid = -np.ones(n, int); newid[keep] = np.arange(keep.sum())
E = E[keep[E.u] & keep[E.v]].assign(u=lambda d: newid[d.u], v=lambda d: newid[d.v])
N = pd.DataFrame({"node": np.arange(keep.sum()), "x": uniq[keep, 0], "y": uniq[keep, 1]})
N.to_parquet(DATA / "b2_graph_nodes.parquet", index=False); E.to_parquet(DATA / "b2_graph_edges.parquet", index=False)
print(f"graph: nodes {len(N):,}, directed edges {len(E):,}, km (two-way counted twice) {E.length_m.sum()/1000:,.0f}")
