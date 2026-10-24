"""B2 step 1b (researcher decision 2026-10-24, "hybrid graph"): Sindh-wide ALL-drive-road graph.

The major-roads graph (02) overstated long trips by up to ~30 % vs OSRM where minor roads give shortcuts.
Source: Geofabrik Pakistan extract 2026-07-12 (path in env GEOFABRIK_DIR, default the researcher's Desktop copy;
read only), gis_osm_roads_free_1.shp, drive classes as paper_a/19_a2_dadu_geofabrik.py (OSMnx "drive"
approximation: track and paths excluded, service kept), clipped to Sindh + 10 km.
Build: one node per distinct line vertex (rounded to 1 cm); OSM ways that share a node share that vertex, so
junctions are real and over/underpasses are NOT joined (unlike shapely noding in 02). One edge per consecutive
vertex pair, RHAI class speeds (`_link` = parent class), one-way F/T/B, largest strongly connected component.
NOTE: the 2026 network postdates the flood (roads built 2022-2026 included); stated as a limitation.
Edges carry `bridge` (Geofabrik T/F) and `seg` (undirected vertex-pair id) for 07_b2_flood_access.py.
Writes data/b2full_graph_nodes.parquet, data/b2full_graph_edges.parquet.
"""
import os

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from _common import DATA, RHAI, log_env, set_seeds

set_seeds(); log_env(__file__)
GDIR = os.environ.get("GEOFABRIK_DIR", r"C:\Users\Abbasi\Desktop\pakistan-260712-free.shp")
DRIVE = {"motorway", "motorway_link", "trunk", "trunk_link", "primary", "primary_link", "secondary",
         "secondary_link", "tertiary", "tertiary_link", "unclassified", "residential", "living_street",
         "service", "road"}
SPEED = {"motorway": 90, "trunk": 70, "primary": 55, "secondary": 45, "tertiary": 35, "unclassified": 25,
         "residential": 20, "track": 15, "path": 8, "service": 20}

sindh = gpd.read_file(RHAI / "data/sindh_boundary.gpkg").to_crs(32642).buffer(10000).union_all()
sindh4326 = gpd.GeoSeries([sindh], crs=32642).to_crs(4326)
roads = gpd.read_file(os.path.join(GDIR, "gis_osm_roads_free_1.shp"), bbox=tuple(sindh4326.total_bounds),
                      engine="pyogrio", columns=["osm_id", "fclass", "oneway", "bridge"])
roads = roads[roads.fclass.isin(DRIVE)].to_crs(32642)
roads = roads[roads.intersects(sindh)].explode(index_parts=False).reset_index(drop=True)
roads["highway"] = (roads.fclass.str.replace("_link", "", regex=False)
                    .replace({"living_street": "residential", "road": "unclassified"}))
print(f"drive lines {len(roads):,}, {roads.length.sum()/1000:,.0f} km; by class:\n"
      f"{(roads.groupby('highway').apply(lambda d: d.length.sum()/1000)).round(0).to_string()}")

coords, offs = shapely.get_coordinates(roads.geometry.to_numpy(), return_index=True)
keys = np.round(coords, 2)
uniq, inv = np.unique(keys, axis=0, return_inverse=True)
same = offs[1:] == offs[:-1]                      # consecutive vertices of the same line
a, b, line = inv[:-1][same], inv[1:][same], offs[:-1][same]
L = np.hypot(*(coords[1:][same] - coords[:-1][same]).T)
spd = roads.highway.map(SPEED).to_numpy(float)[line]
tt = L / (spd * 1000 / 3600)
ow = roads.oneway.astype(str).to_numpy()[line]
cls = roads.fclass.to_numpy()[line]
brg = (roads.bridge.astype(str) == "T").to_numpy()[line]
seg = np.arange(len(a))                            # undirected segment id (for flood exposure, 07)
fw, rv = ow != "T", ow != "F"
E = pd.DataFrame({"u": np.r_[a[fw], b[rv]], "v": np.r_[b[fw], a[rv]], "length_m": np.r_[L[fw], L[rv]],
                  "travel_s": np.r_[tt[fw], tt[rv]], "fclass": np.r_[cls[fw], cls[rv]],
                  "bridge": np.r_[brg[fw], brg[rv]], "seg": np.r_[seg[fw], seg[rv]]})
E = E[(E.u != E.v) & (E.length_m > 0)]
n = len(uniq)
ncomp, lab = connected_components(coo_matrix((np.ones(len(E)), (E.u, E.v)), shape=(n, n)), directed=True,
                                  connection="strong")
keep = lab == np.bincount(lab).argmax()
print(f"vertices {n:,}, strong components {ncomp:,}; largest keeps {keep.mean():.1%}")
newid = -np.ones(n, int); newid[keep] = np.arange(keep.sum())
E = E[keep[E.u] & keep[E.v]].assign(u=lambda d: newid[d.u], v=lambda d: newid[d.v])
N = pd.DataFrame({"node": np.arange(keep.sum()), "x": uniq[keep, 0], "y": uniq[keep, 1]})
N.to_parquet(DATA / "b2full_graph_nodes.parquet", index=False); E.to_parquet(DATA / "b2full_graph_edges.parquet", index=False)
print(f"graph: nodes {len(N):,}, directed edges {len(E):,}")
