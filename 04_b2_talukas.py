"""Taluka (tehsil) units for the pre-registered alternative clustering (deviation: GADM 4.1 has no taluka level).

Source: Geofabrik Pakistan 2026-07-12 gis_osm_adminareas_a_free_1.shp, fclass admin_level7 (tehsil/taluka in
Pakistan OSM; Karachi's level-7 units are towns), polygons whose representative point is inside Sindh. ODbL.
Each cell centroid and village is assigned to the polygon containing it; unassigned points (gaps/overlaps
between OSM polygons) go to the nearest polygon within 2 km, else remain missing (count reported).
The researcher chose OCHA tehsil boundaries; OSM level 7 was found locally and is used pending their
confirmation (decision log 2026-10-24); OCHA can replace it without other changes.
Writes data/b2_talukas.gpkg, data/b2_cell_taluka.parquet, data/b2_village_taluka.parquet.
"""
import os

import geopandas as gpd
import numpy as np
import pandas as pd

from _common import DATA, RHAI, PAPER_A, VID, CID, log_env, set_seeds

set_seeds(); log_env(__file__)
GDIR = os.environ.get("GEOFABRIK_DIR", r"C:\Users\Abbasi\Desktop\pakistan-260712-free.shp")
a = gpd.read_file(os.path.join(GDIR, "gis_osm_adminareas_a_free_1.shp"), engine="pyogrio")
s = gpd.read_file(RHAI / "data/sindh_boundary.gpkg").to_crs(4326).union_all()
t = a[(a.fclass == "admin_level7") & a.representative_point().within(s)].to_crs(32642).reset_index(drop=True)
t["taluka_id"] = "t" + t.osm_id.astype(str)
t[["taluka_id", "osm_id", "name", "geometry"]].to_file(DATA / "b2_talukas.gpkg", driver="GPKG")
ov = gpd.overlay(t[["taluka_id", "geometry"]], t[["taluka_id", "geometry"]].rename(columns={"taluka_id": "t2"}),
                 how="intersection", keep_geom_type=True)
ovk = ov[ov.taluka_id < ov.t2].area.sum() / 1e6
print(f"talukas {len(t)}, area {t.area.sum()/1e6:,.0f} km2, pairwise overlap {ovk:,.1f} km2")


def assign(df, key):
    p = gpd.GeoDataFrame(df[[key]], geometry=gpd.points_from_xy(df.x, df.y), crs=32642)
    j = gpd.sjoin(p, t[["taluka_id", "geometry"]], how="left", predicate="within")
    j = j[~j.index.duplicated(keep="first")]
    miss = j.taluka_id.isna()
    if miss.any():
        nn = gpd.sjoin_nearest(p[miss], t[["taluka_id", "geometry"]], how="left", max_distance=2000)
        nn = nn[~nn.index.duplicated(keep="first")]
        j.loc[miss, "taluka_id"] = nn.taluka_id
    print(f"{key}: assigned {j.taluka_id.notna().mean():.4%}, by nearest <= 2 km {int(miss.sum() - j.taluka_id.isna().sum())}, "
          f"missing {int(j.taluka_id.isna().sum())}")
    return pd.DataFrame({key: df[key].to_numpy(), "taluka_id": j.taluka_id.to_numpy()})


C = pd.read_parquet(DATA / "b1_cells.parquet")
CT = assign(C, CID); CT.to_parquet(DATA / "b2_cell_taluka.parquet", index=False)
V = pd.read_parquet(PAPER_A / "data/village_master.parquet")[[VID, "lon", "lat"]]
vp = gpd.GeoSeries(gpd.points_from_xy(V.lon, V.lat), crs=4326).to_crs(32642)
VT = assign(V.assign(x=vp.x.to_numpy(), y=vp.y.to_numpy()), VID); VT.to_parquet(DATA / "b2_village_taluka.parquet", index=False)
M = pd.read_parquet(DATA.parent.parent / "flood-recovery-sindh/data/processed/master_dataset.parquet",
                    columns=[CID, "matched_pair_id"])
mt = CT.merge(M, on=CID)
mt = mt[mt.matched_pair_id.notna()]
print(f"matched cells in {mt.taluka_id.nunique()} talukas; cells per taluka median "
      f"{mt.groupby('taluka_id').size().median():.0f}; pairs split across talukas "
      f"{(mt.groupby('matched_pair_id').taluka_id.nunique() > 1).mean():.1%}")
