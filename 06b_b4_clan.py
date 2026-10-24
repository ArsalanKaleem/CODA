"""B4 descriptive follow-up: modifier means by CATE quintile (CLAN-style, Chernozhukov et al.) on the held-out
CATEs of 06 (full_t2000), and Spearman correlation of CATE with each modifier. Descriptive: modifiers are
correlated and non-random, so this says WHERE effects are larger, not WHY.
Writes results/tables/b4_clan.csv.
"""
import numpy as np
import pandas as pd

from _common import DATA, FLOOD, TABLES, CID, log_env

log_env(__file__)
C = pd.read_parquet(DATA / "b4_cate_full_t2000.parquet")
M = pd.read_parquet(FLOOD / "data/processed/master_dataset.parquet").set_index(CID)
cm = pd.read_parquet(DATA / "b2_cell_market.parquet").set_index(CID)
cf = pd.read_parquet(DATA / "b2_cell_market_flood.parquet").set_index(CID)
fl = pd.read_parquet(DATA / "b3_cell_flood_flags.parquet").set_index(CID)
X = C.set_index(CID).join(M[["district", "treated_t3", "elevation_m", "sbi_persist_months", "sm_persist_months",
                              "flood_duration", "dist_to_road_m", "ndwi_persist_months"]])
X = X.join(cm[["M_hq_osm"]]).join(cf[["isolated_F", "dM_hq_osm"]]).join(fl[["flag_2025"]])
X["quintile"] = X.groupby("fit_half").cate.transform(lambda s: pd.qcut(s, 5, labels=False) + 1)
cols = ["cate", "M_hq_osm", "elevation_m", "sbi_persist_months", "sm_persist_months", "dist_to_road_m",
        "treated_t3", "flood_duration", "ndwi_persist_months", "isolated_F", "flag_2025"]
T = X.groupby("quintile")[cols].mean()
T.loc["spearman_with_cate"] = [X[c].astype(float).corr(X.cate, method="spearman") for c in cols]
T.to_csv(TABLES / "b4_clan.csv")
top = X[X.quintile == 1].district.value_counts(normalize=True).head(6)
pd.set_option("display.width", 220)
print(T.round(3).to_string()); print("\nmost negative CATE quintile, district shares:\n" + top.round(3).to_string())
