"""B5 (prereg §4, H3): do water, health and crop-loss deficits overlap more than spatial clustering alone produces?

Unit: villages with >= 10 cropland cells within 2 km (B1). Deficits (indicator 1):
  water  = bottom quintile of S = max(R, W) (SWAI static components; full SWAI score as a variant);
  health = top quintile of RHAI Scenario B min_any (variants: Scenario A; Paper A ensemble median); unreachable
           villages count as the worst value;
  crop   = bottom quintile of L (descriptive persistent crop change, Kharif 2025-26 vs 2021; variant: flood-attributable
           pair gap gap_P in villages with >= 10 treated cells).
Statistic T3 = villages with all three; also each pair and ">= 2". Nulls (9,999 each): (a) each layer permuted
independently within 25 km blocks; (b) Moran spectral randomisation surrogates of each continuous layer,
re-thresholded at the same quintile. H3 supported iff T3 > 95th percentile of BOTH nulls (prereg).
Rule vs multidimensional: rule = deficit on >= 2 of 3; multidimensional = Gaussian mixture on the three deficit
percentiles, k = 2..8 by BIC, bootstrap ARI stability (200); flag = component with the highest mean deficit
percentile (clarification 4). Flagged-by-one-only villages compared on the independent outcome (Paper A newly
isolated under the 2022 flood). All of this is descriptive co-location, not one deficit causing another.
Writes results/tables/b5_*.csv.
"""
import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture

from _common import DATA, PAPER_A, SWAI, TABLES, VID, SEED, log_env, set_seeds

set_seeds(); log_env(__file__)
B = 9999
rng = np.random.default_rng(SEED)
VR = pd.read_parquet(DATA / "b1_village_recovery.parquet").set_index(VID)
vm = pd.read_parquet(PAPER_A / "data/village_master.parquet").set_index(VID)
sw = pd.read_csv(SWAI / "data/swai_final.csv").set_index(VID)
a4 = pd.read_parquet(PAPER_A / "data/a4_village_distributions.parquet").query("config == 'normal'").set_index(VID)
pts = gpd.GeoSeries(gpd.points_from_xy(vm.lon, vm.lat), crs=4326, index=vm.index).to_crs(32642)
X = pd.DataFrame({"S": sw[["R", "W"]].max(axis=1), "SWAI": sw.water_access_score, "H_B": vm.min_any_B, "H_A": vm.min_any_A,
                  "H_ens": a4.p50_min, "L": VR.L, "gapP": VR.gap_P, "n_cells": VR.n_cells, "n_treated": VR.n_treated,
                  "x": pts.x, "y": pts.y, "district": vm.district})
fl = pd.concat([pd.read_parquet(p) for p in sorted((PAPER_A / "data/a2_flood").glob("*.parquet"))])
fl = fl[(fl.scenario == "B") & (fl["mask"] == "v1")]
nrm = fl[fl.event == "normal"].set_index(VID); f22 = fl[(fl.event == "2022") & (fl.D == 3) & (fl["T"] == 0.6)].set_index(VID)
X["newly_isolated_2022"] = (f22.min_any_r1000.isna() & nrm.min_any_r1000.notna() & ~f22.no_road_r1000.astype(bool)).reindex(X.index)
X = X[X.n_cells >= 10].copy()
print(f"villages with >= 10 cells: {len(X):,}; missing S {X.S.isna().sum()}, H_B {X.H_B.isna().sum()}, L {X.L.isna().sum()}; "
      f"newly isolated 2022: {X.newly_isolated_2022.mean():.3%}")


def deficits(df, w="S", h="H_B", c="L"):
    d = df[[w, h, c]].copy()
    hv = d[h].replace(np.inf, np.nan)
    d[h] = hv.fillna(hv.max() + 1)                     # unreachable = worst health access
    ok = d.notna().all(axis=1)
    d = d[ok]
    F = pd.DataFrame({"water": d[w] <= d[w].quantile(.2), "health": d[h] >= d[h].quantile(.8),
                      "crop": d[c] <= d[c].quantile(.2)})
    return F, d


def stats(F):
    a = F.to_numpy().astype(int); s = a.sum(1)
    return {"T3": int((s == 3).sum()), "ge2": int((s >= 2).sum()), "WH": int((a[:, 0] & a[:, 1]).sum()),
            "WC": int((a[:, 0] & a[:, 2]).sum()), "HC": int((a[:, 1] & a[:, 2]).sum())}


def msr_basis(xy):
    n = len(xy); k = 8
    _, nn = cKDTree(xy).query(xy, k + 1)
    Wm = csr_matrix((np.ones(n * k), (np.repeat(np.arange(n), k), nn[:, 1:].ravel())), shape=(n, n)).toarray()
    Wm = Wm / Wm.sum(1, keepdims=True); Wm = (Wm + Wm.T) / 2
    Hc = np.eye(n) - 1 / n
    val, vec = np.linalg.eigh(Hc @ Wm @ Hc)
    keep = np.abs(vec.sum(0)) < 1e-6 * np.sqrt(n)     # drop the constant direction
    return vec[:, keep]


def run(label, w="S", h="H_B", c="L", sub=None):
    df = X if sub is None else X[sub]
    F, Dv = deficits(df, w, h, c)
    obs = stats(F)
    xy = df.loc[F.index, ["x", "y"]].to_numpy()
    blk = pd.factorize(pd.Series(list(zip((xy[:, 0] // 25000).astype(int), (xy[:, 1] // 25000).astype(int)))))[0]
    groups = [np.flatnonzero(blk == g) for g in np.unique(blk) if (blk == g).sum() > 1]
    a = F.to_numpy().astype(int)
    na = np.empty(B, int)
    for b in range(B):
        cols = a.copy()
        for j in range(3):
            for gix in groups:
                cols[gix, j] = a[rng.permutation(gix), j]
        na[b] = int((cols.sum(1) == 3).sum())
    E = msr_basis(xy)
    vals = Dv.to_numpy(float); q = [0.2, 0.8, 0.2]; low = [True, False, True]
    coef = [E.T @ (vals[:, j] - vals[:, j].mean()) for j in range(3)]
    nb = np.empty(B, int)
    for b in range(B):
        ind = []
        for j in range(3):
            s = E @ (coef[j] * rng.choice([-1.0, 1.0], len(coef[j])))
            thr = np.quantile(s, q[j])
            ind.append(s <= thr if low[j] else s >= thr)
        nb[b] = int(np.all(ind, 0).sum())
    r = {"variant": label, "n_villages": len(F), **obs, "expected_if_independent": len(F) * 0.2 ** 3,
         "null_block_mean": na.mean(), "null_block_p95": np.percentile(na, 95), "p_block": (1 + (na >= obs["T3"]).sum()) / (B + 1),
         "null_msr_mean": nb.mean(), "null_msr_p95": np.percentile(nb, 95), "p_msr": (1 + (nb >= obs["T3"]).sum()) / (B + 1)}
    r["H3_supported"] = bool(obs["T3"] > r["null_block_p95"] and obs["T3"] > r["null_msr_p95"])
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    return r, F, Dv


rows = []
prim, Fp, Dp = run("primary (S, RHAI B, L)")
rows.append(prim)
for lab_, kw in (("water = full SWAI", {"w": "SWAI"}), ("health = RHAI A (open data)", {"h": "H_A"}),
                 ("health = Paper A ensemble median", {"h": "H_ens"})):
    rows.append(run(lab_, **kw)[0])
rows.append(run("crop = flood-attributable gap (villages >= 10 treated cells)", c="gapP", sub=X.n_treated >= 10)[0])
R = pd.DataFrame(rows); R.to_csv(TABLES / "b5_overlap_tests.csv", index=False)

# ---- rule vs multidimensional (primary layers) ----------------------------------------------------------
pct = pd.DataFrame({"water": 1 - Dp["S"].rank(pct=True), "health": Dp["H_B"].rank(pct=True), "crop": 1 - Dp["L"].rank(pct=True)})
Z = pct.to_numpy()
bic = {k: GaussianMixture(k, n_init=5, random_state=SEED).fit(Z).bic(Z) for k in range(2, 9)}
kb = min(bic, key=bic.get)
lab = GaussianMixture(kb, n_init=5, random_state=SEED).fit(Z).predict(Z)
aris = []
for b in range(200):
    i = rng.choice(len(Z), len(Z), replace=True)
    aris.append(adjusted_rand_score(lab, GaussianMixture(kb, n_init=2, random_state=SEED + b).fit(Z[i]).predict(Z)))
prof = pct.assign(comp=lab).groupby("comp").agg(n=("water", "size"), water=("water", "mean"), health=("health", "mean"), crop=("crop", "mean"))
prof["mean_deficit_pct"] = prof[["water", "health", "crop"]].mean(axis=1)
worst = prof.mean_deficit_pct.idxmax()
rule = Fp.sum(axis=1) >= 2; multi = pd.Series(lab == worst, index=Fp.index)
out = X.loc[Fp.index, "newly_isolated_2022"].astype(float)
cmp = [{"group": name, "n": int(m.sum()), "newly_isolated_2022_share": out[m].mean()}
       for name, m in (("both", rule & multi), ("rule only", rule & ~multi), ("multidimensional only", ~rule & multi),
                       ("neither", ~rule & ~multi))]
CMP = pd.DataFrame(cmp)
prof.to_csv(TABLES / "b5_profiles.csv"); CMP.to_csv(TABLES / "b5_rule_vs_multi.csv", index=False)
pd.DataFrame({"k": list(bic), "bic": list(bic.values())}).to_csv(TABLES / "b5_gmm_bic.csv", index=False)
pd.DataFrame({"ari": aris}).describe().to_csv(TABLES / "b5_gmm_stability.csv")
pd.set_option("display.width", 250)
print(R.round(4).to_string(index=False)); print(f"GMM k = {kb}; ARI mean {np.mean(aris):.3f} (p5 {np.percentile(aris, 5):.3f})")
print(prof.round(3).to_string()); print(CMP.round(4).to_string(index=False))
