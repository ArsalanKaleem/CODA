"""B4 (prereg §3, H7, exploratory): heterogeneity of the flood effect with CausalForestDML.

Modifiers are NOT randomly assigned and are correlated; CATE variation along a modifier is effect modification,
not the causal effect of that modifier. `sbi_persist_months` and `sm_persist_months` are measured after the
flood (post-treatment), so heterogeneity along them is descriptive only.

Outcome Y = mean ndvi_anom over k 6-8 minus mean over k -2, -1. Treatment D = treated_t3 (matched sample).
X (modifiers): Mz (log1p M_hq_osm, standardised), water S = max(R, W) (SWAI static components), health H =
RHAI B min_any (log1p), elevation_m, sbi_persist_months, sm_persist_months. W/H from the nearest village within
5 km; cells without one are dropped (count reported). W_ctl: 7 base matching covariates + district indicators.
Honest sample splitting: matched pairs split 50/50 (seed); forest fit on one half (cross-fitting cv = 5 grouped by
pair, nuisance = HistGradientBoosting), evaluated on the other; then roles swapped. Evaluation on the held-out
half: AIPW scores from 5-fold pair-grouped cross-fitted nuisances (propensity clipped to [0.01, 0.99]);
BLP: regress score on [1, CATE - mean CATE], district-clustered; "real heterogeneity" iff slope > 0 with
p < 0.05 (prereg). GATES by CATE quintile. Overlap: cross-fitted propensity share outside [0.05, 0.95].
Usage: --n N (subsample cells for timing) --trees T --stage {fit, lodo}.
Writes data/b4_*.parquet, results/tables/b4_*.csv.
"""
import argparse
import time

import numpy as np
import pandas as pd
import statsmodels.api as sm
from econml.dml import CausalForestDML
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold

from _common import DATA, FLOOD, PAPER_A, SWAI, TABLES, CID, VID, SEED, log_env, set_seeds

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=0); ap.add_argument("--trees", type=int, default=2000)
ap.add_argument("--stage", default="fit", choices=["fit", "lodo"])
args = ap.parse_args()
set_seeds(); log_env(__file__)
rng = np.random.default_rng(SEED)
MATCH = ["elevation_m", "slope_deg", "twi_approx", "dist_to_drainage_m", "cropland_frac", "preflood_ndvi_level",
         "preflood_ndvi_trend"]
XN = ["Mz", "water_S", "health_logmin", "elevation_m", "sbi_persist_months", "sm_persist_months"]
tag = f"n{args.n}_t{args.trees}" if args.n else f"full_t{args.trees}"

# ---- cell table ---------------------------------------------------------------------------------------------
P = pd.read_parquet(FLOOD / "data/processed/did_panel.parquet", columns=[CID, "event_time", "ndvi_anom"])
post = P[P.event_time.isin([6, 7, 8])].groupby(CID).ndvi_anom.mean()
pre = P[P.event_time.isin([-2, -1])].groupby(CID).ndvi_anom.mean()
M = pd.read_parquet(FLOOD / "data/processed/master_dataset.parquet").set_index(CID)
M = M[M.matched_pair_id.notna()]
T = M[["district", "treated_t3", "matched_pair_id", *MATCH, "sbi_persist_months", "sm_persist_months"]].copy()
T["Y"] = (post - pre).reindex(T.index)
cm = pd.read_parquet(DATA / "b2_cell_market.parquet").set_index(CID)
T["Mz"] = np.log1p(cm.M_hq_osm.reindex(T.index)); T["Mz"] = (T.Mz - T.Mz.mean()) / T.Mz.std()
nv = pd.read_parquet(DATA / "b1_cell_nearest_village.parquet").set_index(CID).reindex(T.index)
sw = pd.read_csv(SWAI / "data/swai_final.csv").set_index(VID)
vm = pd.read_parquet(PAPER_A / "data/village_master.parquet").set_index(VID)
ok = nv.dist_m <= 5000
T["water_S"] = np.where(ok, nv[VID].map(sw[["R", "W"]].max(axis=1)), np.nan)
T["health_logmin"] = np.where(ok, np.log1p(nv[VID].map(vm.min_any_B)), np.nan)
n0 = len(T)
T = T.dropna(subset=["Y", *XN, *MATCH])
# keep complete pairs only (honest split by pair)
T = T[T.groupby("matched_pair_id").treated_t3.transform("size") == 2]
print(f"matched cells {n0:,}; with W and H (village <= 5 km) and complete pairs: {len(T):,} "
      f"({len(T)/n0:.1%}); pairs {T.matched_pair_id.nunique():,}")
if args.n:
    keep = rng.choice(T.matched_pair_id.unique(), args.n // 2, replace=False)
    T = T[T.matched_pair_id.isin(keep)]
T = T.reset_index()
dist_dum = pd.get_dummies(T.district, drop_first=True, dtype=float)
Xm = T[XN].to_numpy(float); Wm = np.c_[T[MATCH].to_numpy(float), dist_dum.to_numpy()]
Dv = T.treated_t3.to_numpy(int); Yv = T.Y.to_numpy(float); grp = pd.factorize(T.matched_pair_id)[0]
distc = pd.factorize(T.district)[0]


def cf(n_trees):
    return CausalForestDML(model_y=HistGradientBoostingRegressor(random_state=SEED),
                           model_t=HistGradientBoostingClassifier(random_state=SEED), discrete_treatment=True,
                           cv=5, n_estimators=n_trees, honest=True, min_samples_leaf=20, random_state=SEED,
                           n_jobs=-1)


def aipw(idx):
    """Cross-fitted AIPW scores and propensities on the evaluation half (5 folds grouped by pair)."""
    Z = np.c_[Xm[idx], Wm[idx]]; y = Yv[idx]; d = Dv[idx]
    mu1, mu0, e = np.empty(len(idx)), np.empty(len(idx)), np.empty(len(idx))
    for tr, te in GroupKFold(5).split(Z, y, grp[idx]):
        my = HistGradientBoostingRegressor(random_state=SEED).fit(np.c_[Z[tr], d[tr]], y[tr])
        mu1[te] = my.predict(np.c_[Z[te], np.ones(len(te))]); mu0[te] = my.predict(np.c_[Z[te], np.zeros(len(te))])
        e[te] = HistGradientBoostingClassifier(random_state=SEED).fit(Z[tr], d[tr]).predict_proba(Z[te])[:, 1]
    ec = np.clip(e, .01, .99)
    return mu1 - mu0 + d * (y - mu1) / ec - (1 - d) * (y - mu0) / (1 - ec), e


def blp(score, cate, cl):
    Xb = sm.add_constant(cate - cate.mean())
    r = sm.OLS(score, Xb).fit(cov_type="cluster", cov_kwds={"groups": cl})
    return {"ate_const": r.params[0], "ate_se": r.bse[0], "blp_slope": r.params[1], "blp_slope_se": r.bse[1],
            "blp_slope_p_one_sided": float(r.pvalues[1] / 2 if r.params[1] > 0 else 1 - r.pvalues[1] / 2)}


if args.stage == "fit":
    pairs = np.unique(grp); half = rng.permutation(pairs)[: len(pairs) // 2]
    A = np.isin(grp, half); halves = {"A": np.flatnonzero(A), "B": np.flatnonzero(~A)}
    rows, gates, imp, cates = [], [], [], []
    for fit_h, ev_h in (("A", "B"), ("B", "A")):
        t0 = time.time()
        fi, ei = halves[fit_h], halves[ev_h]
        est = cf(args.trees).fit(Yv[fi], Dv[fi], X=Xm[fi], W=Wm[fi], groups=grp[fi])
        tfit = time.time() - t0
        cate = est.effect(Xm[ei]).ravel()
        score, e = aipw(ei)
        b = blp(score, cate, distc[ei])
        tm = (e >= .05) & (e <= .95)                     # prereg: trimmed-overlap sensitivity
        b |= {f"trim_{k}": v for k, v in blp(score[tm], cate[tm], distc[ei][tm]).items()} | {"trim_n": int(tm.sum())}
        q = pd.qcut(cate, 5, labels=False, duplicates="drop")
        for k in np.unique(q):
            m = q == k
            gates.append({"fit_half": fit_h, "quintile": int(k) + 1, "mean_cate": cate[m].mean(),
                          "gate_aipw": score[m].mean(), "gate_se": score[m].std() / np.sqrt(m.sum()), "n": int(m.sum())})
        rows.append({"tag": tag, "fit_half": fit_h, "n_fit": len(fi), "n_eval": len(ei), "fit_seconds": tfit,
                     "cate_mean": cate.mean(), "cate_sd": cate.std(), "cate_p10": np.percentile(cate, 10),
                     "cate_p90": np.percentile(cate, 90), "ate_forest": float(est.ate(Xm[ei])), **b,
                     "prop_outside_05_95": float(((e < .05) | (e > .95)).mean()),
                     "prop_outside_01_99": float(((e < .01) | (e > .99)).mean())})
        imp.append(pd.Series(est.feature_importances_, XN, name=f"importance_fit{fit_h}"))
        cates.append(pd.DataFrame({CID: T[CID].to_numpy()[ei], "cate": cate, "aipw": score, "propensity": e, "fit_half": fit_h}))
        print(rows[-1], flush=True)
    R = pd.DataFrame(rows); G = pd.DataFrame(gates); I = pd.concat(imp, axis=1)
    R.to_csv(TABLES / f"b4_blp_{tag}.csv", index=False); G.to_csv(TABLES / f"b4_gates_{tag}.csv", index=False)
    I.to_csv(TABLES / f"b4_importance_{tag}.csv")
    pd.concat(cates).to_parquet(DATA / f"b4_cate_{tag}.parquet", index=False)
    print(R.round(4).T.to_string()); print(G.round(4).to_string(index=False)); print(I.round(3).to_string())

    if not args.n:      # linear interaction benchmark on the same cells
        Z = pd.DataFrame(Xm, columns=XN); Z = (Z - Z.mean()) / Z.std()
        Xo = pd.concat([pd.Series(Dv, name="D"), Z, Z.mul(Dv, axis=0).add_prefix("D_x_"),
                        pd.DataFrame(T[MATCH].to_numpy(float), columns=MATCH), dist_dum.reset_index(drop=True)], axis=1)
        r = sm.OLS(Yv, sm.add_constant(Xo)).fit(cov_type="cluster", cov_kwds={"groups": distc})
        ols = pd.DataFrame({"coef": r.params, "se": r.bse, "p": r.pvalues}).loc[["D"] + [f"D_x_{c}" for c in XN]]
        ols.to_csv(TABLES / "b4_ols_benchmark.csv"); print(ols.round(4).to_string())
else:   # leave-one-district-out BLP: fit on 21 districts, evaluate on the held-out one
    rows = []
    for k, dname in enumerate(sorted(T.district.unique())):
        t0 = time.time()
        fi, ei = np.flatnonzero(T.district != dname), np.flatnonzero(T.district == dname)
        Wl = np.delete(Wm, np.s_[len(MATCH):], axis=1)          # district dummies dropped (held-out level unseen)
        est = cf(args.trees).fit(Yv[fi], Dv[fi], X=Xm[fi], W=Wl[fi], groups=grp[fi])
        cate = est.effect(Xm[ei]).ravel()
        if (Dv[ei] == 1).sum() < 50 or (Dv[ei] == 0).sum() < 50:
            rows.append({"district": dname, "n_eval": len(ei), "note": "too few treated or controls"}); continue
        Z = np.c_[Xm[ei], Wl[ei]]; y = Yv[ei]; d = Dv[ei]
        my = HistGradientBoostingRegressor(random_state=SEED); mu1, mu0, e = (np.empty(len(ei)) for _ in range(3))
        for tr, te in GroupKFold(5).split(Z, y, grp[ei]):
            mm = HistGradientBoostingRegressor(random_state=SEED).fit(np.c_[Z[tr], d[tr]], y[tr])
            mu1[te] = mm.predict(np.c_[Z[te], np.ones(len(te))]); mu0[te] = mm.predict(np.c_[Z[te], np.zeros(len(te))])
            e[te] = HistGradientBoostingClassifier(random_state=SEED).fit(Z[tr], d[tr]).predict_proba(Z[te])[:, 1]
        ec = np.clip(e, .01, .99); s = mu1 - mu0 + d * (y - mu1) / ec - (1 - d) * (y - mu0) / (1 - ec)
        r = sm.OLS(s, sm.add_constant(cate - cate.mean())).fit(cov_type="HC1")
        rows.append({"district": dname, "n_eval": len(ei), "blp_slope": r.params[1], "blp_se": r.bse[1],
                     "blp_p": r.pvalues[1], "cate_sd": cate.std(), "seconds": time.time() - t0})
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(TABLES / f"b4_lodo_{tag}.csv", index=False)
