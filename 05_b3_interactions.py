"""B3 (prereg §2, H6): flood x event time x market access on the base matched panel.

Access is NOT randomly assigned: g_k measures how the flood gap varies with market access (effect modification),
not the effect of access. No causal reading of g_k.

Model (cell i, event time k, reference k = -2):
  y_it = a_i + l_t + sum_k b_k D 1[k] + sum_k g_k D 1[k] Mz + sum_k d_k 1[k] Mz + sum_k r_k D 1[k] Rz + sum_k s_k 1[k] Rz
  Mz = standardised log1p(M) over matched cells, M = M_hq_osm (full graph, deviations D-B1, D-B3);
  Rz = standardised log1p(dist_to_road_m) (base covariate).
Estimation: the panel is balanced (11 seasons per cell), so the fixed effects are removed exactly by two-way
demeaning (cell and season); post hoc variant (D-B4): cell and district x season effects (cells are nested in
districts: z - mean_i - mean_dt + mean_d). OLS on demeaned data.
Inference: CR1 by district (G = 22), taluka (OSM level 7, D-B2) and matched pair. PRIMARY p-value: wild cluster
restricted bootstrap by district (Webb 6-point weights, 9,999 draws, seed 20221001), for g_P = mean(g_6, g_7, g_8)
and the pre-trend g_-1, computed exactly from cluster-level sufficient statistics.
Acceptance check: with interactions omitted, coefficients equal B0 (b0_eventstudy_ndvi.csv) within 1e-6.
Specs: primary; no road control; Kharif-only g_PK; Rabi-only g_7; terciles (b_k by M tercile of treated cells);
2025 re-flood exclusion (pairs with either member flagged); M sensitivity sets (hq_pop50k, hq, major graph);
EVI; negative control (preflood_ndvi_trend, cross-section); all also under district x season effects where
applicable. Writes results/tables/b3_*.csv.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp

from _common import DATA, FLOOD, TABLES, CID, SEED, log_env, set_seeds

set_seeds(); log_env(__file__)
B = 9999
KS = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]
rng = np.random.default_rng(SEED)
WEBB = np.array([-np.sqrt(1.5), -1, -np.sqrt(.5), np.sqrt(.5), 1, np.sqrt(1.5)])

# ---- data (cells x 11 seasons, sorted) ----------------------------------------------------------------
P = pd.read_parquet(FLOOD / "data/processed/did_panel.parquet")
P["event_time"] = P.event_time.astype(int)
P = P.sort_values([CID, "event_time"]).reset_index(drop=True)
cells = P[CID].to_numpy()[::11]
assert (P.groupby(CID).size() == 11).all() and (P[CID].to_numpy().reshape(-1, 11) == cells[:, None]).all()
ET = P.event_time.to_numpy().reshape(-1, 11)[0]
assert (P.event_time.to_numpy().reshape(-1, 11) == ET).all()
Y = {"ndvi": P.ndvi_anom.to_numpy().reshape(-1, 11), "evi": P.evi_anom.to_numpy().reshape(-1, 11)}
cs = P.groupby(CID, sort=True).first()
M0 = pd.read_parquet(FLOOD / "data/processed/master_dataset.parquet",
                     columns=[CID, "dist_to_road_m", "preflood_ndvi_trend", "district"]).set_index(CID).loc[cells]
cm = pd.read_parquet(DATA / "b2_cell_market.parquet").set_index(CID).loc[cells]
cmj = pd.read_parquet(DATA / "b2_cell_market_major.parquet").set_index(CID).loc[cells]
tal = pd.read_parquet(DATA / "b2_cell_taluka.parquet").set_index(CID).loc[cells, "taluka_id"]
fl = pd.read_parquet(DATA / "b3_cell_flood_flags.parquet").set_index(CID).loc[cells]
D = cs.loc[cells, "treated_t3"].to_numpy(float)
pair = cs.loc[cells, "matched_pair_id"].to_numpy()
dist = cs.loc[cells, "district"].astype(str).to_numpy()
z = lambda x: (np.log1p(x) - np.log1p(x).mean()) / np.log1p(x).std()
ACC = {"M_hq_osm": z(cm.M_hq_osm.to_numpy()), "M_hq_pop50k": z(cm.M_hq_pop50k.to_numpy()),
       "M_hq": z(cm.M_hq.to_numpy()), "M_hq_osm_major": z(cmj.M_hq_osm.to_numpy())}
Rz = z(M0.dist_to_road_m.to_numpy())
pair_reflood = pd.Series(fl.flag_2025.to_numpy(), index=pair).groupby(level=0).any()
keep_noreflood = ~pd.Series(pair).map(pair_reflood).to_numpy()
print(f"cells {len(cells):,}; talukas missing {int(pd.isna(tal).sum())}; pairs re-flooded 2025 {pair_reflood.mean():.1%}")
print("Mz (hq_osm) by treatment: treated mean %.3f, control mean %.3f" % (ACC["M_hq_osm"][D == 1].mean(), ACC["M_hq_osm"][D == 0].mean()))


def demean(Z, grp, mode):
    """Z: (n, 11). mode 'twfe' (cell + season) or 'dsfe' (cell + district x season)."""
    Z = Z - Z.mean(1, keepdims=True)
    if mode == "twfe":
        return Z - Z.mean(0, keepdims=True)
    out = np.empty_like(Z)
    for g in np.unique(grp):
        m = grp == g
        out[m] = Z[m] - Z[m].mean(0, keepdims=True)
    return out


def design(mz, mode, idx, terciles=None, road=True):
    cols, names = [], []
    T1 = (ET[None, :] == np.array(KS)[:, None])                      # (10, 11)
    if terciles is None:
        for k, t in zip(KS, T1):
            cols.append(D[idx, None] * t); names.append(f"b_{k}")
        for k, t in zip(KS, T1):
            cols.append((D[idx] * mz[idx])[:, None] * t); names.append(f"g_{k}")
        for k, t in zip(KS, T1):
            cols.append(mz[idx, None] * t); names.append(f"d_{k}")
    else:
        for q in (1, 2, 3):
            tq = (terciles[idx] == q).astype(float)
            for k, t in zip(KS, T1):
                cols.append((D[idx] * tq)[:, None] * t); names.append(f"b{q}_{k}")
                if q > 1:
                    cols.append(tq[:, None] * t); names.append(f"q{q}_{k}")
    if road:
        for k, t in zip(KS, T1):
            cols.append((D[idx] * Rz[idx])[:, None] * t); names.append(f"r_{k}")
            cols.append(Rz[idx, None] * t); names.append(f"s_{k}")
    grp = dist[idx]
    X = np.stack([demean(np.broadcast_to(c, (idx.sum(), 11)).astype(float), grp, mode).ravel() for c in cols], 1)
    return X, names


def cluster_ids(name, idx):
    v = {"district": dist, "taluka": tal.to_numpy().astype(str), "pair": pair.astype(str)}[name][idx]
    return pd.factorize(np.repeat(v, 11))[0]


def fit(y2, X, names, idx, mode, combos, boot=True):
    y = demean(y2[idx], dist[idx], mode).ravel()
    XtX = X.T @ X; XtXi = np.linalg.inv(XtX)
    beta = XtXi @ (X.T @ y); u = y - X @ beta
    N, K = X.shape
    res = {"beta": pd.Series(beta, names)}
    se = {}
    for cl in ("district", "taluka", "pair"):
        g = cluster_ids(cl, idx); G = g.max() + 1
        S = sp.csr_matrix((np.ones(N), (g, np.arange(N))), shape=(G, N)) @ (X * u[:, None])
        V = XtXi @ (S.T @ S) @ XtXi * (G / (G - 1)) * ((N - 1) / (N - K))
        se[cl] = (V, G)
    out = []
    for cname, w in combos.items():
        R = np.array([w.get(n, 0.0) for n in names]); est = R @ beta
        row = {"combo": cname, "estimate": est}
        for cl, (V, G) in se.items():
            s = np.sqrt(R @ V @ R); row[f"se_{cl}"] = s
            from scipy import stats
            row[f"p_{cl}"] = 2 * stats.t.sf(abs(est / s), G - 1)
        if boot:
            row["p_wcb_district"] = wcb(y, X, XtXi, beta, R, idx)
        out.append(row)
    return res, pd.DataFrame(out)


def wcb(y, X, XtXi, beta, R, idx):
    """Wild cluster restricted bootstrap (district), Webb weights; t-stat with CR1 SE."""
    g = cluster_ids("district", idx); G = g.max() + 1; N, K = X.shape
    a = XtXi @ R; br = beta - a * (R @ beta) / (R @ a)                # restricted estimate (R b = 0)
    ur = y - X @ br
    Gm = sp.csr_matrix((np.ones(N), (g, np.arange(N))), shape=(G, N))
    Su = Gm @ (X * ur[:, None])                                          # (G, K) X_g' u_r,g
    XgXg = [X[g == j].T @ X[g == j] for j in range(G)]
    c = (G / (G - 1)) * ((N - 1) / (N - K))
    W = WEBB[rng.integers(0, 6, (B, G))]
    # beta* - br = XtXi sum_g w_g X_g'u_g ; score_g* = w_g X_g'u_g - X_g'X_g (beta* - br)
    db = (W @ Su) @ XtXi.T                                               # (B, K)
    num = db @ R                                                         # R beta*  (R br = 0)
    Ra = XtXi @ R
    vstar = np.zeros(B)
    for j in range(G):
        sj = W[:, [j]] * Su[j][None, :] - db @ XgXg[j]                   # (B, K)
        vstar += (sj @ Ra) ** 2
    tstar = num / np.sqrt(c * vstar)
    # observed t with unrestricted CR1
    u = y - X @ beta
    S = Gm @ (X * u[:, None]); t0 = (R @ beta) / np.sqrt(c * ((S @ Ra) ** 2).sum())
    return float((np.abs(tstar) >= abs(t0)).mean())


def combos_for(names, prefix="g"):
    c = {f"{prefix}_P (mean k6-8)": {f"{prefix}_{k}": 1 / 3 for k in (6, 7, 8)},
         f"{prefix}_PK (Kharif k6,8)": {f"{prefix}_{k}": 1 / 2 for k in (6, 8)},
         f"{prefix}_R (Rabi k7)": {f"{prefix}_7": 1.0},
         f"{prefix}_pre (k-1)": {f"{prefix}_-1": 1.0},
         f"{prefix}_0 (shock)": {f"{prefix}_0": 1.0}}
    return c


ALL = np.ones(len(cells), bool)
summary, coefs = [], []


def run(label, y="ndvi", acc="M_hq_osm", mode="twfe", idx=ALL, road=True, boot=True):
    X, names = design(ACC[acc], mode, idx, road=road)
    res, tab = fit(Y[y], X, names, idx, mode, combos_for(names) | {"b_P (avg gap, Mz=0)": {f"b_{k}": 1 / 3 for k in (6, 7, 8)}}, boot)
    tab.insert(0, "spec", label); tab.insert(1, "fe", mode); tab.insert(2, "n_cells", int(idx.sum()))
    summary.append(tab)
    coefs.append(res["beta"].rename(f"{label}|{mode}"))
    print(f"\n== {label} [{mode}]"); print(tab.drop(columns=["spec"]).round(5).to_string(index=False), flush=True)


# ---- acceptance: interactions omitted == B0 --------------------------------------------------------------
T1 = (ET[None, :] == np.array(KS)[:, None])
X0 = np.stack([demean(np.broadcast_to(D[:, None] * t, (len(D), 11)).astype(float), dist, "twfe").ravel() for t in T1], 1)
y0 = demean(Y["ndvi"], dist, "twfe").ravel()
b0 = np.linalg.solve(X0.T @ X0, X0.T @ y0)
base = pd.read_csv(TABLES / "b0_eventstudy_ndvi.csv").set_index("event_time").loc[KS, "coef"].to_numpy()
acc_diff = np.abs(b0 - base).max()
print(f"ACCEPTANCE (interactions omitted vs B0): max |diff| = {acc_diff:.2e} -> {'PASS' if acc_diff < 1e-6 else 'FAIL'}")
pd.DataFrame({"event_time": KS, "coef_demeaned_ols": b0, "coef_b0": base}).to_csv(TABLES / "b3_acceptance.csv", index=False)
if acc_diff >= 1e-6:
    raise SystemExit("acceptance failed")

for mode in ("twfe", "dsfe"):
    run("primary", mode=mode)
run("no road control", road=False)
run("EVI", y="evi")
run("excl. 2025 re-flooded pairs", idx=keep_noreflood)
run("excl. 2025 re-flooded pairs", idx=keep_noreflood, mode="dsfe")
for a in ("M_hq_pop50k", "M_hq", "M_hq_osm_major"):
    run(f"access {a}", acc=a, boot=False)

# ---- terciles ------------------------------------------------------------------------------------------
cut = np.quantile(cm.M_hq_osm.to_numpy()[D == 1], [1 / 3, 2 / 3])
terc = 1 + (cm.M_hq_osm.to_numpy() > cut[0]) + (cm.M_hq_osm.to_numpy() > cut[1])
for mode in ("twfe", "dsfe"):
    X, names = design(None, mode, ALL, terciles=terc)
    cmb = {f"b{q}_P": {f"b{q}_{k}": 1 / 3 for k in (6, 7, 8)} for q in (1, 2, 3)} | \
          {f"b{q}_pre": {f"b{q}_-1": 1.0} for q in (1, 2, 3)} | \
          {"b3_P - b1_P": {**{f"b3_{k}": 1 / 3 for k in (6, 7, 8)}, **{f"b1_{k}": -1 / 3 for k in (6, 7, 8)}}}
    res, tab = fit(Y["ndvi"], X, names, ALL, mode, cmb, boot=True)
    tab.insert(0, "spec", f"terciles (cut {cut[0]:.1f}, {cut[1]:.1f} min)"); tab.insert(1, "fe", mode); tab.insert(2, "n_cells", len(cells))
    summary.append(tab); coefs.append(res["beta"].rename(f"terciles|{mode}"))
    print(f"\n== terciles [{mode}]"); print(tab.round(5).to_string(index=False), flush=True)

# ---- negative control (cross-section): preflood trend on D x Mz, district FE, district clusters ---------
import statsmodels.formula.api as smf
nc = pd.DataFrame({"y": M0.preflood_ndvi_trend.to_numpy(), "D": D, "Mz": ACC["M_hq_osm"], "Rz": Rz, "district": dist})
r = smf.ols("y ~ D * Mz + D * Rz + C(district)", nc).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(nc.district)[0]})
negc = pd.DataFrame([{"spec": "negative control: preflood_ndvi_trend", "term": "D:Mz", "estimate": r.params["D:Mz"],
                      "se_district": r.bse["D:Mz"], "p_district": r.pvalues["D:Mz"], "sd_outcome": nc.y.std()}])
print("\n== negative control\n" + negc.round(6).to_string(index=False))

S_ = pd.concat(summary, ignore_index=True)
S_.to_csv(TABLES / "b3_interaction_summary.csv", index=False)
pd.concat(coefs, axis=1).to_csv(TABLES / "b3_coefficients.csv")
negc.to_csv(TABLES / "b3_negative_control.csv", index=False)
print("\nwritten b3_interaction_summary.csv, b3_coefficients.csv, b3_negative_control.csv")
