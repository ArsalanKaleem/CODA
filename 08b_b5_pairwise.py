"""B5 secondary (prereg §4 "each pairwise count"): pairwise overlaps (water-health, water-crop, health-crop) and
">= 2 of 3" against the same two spatial nulls as 08_b5_overlap.py (25 km block permutation; Moran spectral
randomisation), primary layers only. Also the observed / expected ratio and Spearman correlations of the three
continuous layers. Descriptive co-location; S (distance to perennial water, dry-season rain) and health travel
time can both track remoteness, which is part of what any overlap means.
Writes results/tables/b5_pairwise.csv.
"""
import importlib.util

import numpy as np
import pandas as pd

from _common import PAPER, TABLES, SEED, log_env

log_env(__file__)
spec = importlib.util.spec_from_file_location("b5", PAPER / "08_b5_overlap.py")
src = (PAPER / "08_b5_overlap.py").read_text(encoding="utf-8").split("rows = []")[0]   # definitions only, no runs
ns = {"__file__": str(PAPER / "08_b5_overlap.py"), "__name__": "b5defs"}
exec(compile(src, "08_b5_overlap.py", "exec"), ns)
X, deficits, msr_basis = ns["X"], ns["deficits"], ns["msr_basis"]
B = 9999; rng = np.random.default_rng(SEED + 1)
F, Dv = deficits(X)
a = F.to_numpy().astype(int)
PAIRS = {"water-health": (0, 1), "water-crop": (0, 2), "health-crop": (1, 2)}
stat = lambda m: {**{k: int((m[:, i] & m[:, j]).sum()) for k, (i, j) in PAIRS.items()}, "ge2": int((m.sum(1) >= 2).sum())}
obs = stat(a)
xy = X.loc[F.index, ["x", "y"]].to_numpy()
blk = pd.factorize(pd.Series(list(zip((xy[:, 0] // 25000).astype(int), (xy[:, 1] // 25000).astype(int)))))[0]
groups = [np.flatnonzero(blk == g) for g in np.unique(blk) if (blk == g).sum() > 1]
NA = {k: np.empty(B, int) for k in obs}
for b in range(B):
    c = a.copy()
    for j in range(3):
        for gix in groups:
            c[gix, j] = a[rng.permutation(gix), j]
    for k, v in stat(c).items():
        NA[k][b] = v
E = msr_basis(xy); vals = Dv.to_numpy(float); q = [0.2, 0.8, 0.2]; low = [True, False, True]
coef = [E.T @ (vals[:, j] - vals[:, j].mean()) for j in range(3)]
NB = {k: np.empty(B, int) for k in obs}
for b in range(B):
    ind = []
    for j in range(3):
        s = E @ (coef[j] * rng.choice([-1.0, 1.0], len(coef[j]))); thr = np.quantile(s, q[j])
        ind.append(s <= thr if low[j] else s >= thr)
    for k, v in stat(np.column_stack(ind).astype(int)).items():
        NB[k][b] = v
n = len(F); expct = {"water-health": n * .04, "water-crop": n * .04, "health-crop": n * .04, "ge2": n * (3 * .04 * .8 + .008)}
rows = [{"stat": k, "observed": obs[k], "expected_if_independent": expct[k], "obs_over_expected": obs[k] / expct[k],
         "null_block_mean": NA[k].mean(), "null_block_p95": np.percentile(NA[k], 95), "p_block": (1 + (NA[k] >= obs[k]).sum()) / (B + 1),
         "null_msr_mean": NB[k].mean(), "null_msr_p95": np.percentile(NB[k], 95), "p_msr": (1 + (NB[k] >= obs[k]).sum()) / (B + 1)}
        for k in obs]
R = pd.DataFrame(rows)
sp = Dv.rank().corr(method="pearson")
R.to_csv(TABLES / "b5_pairwise.csv", index=False); sp.to_csv(TABLES / "b5_layer_spearman.csv")
pd.set_option("display.width", 220); print(R.round(4).to_string(index=False)); print(sp.round(3).to_string())
