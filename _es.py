"""Event-study estimator shared by Paper B scripts (ported from flood-recovery-sindh/robustness_v2/_rv2.py, where it
was validated against linearmodels PanelOLS and the published table: coefficients exact, SEs within 0.11 %).
Balanced panel, exact within transform; cluster-robust SEs in the published (linearmodels) convention; wild cluster
restricted bootstrap with Webb 6-point weights.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from _common import SEED

WEBB = np.array([-np.sqrt(1.5), -1, -np.sqrt(.5), np.sqrt(.5), 1, np.sqrt(1.5)])


def demean(Z, grp=None):
    """Exact within transform for a balanced panel: cell + season (grp None) or cell + group x season effects
    (cells nested in grp)."""
    Z = Z - Z.mean(1, keepdims=True)
    if grp is None:
        return Z - Z.mean(0, keepdims=True)
    out = np.empty_like(Z)
    for g in np.unique(grp):
        m = grp == g
        out[m] = Z[m] - Z[m].mean(0, keepdims=True)
    return out


class EventStudy:
    """y_it = a_i + l_t (or l_{g,t}) + sum_k b_k D_i 1[t = k] (+ extra columns), OLS on demeaned data."""

    def __init__(self, Y, D, ets, ref, extra=None, fe_group=None):
        self.ets = list(ets); self.ks = [k for k in self.ets if k != ref]
        n, T = Y.shape
        cols, names = [], []
        for k in self.ks:
            c = np.zeros((n, T)); c[:, self.ets.index(k)] = D; cols.append(c); names.append(f"b_{k}")
        for nm, c in (extra or {}).items():
            cols.append(c); names.append(nm)
        self.names = names
        self.X = np.stack([demean(c, fe_group).ravel() for c in cols], 1)
        self.y = demean(Y, fe_group).ravel()
        self.XtXi = np.linalg.inv(self.X.T @ self.X)
        self.beta = self.XtXi @ (self.X.T @ self.y)
        self.u = self.y - self.X @ self.beta
        self.T = T

    def vcov(self, clusters):
        g = pd.factorize(np.repeat(clusters, self.T))[0]; G = g.max() + 1; N, K = self.X.shape
        S = sp.csr_matrix((np.ones(N), (g, np.arange(N))), shape=(G, N)) @ (self.X * self.u[:, None])
        # published convention (linearmodels PanelOLS clustered, as did_day13.py): equals CR1 x G/(G-1); verified in T0
        return self.XtXi @ (S.T @ S) @ self.XtXi * (G / (G - 1)) ** 2 * ((N - 1) / (N - K)), G

    def test(self, R, clusters, label=""):
        """Linear combination(s) R (q x K): estimate, CR1 SE, t/F test with G-1 df."""
        R = np.atleast_2d(R); V, G = self.vcov(clusters); est = R @ self.beta
        if R.shape[0] == 1:
            se = float(np.sqrt((R @ V @ R.T).item())); t = est[0] / se
            return {"label": label, "estimate": float(est[0]), "se": se, "p": float(2 * stats.t.sf(abs(t), G - 1)), "G": G}
        W = float(est @ np.linalg.solve(R @ V @ R.T, est)); q = R.shape[0]
        return {"label": label, "F": W / q, "p": float(stats.f.sf(W / q, q, G - 1)), "q": q, "G": G}

    def wcb(self, R, clusters, B=9999, seed=SEED):
        """Wild cluster restricted bootstrap p-value and percentile-t 95 % interval for scalar R beta."""
        rng = np.random.default_rng(seed)
        R = np.asarray(R, float).ravel(); X, y = self.X, self.y
        g = pd.factorize(np.repeat(clusters, self.T))[0]; G = g.max() + 1; N, K = X.shape
        c = (G / (G - 1)) ** 2 * ((N - 1) / (N - K))      # same convention as vcov
        Ra = self.XtXi @ R
        Gm = sp.csr_matrix((np.ones(N), (g, np.arange(N))), shape=(G, N))
        XgXg = [X[g == j].T @ X[g == j] for j in range(G)]
        S = Gm @ (X * self.u[:, None])
        se0 = np.sqrt(c * ((S @ Ra) ** 2).sum()); est0 = R @ self.beta; t0 = est0 / se0

        def draws(br, ur):
            Su = Gm @ (X * ur[:, None])
            W = WEBB[rng.integers(0, 6, (B, G))]
            db = (W @ Su) @ self.XtXi.T
            v = np.zeros(B)
            for j in range(G):
                v += ((W[:, [j]] * Su[j][None, :] - db @ XgXg[j]) @ Ra) ** 2
            return db @ R, np.sqrt(c * v)

        a = self.XtXi @ R; br = self.beta - a * est0 / (R @ a); ur = y - X @ br       # restricted (R b = 0)
        num, se = draws(br, ur)
        p = float((np.abs(num / se) >= abs(t0)).mean())
        num_u, se_u = draws(self.beta, self.u)                                         # unrestricted, for interval
        tq = np.percentile(np.abs(num_u / se_u), 95)
        return {"p_wcb": p, "ci_wcb_lo": float(est0 - tq * se0), "ci_wcb_hi": float(est0 + tq * se0), "B": B,
                "weights": "Webb 6-point"}


