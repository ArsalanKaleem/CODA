"""Paper B figures (static, for the manuscript). Reads cached tables / data only.

Palette: first three categorical slots of the dataviz reference palette (blue #2a78d6, orange #eb6834, aqua
#1baf7a; documented as passing all-pairs CVD and normal-vision checks in light mode), grey for context marks.
Thin marks, recessive grid, legends for >= 2 series, text in ink colours.
  fig1_baseline_event_study.png   base event study (published FE) vs district x season FE (B0, B0b)
  fig2_access_terciles.png        persistent-window gap by market-access tercile (B3), point estimates by season
  fig4_decoupling.png             pairwise deficit overlaps: observed / expected vs two spatial nulls (B5)
  fig5_maps.png                   (a) villages with water+health vs crop-loss deficits; (b) flooded cells cut off
                                  from every market town in the 2022 closures (M_F)
  fig3_heterogeneity.png          (a) GATES by CATE quintile, both halves; (b) leave-district-out BLP slopes (B4)
"""
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _common import DATA, FIGS, FLOOD, PAPER, PAPER_A, TABLES, CID, VID, log_env

log_env(__file__)
C1, C2, C3, GREY, INK, INK2 = "#2a78d6", "#eb6834", "#1baf7a", "#b8b6ad", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": "#e6e5e0",
                     "grid.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2})


def save(fig, name):
    fig.tight_layout(); fig.savefig(FIGS / name, dpi=600); plt.close(fig); print("wrote", name)


# ---- fig 1 --------------------------------------------------------------------------------------------
a = pd.read_csv(TABLES / "b0_eventstudy_ndvi.csv"); b = pd.read_csv(TABLES / "b0b_eventstudy_ndvi_district_season_fe.csv")
fig, ax = plt.subplots(figsize=(6.5, 3))
ax.axhline(0, color=INK2, lw=0.6); ax.axvline(-0.5, color=GREY, lw=0.8, ls="--")
ax.errorbar(a.event_time - 0.1, a.coef, yerr=1.96 * a.se, fmt="o-", color=C1, ms=4, capsize=2, label="Cell + season FE (published)")
b2 = pd.concat([b, pd.DataFrame([{"event_time": -2, "coef": 0, "se": 0}])]).sort_values("event_time")
ax.errorbar(b2.event_time + 0.1, b2.coef, yerr=1.96 * b2.se, fmt="s-", color=C2, ms=4, capsize=2, label="Cell + district x season FE")
ax.set_xlabel("Season relative to Kharif 2022 (k)"); ax.set_ylabel("NDVI anomaly, flooded minus matched control")
ax.legend(frameon=False, loc="lower right"); save(fig, "fig1_baseline_event_study.png")

# ---- fig 2 --------------------------------------------------------------------------------------------
co = pd.read_csv(TABLES / "b3_coefficients.csv", index_col=0)["terciles|twfe"]
ks = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]
fig, ax = plt.subplots(figsize=(6.5, 3))
ax.axhline(0, color=INK2, lw=0.6)
for q, c, lab, m in ((1, C1, "Nearest third (< 16.9 min)", "o"), (2, C3, "Middle third", "^"), (3, C2, "Most remote third (> 28.6 min)", "s")):
    ax.plot(ks, [co[f"b{q}_{k}"] for k in ks], marker=m, color=c, ms=4, label=lab)
ax.set_xlabel("Season relative to Kharif 2022 (k)"); ax.set_ylabel("NDVI gap (point estimates)")
ax.legend(frameon=False, loc="lower right", title="Travel time to nearest market town", title_fontsize=8)
save(fig, "fig2_access_terciles.png")

# ---- fig 3 --------------------------------------------------------------------------------------------
pw = pd.read_csv(TABLES / "b5_pairwise.csv").set_index("stat").loc[["water-health", "water-crop", "health-crop"]]
fig, ax = plt.subplots(figsize=(5.5, 3.3))
y = np.arange(3)[::-1]
ax.barh(y, pw.obs_over_expected, color=GREY, height=0.5, label="Observed")
for col, c, m, lab in (("null_block_mean", C1, "o", "Null: 25 km block permutation (mean)"),
                       ("null_msr_mean", C2, "D", "Null: Moran spectral randomisation (mean)")):
    ax.scatter(pw[col] / pw.expected_if_independent, y, color=c, marker=m, s=36, zorder=3, label=lab)
ax.axvline(1, color=INK2, lw=0.8, ls="--")
ax.set_yticks(y); ax.set_yticklabels(["Water and health", "Water and crop loss", "Health and crop loss"])
ax.set_xlabel("Villages with both deficits / number expected if independent")
for yi, (o, n) in zip(y, zip(pw.observed, pw.expected_if_independent)):
    ax.text(0.03, yi + 0.3, f"{int(o)} villages (vs {n:.0f})", color=INK2, fontsize=7)
ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.45, -0.22), ncol=2); save(fig, "fig4_decoupling.png")

# ---- fig 4 --------------------------------------------------------------------------------------------
src = (PAPER / "08_b5_overlap.py").read_text(encoding="utf-8").split("rows = []")[0]
ns = {"__file__": str(PAPER / "08_b5_overlap.py"), "__name__": "b5defs"}; exec(compile(src, "08_b5_overlap.py", "exec"), ns)
F, _ = ns["deficits"](ns["X"]); X = ns["X"].loc[F.index]
wh = F.water & F.health; cr = F.crop
cells = pd.read_parquet(DATA / "b1_cells.parquet").set_index(CID)
cf = pd.read_parquet(DATA / "b2_cell_market_flood.parquet").set_index(CID)
M = pd.read_parquet(FLOOD / "data/processed/master_dataset.parquet", columns=[CID, "treated_t3", "matched_pair_id"]).set_index(CID)
tr = M[(M.treated_t3 == 1) & M.matched_pair_id.notna()].join(cells[["x", "y"]]).join(cf.isolated_F)
fig, axs = plt.subplots(1, 2, figsize=(7, 4.8))
for ax in axs:
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
ax = axs[0]
ax.scatter(X.x, X.y, s=1, color=GREY, label="Other villages")
ax.scatter(X.x[wh & ~cr], X.y[wh & ~cr], s=4, color=C1, label=f"Water and health deficit ({int((wh & ~cr).sum())})")
ax.scatter(X.x[cr & ~wh], X.y[cr & ~wh], s=4, color=C2, label=f"Crop-loss deficit ({int((cr & ~wh).sum())})")
ax.scatter(X.x[cr & wh], X.y[cr & wh], s=10, color=C3, marker="^", label=f"All three ({int((cr & wh).sum())})")
ax.set_title("(a) Villages by deficit (bottom or top quintile)", fontsize=8, color=INK, loc="left")
ax.legend(frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(0, 0), markerscale=2)
ax = axs[1]
ax.scatter(tr.x[~tr.isolated_F], tr.y[~tr.isolated_F], s=0.2, color=GREY, label="Flooded cell, still linked")
ax.scatter(tr.x[tr.isolated_F], tr.y[tr.isolated_F], s=0.4, color=C2, label=f"Flooded cell cut off from every market town ({tr.isolated_F.mean():.0%})")
ax.set_title("(b) Flooded cropland, 2022 road closures", fontsize=8, color=INK, loc="left")
ax.legend(frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(0, 0), markerscale=12)
save(fig, "fig5_maps.png")

# ---- fig 5 --------------------------------------------------------------------------------------------
g = pd.read_csv(TABLES / "b4_gates_full_t2000.csv"); lo = pd.read_csv(TABLES / "b4_lodo_full_t500.csv").dropna(subset=["blp_slope"])
fig, axs = plt.subplots(1, 2, figsize=(7, 3), gridspec_kw={"width_ratios": [1, 1.2]})
ax = axs[0]; ax.axhline(0, color=INK2, lw=0.6)
for h, c, m in (("A", C1, "o"), ("B", C2, "s")):
    d = g[g.fit_half == h]
    ax.errorbar(d.quintile + (-.08 if h == "A" else .08), d.gate_aipw, yerr=1.96 * d.gate_se, fmt=m, color=c, ms=4, capsize=2,
                label=f"Forest fit on half {h}, evaluated on the other")
ax.set_xlabel("Predicted-effect quintile (1 = most negative)"); ax.set_ylabel("Group effect on persistent NDVI change")
ax.set_title("(a) In-sample heterogeneity", fontsize=8, loc="left"); ax.legend(frameon=False, fontsize=6.5)
ax = axs[1]; lo = lo.sort_values("blp_slope"); yy = np.arange(len(lo))
ax.errorbar(lo.blp_slope, yy, xerr=1.96 * lo.blp_se, fmt="o", color=C1, ms=3, capsize=1.5)
ax.axvline(0, color=INK2, lw=0.6); ax.axvline(1, color=GREY, lw=0.8, ls="--")
ax.set_yticks(yy); ax.set_yticklabels(lo.district, fontsize=6); ax.set_xlim(-8, 8)
ax.set_xlabel("BLP slope in the held-out district (1 = calibrated)"); ax.set_title("(b) Leave-one-district-out", fontsize=8, loc="left")
save(fig, "fig3_heterogeneity.png")
