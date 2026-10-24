"""Audit: every number cited in docs/manuscript_skeleton.md, read from the cached tables / data (no new estimation).
Prints KEY = value and writes results/paper_b_numbers.txt / .json. The skeleton cites numbers by these keys.
Baseline-audit numbers come from flood-recovery-sindh/robustness_v2/results/robustness_v2_numbers.json.
"""
import json

import numpy as np
import pandas as pd

from _common import DATA, FLOOD, PAPER, TABLES, CID, log_env

log_env(__file__)
t = lambda n: pd.read_csv(TABLES / n)
N = {}
# B0 / B0b
b0 = t("b0_reproduction_check.csv"); N["b0.pass"] = bool(b0["pass"].all())
e = t("b0_eventstudy_ndvi.csv").set_index("event_time")
for k in (-1, 0, 3, 8):
    N[f"b0.k{k}"] = float(e.loc[k, "coef"])
N["b0.mean6_8"] = float(e.loc[[6, 7, 8], "coef"].mean())
bb = t("b0b_eventstudy_ndvi_district_season_fe.csv").set_index("event_time")
for k in (-1, 0, 8):
    N[f"b0b.k{k}"] = float(bb.loc[k, "coef"]); N[f"b0b.k{k}_p"] = float(bb.loc[k, "pvalue"])
# B1
cov = t("b1_coverage_by_district.csv").set_index("district").loc["TOTAL"]
N["b1.villages"] = int(cov.villages); N["b1.villages_ge10_cells"] = int(cov.cells_ge10); N["b1.villages_ge10_treated"] = int(cov.treated_ge10)
nv = pd.read_parquet(DATA / "b1_cell_nearest_village.parquet"); N["b1.cells_village_le5km"] = float((nv.dist_m <= 5000).mean())
# B2
cm = pd.read_parquet(DATA / "b2_cell_market.parquet"); vmk = pd.read_parquet(DATA / "b2_village_market.parquet")
N["b2.cell_M_median"] = float(cm.M_hq_osm.median()); N["b2.cell_M_p90"] = float(cm.M_hq_osm.quantile(.9))
N["b2.village_M_median"] = float(vmk.M_hq_osm.median()); N["b2.cell_snap_median_m"] = float(cm.snap_m.median())
towns = pd.read_parquet(DATA / "b2_towns.parquet"); N["b2.n_towns_primary"] = int(towns.in_hq_osm.sum()); N["b2.n_hq"] = int(towns.hq.sum())
sc = t("b2_spotcheck_osrm.csv"); N["b2.spot_ratio_min"] = float(sc.ratio_full_osrm.min()); N["b2.spot_ratio_max"] = float(sc.ratio_full_osrm.max())
fa = t("b2_flood_access_summary.csv").set_index("group")
N["mf.iso_treated"] = float(fa.loc["matched treated", "isolated_share"]); N["mf.iso_control"] = float(fa.loc["matched controls", "isolated_share"])
N["mf.dM_gt10_treated"] = float(fa.loc["matched treated", "dM_gt10min_share"])
cs = pd.read_parquet(DATA / "b2_flood_closed_segs.parquet"); N["mf.closed_piece_share"] = float(cs.closed.mean())
# re-flood flags
rf = t("b3_reflood_check.csv").set_index("check")["value"]
N["rf.agree_2022"] = float(rf["2022 rule vs treated_t3, agreement (all cells)"])
N["rf.treated_2025"] = float(rf["2025 flagged share, matched treated"]); N["rf.control_2025"] = float(rf["2025 flagged share, matched controls"])
# B3
s3 = t("b3_interaction_summary.csv")
def g3(spec, fe, combo, col):
    return float(s3[(s3.spec == spec) & (s3.fe == fe) & (s3.combo == combo)][col].iloc[0])
for spec, tag in (("primary", "prim"), ("no road control", "noroad"), ("excl. 2025 re-flooded pairs", "noreflood"), ("EVI", "evi")):
    for fe in ("twfe", "dsfe"):
        if ((s3.spec == spec) & (s3.fe == fe)).any():
            N[f"b3.{tag}.{fe}.gP"] = g3(spec, fe, "g_P (mean k6-8)", "estimate")
            N[f"b3.{tag}.{fe}.gP_pwcb"] = g3(spec, fe, "g_P (mean k6-8)", "p_wcb_district")
            N[f"b3.{tag}.{fe}.gpre_pwcb"] = g3(spec, fe, "g_pre (k-1)", "p_wcb_district")
N["b3.prim.twfe.gP_se"] = g3("primary", "twfe", "g_P (mean k6-8)", "se_district")
for fe in ("twfe", "dsfe"):
    tq = s3[(s3.spec.str.startswith("terciles")) & (s3.fe == fe)].set_index("combo")
    for c in ("b1_P", "b2_P", "b3_P", "b3_P - b1_P"):
        N[f"b3.terc.{fe}.{c.replace(' ', '')}"] = float(tq.loc[c, "estimate"])
    N[f"b3.terc.{fe}.diff_pwcb"] = float(tq.loc["b3_P - b1_P", "p_wcb_district"])
nc = t("b3_negative_control.csv"); N["b3.negctrl"] = float(nc.estimate.iloc[0]); N["b3.negctrl_p"] = float(nc.p_district.iloc[0])
# B4
bl = t("b4_blp_full_t2000.csv")
N["b4.n_fit_total"] = int(bl.n_fit.sum()); N["b4.blp_slope"] = [float(x) for x in bl.blp_slope]; N["b4.blp_p"] = [float(x) for x in bl.blp_slope_p_one_sided]
N["b4.ate"] = [float(x) for x in bl.ate_const]; N["b4.prop_outside"] = [float(x) for x in bl.prop_outside_05_95]
gt = t("b4_gates_full_t2000.csv"); N["b4.gate_q1"] = float(gt[gt.quintile == 1].gate_aipw.mean()); N["b4.gate_q5"] = float(gt[gt.quintile == 5].gate_aipw.mean())
lo = t("b4_lodo_full_t500.csv").dropna(subset=["blp_slope"]); w = 1 / lo.blp_se ** 2
N["b4.lodo_n"] = int(len(lo)); N["b4.lodo_pooled"] = float((w * lo.blp_slope).sum() / w.sum()); N["b4.lodo_pooled_se"] = float(np.sqrt(1 / w.sum()))
N["b4.lodo_pos_sig"] = int(((lo.blp_slope > 0) & (lo.blp_p < .05)).sum()); N["b4.lodo_neg_sig"] = int(((lo.blp_slope < 0) & (lo.blp_p < .05)).sum())
ols = t("b4_ols_benchmark.csv").set_index("Unnamed: 0")
N["b4.ols_DxMz"] = float(ols.loc["D_x_Mz", "coef"]); N["b4.ols_DxMz_p"] = float(ols.loc["D_x_Mz", "p"])
M = pd.read_parquet(FLOOD / "data/processed/master_dataset.parquet", columns=["district", "treated_t3", "matched_pair_id"])
mm = M[M.matched_pair_id.notna()]; N["b4.mithi_treated"] = int(mm[mm.district == "Mithi"].treated_t3.sum()); N["b4.mithi_cells"] = int((mm.district == "Mithi").sum())
# B5
ov = t("b5_overlap_tests.csv").set_index("variant"); p = ov.loc["primary (S, RHAI B, L)"]
N["b5.n"] = int(p.n_villages); N["b5.T3"] = int(p.T3); N["b5.T3_expected"] = float(p.expected_if_independent)
N["b5.p_block"] = float(p.p_block); N["b5.p_msr"] = float(p.p_msr); N["b5.supported_any"] = bool(ov.H3_supported.any())
pw = t("b5_pairwise.csv").set_index("stat")
for s in ("water-health", "water-crop", "health-crop"):
    N[f"b5.{s}.obs"] = int(pw.loc[s, "observed"]); N[f"b5.{s}.ratio"] = float(pw.loc[s, "obs_over_expected"])
    N[f"b5.{s}.p_block"] = float(pw.loc[s, "p_block"]); N[f"b5.{s}.p_msr"] = float(pw.loc[s, "p_msr"])
N["b5.expected_pair"] = float(pw.loc["water-health", "expected_if_independent"])
fv = ov.loc["crop = flood-attributable gap (villages >= 10 treated cells)"]; N["b5.flooded.HC"] = int(fv.HC); N["b5.flooded.n"] = int(fv.n_villages)
N["b5.flooded.HC_expected"] = float(fv.n_villages * 0.04)
sp = pd.read_csv(TABLES / "b5_layer_spearman.csv", index_col=0); N["b5.rho_H_L"] = float(sp.loc["H_B", "L"])
st = pd.read_csv(TABLES / "b5_gmm_stability.csv", index_col=0); N["b5.ari_mean"] = float(st.loc["mean", "ari"])
bic = t("b5_gmm_bic.csv"); N["b5.k_bic"] = int(bic.loc[bic.bic.idxmin(), "k"])
rv = t("b5_rule_vs_multi.csv").set_index("group"); N["b5.multi_only"] = int(rv.loc["multidimensional only", "n"])
# H8
h8 = t("b8_flood_cut_access.csv")
def g8(spec, fe, combo, col):
    return float(h8[(h8.spec == spec) & (h8.fe == fe) & (h8.combo == combo)][col].iloc[0])
N["h8.hP"] = g8("primary: Iso_F", "twfe", "h_P (mean k6-8)", "estimate"); N["h8.hP_pwcb"] = g8("primary: Iso_F", "twfe", "h_P (mean k6-8)", "p_wcb")
N["h8.hpre"] = g8("primary: Iso_F", "twfe", "h_pre (k-1)", "estimate"); N["h8.hpre_pwcb"] = g8("primary: Iso_F", "twfe", "h_pre (k-1)", "p_wcb")
N["h8.dM.hP"] = g8("secondary: dMz (added minutes)", "twfe", "h_P (mean k6-8)", "estimate"); N["h8.dM.hP_pwcb"] = g8("secondary: dMz (added minutes)", "twfe", "h_P (mean k6-8)", "p_wcb")
# baseline audit (robustness_v2)
rv2 = json.loads((FLOOD / "robustness_v2/results/robustness_v2_numbers.json").read_text())
for k in ("range.persistent.all", "range.shock.all", "t4.a.v2.thr0.5.share_of_persistent_gap_linked", "t4.a.v1.thr0.5.share_of_persistent_gap_linked",
          "t2c.ndvi_pub_twfe.mean6_8.breakdown_Mbar", "t2c.ndvi_pub_dsfe.mean6_8.breakdown_Mbar", "t1.ext.joint_p", "t1.dsfe.lead_m1_pwcb",
          "t3.published.pair_km_50"):
    N[f"audit.{k}"] = rv2.get(k)
(PAPER / "results/paper_b_numbers.json").write_text(json.dumps(N, indent=1, default=str))
(PAPER / "results/paper_b_numbers.txt").write_text("\n".join(f"{k} = {v}" for k, v in N.items()))
print("\n".join(f"{k} = {v}" for k, v in N.items()))
