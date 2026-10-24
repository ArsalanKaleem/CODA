# Paper B pre-registration: H6 (B3), H7 (B4), H3 (B5)

Status: written 2026-10-23, **before** any access measure is computed or any interaction, heterogeneity or
overlap statistic is estimated. Only B0 (reproduction of the base event study, no new estimand) may run before
freezing. The frozen copy is `docs/prereg_B_<date>.md`; its SHA-256 goes into `docs/decision_log.md`. Any later
change is a logged deviation, labelled in the manuscript.

Defaults agreed by the researcher on 2026-10-23 (plan.md questions 1-4, 6, 7): water access = static SWAI
components; persistent window k = 6 to 8 with a Kharif-only variant; market towns = district and taluka HQ plus
places with population >= 50,000; Sindh-wide major-road graph; gravity decay from a sourced literature range.

Design caveat that applies throughout: market, water and health access are **not randomly assigned**. All
modifier results are effect modification (where the flood gap is larger or smaller), never the causal effect
of access. No causal language about modifiers in text or comments.

## 1. Common definitions

- **Sample**: the base matched sample, `did_panel.parquet` (87,447 pairs, 174,894 cells, 11 seasons,
  k = -2 to 8, reference k = -2 = Kharif 2021). Treated = `treated_t3`. No re-matching.
- **Persistent window** P = {6, 7, 8} (Kharif 2025, Rabi 2025-26, Kharif 2026). Kharif-only variant
  P_K = {6, 8}. Pre-period = {-2, -1}.
- **Market access M** (primary): road-network travel time, minutes, from the cell centroid (or village point) to the
  nearest market town. Towns = GeoNames feature codes PPLA, PPLA2, PPLA3 in Sindh (provincial, district, taluka
  HQ) plus GeoNames places with population >= 50,000. Graph: Sindh-wide major roads (Geofabrik), RHAI speed
  assumptions; off-network leg at the RHAI convention (to be checked, PB7). Sensitivities: thresholds 20,000 and
  100,000; gravity access G = sum_j pop_j exp(-beta t_ij) over a sourced beta range (PB2). Fixed after B2 passes
  its spot-check, before B3.
- **Flood-conditioned market access** M_F (secondary): same, with roads closed under Paper A's 2022 base rule
  (v1 mask, D3, T0.6, bridges open). Delta M = M_F - M.
- **Water access W** (primary): S = max(R, W_dry) from `SWAI/data/swai_final.csv` (distance to perennial
  water and dry-season rainfall, static or pre-flood). Full SWAI score is a secondary variant (its N component is a
  2023 post-flood NDWI composite and E rewards low elevation; both track flooding).
- **Health access H** (primary): RHAI Scenario B `min_any` (minutes). Variants: Scenario A; Paper A ensemble
  median (`paper_a/data/a4_village_distributions.parquet`).
- **Cell <-> village link**: B3/B4 use cell-level M. W and H are cell-assigned from the nearest village within
  5 km (cells beyond 5 km: missing, excluded from models using W or H; count reported). B5 uses villages with
  >= 10 cropland cells within 2 km (sensitivity 5 and 25; PB6).
- **2025 re-flooding**: cell flagged if the Paper A 2025 v1 mask has >= 3 detections over >= 50 % of its 30 m
  pixels (the >= 3 threshold follows the base treatment; the 50 % pixel-share aggregation is a Paper B choice,
  because the base's 30 m to 300 m aggregation rule has not been verified yet). Primary analyses keep all cells;
  robustness excludes pairs with either member reflooded.

## 2. H6 (B3): market access and the flood gap

Model (cell i, season t, event time k):

y_it = a_i + l_t + sum_{k != -2} b_k D_i 1[k] + sum_{k != -2} g_k D_i 1[k] Mz_i + sum_{k != -2} d_k 1[k] Mz_i
       + sum_{k != -2} r_k D_i 1[k] Rz_i + sum_{k != -2} s_k 1[k] Rz_i + e_it

y = `ndvi_anom` (EVI secondary). D = treated. Mz = log1p(M) standardised over the matched sample. Rz = log1p
(`dist_to_road_m`) standardised (the base paper's road covariate, so M is tested beyond it).

- **Primary estimand** g_P = mean of g_k over k in P. **H6 supported** if g_P < 0 (longer travel time, larger
  persistent deficit) with wild cluster bootstrap p < 0.05 **and** the pre-trend coefficient g_-1 has p >= 0.10.
  A positive or null g_P is reported as such.
- **Inference**: CR1 clustered by district (22); **primary p-value = wild cluster restricted bootstrap**, Webb
  6-point weights, 9,999 draws, seed 20221001; robustness clustering by taluka (GADM L3) and by matched pair.
- **Terciles**: M terciles from the treated cells' distribution; event study per tercile (same FE).
- **Checks**: B0 coefficients unchanged when interactions are omitted (within 1e-6); negative control
  `preflood_ndvi_trend` regressed on D x Mz (cross-section, district FE), expected null; Kharif-only g_PK; 2025
  re-flood exclusion; M thresholds and G; M_F and Delta M in place of M (secondary).
- Multiple outcomes: NDVI primary; everything else is labelled secondary or robustness, no correction applied to them.

## 3. H7 (B4): heterogeneity of the flood effect (exploratory)

- **Outcome** Y_i = mean `ndvi_anom` over P minus mean over the pre-period (cell-level change).
- **Treatment** D_i; **modifiers X**: Mz, W (S), H, elevation, `sbi_persist_months`, `sm_persist_months`,
  `flood_duration` is excluded from X (it is zero for controls by construction). **Controls W_ctl**: the seven
  matching covariates + district indicators.
- **Estimator**: EconML `CausalForestDML`, honest, cross-fitting 5 folds grouped by matched pair, nuisance models
  gradient boosting, 2,000 trees (reduced to 500 if runtime > 30 min, logged), seed 20221001. Cells missing W or
  H dropped (count reported).
- **Heterogeneity test** (decides whether to interpret CATEs): BLP calibration on held-out folds (Chernozhukov et al.
  style; `econml.validate.DRTester` or equivalent): **real heterogeneity** if the BLP slope on the CATE is > 0
  with p < 0.05. Also GATES by CATE quintile and RATE/AUTOC. If the test fails: CATEs not interpreted; the
  researcher is asked whether to reframe.
- **Overlap**: propensity of D given X and W_ctl, share of cells outside [0.05, 0.95] reported, trimmed sensitivity.
- **Benchmark**: OLS of Y on D, X, D x X, W_ctl, district-clustered; compare sign and rank of modifiers.
- **Predictive claims**: leave-one-district-out (fit on 21 districts, BLP on the held-out one); only effects stable
  across folds are described.
- Variable importance is reported with the caveat that modifiers are correlated and non-random.

## 4. H3 (B5): overlap of deficits

- **Unit**: villages meeting the cell rule (§1).
- **Deficits** (indicator = 1):
  - water: bottom quintile of S;
  - health: top quintile of H (Scenario B);
  - crop: top quintile of village persistent crop loss L = mean over cropland cells within 2 km of
    [mean integrated NDVI in Kharif 2025 and 2026] minus [integrated NDVI in Kharif 2021]
    (`phenology.parquet`; all cells, not only matched, so every village has a value; descriptive, not flood-attributable).
- **Statistic**: count of villages with all three deficits (T3); secondary: each pairwise count and the number with
  >= 2.
- **Nulls** (9,999 replicates each): (a) permutation of each layer within 25 km x 25 km blocks (keeps block-level
  structure); (b) Moran spectral randomisation surrogates of each continuous layer (Wagner and Dray 2015),
  re-thresholded at the same quintiles, preserving each layer's spatial autocorrelation. **H3 supported** if T3
  exceeds the 95th percentile of **both** nulls. Reported either way.
- **Variants**: full SWAI; health A and Paper A ensemble (shows whether open facility data change who is flagged);
  flood-attributable crop loss in flooded villages (mean pair-differenced gap over P for villages with >= 10
  treated matched cells).
- **Rule vs multidimensional**: rule = bottom quintile on >= 2 of 3. Multidimensional = Gaussian mixture on the
  three standardised ranks, k = 2 to 8 by BIC, stability by bootstrap ARI (200 resamples); profiles with mean ARI
  < 0.6 are reported as unstable, k is not presupposed. Villages flagged by one method and not the other are
  compared on an **independent outcome**: Paper A newly isolated under the 2022 flood (1 km rule), which uses
  none of the three deficit inputs.

## 5. Not pre-registered

B6 (income-loss range) is a conversion with user-supplied yields and prices, not a test. Any analysis not
listed here is labelled post hoc.
