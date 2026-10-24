# Paper B: novelty delta against the base flood manuscript

Written 2026-10-23, before any Paper B code. Sources read: `flood-recovery-sindh/docs/results_summary.md`
(Day 20 freeze), `docs/data_dictionary.md`, `src/` file list, `results/tables/` list. The base manuscript text
itself is not in `core-data` (`manuscript/` holds only `references.bib`), so "claims" below are the frozen
results the manuscript is built on. **The researcher should confirm against the manuscript draft.**

## What the base already claims (Paper B must not repeat as findings)

| Base result | Source table | Paper B use |
|---|---|---|
| 555,309 cropland cells (300 m), 22 districts; treated = flood_duration >= 3 (87,948, 15.8 %) | master_dataset | reused input |
| 87,447 matched pairs, caliper 0.2 SD logit, 7 covariates |SMD| < 0.03 | day12_balance_* | reused sample (B0 reproduces) |
| Event study: shock -0.2015 (k=0), smallest gap -0.0576 (k=3, Rabi 2023-24), -0.1167 at k=8 (Kharif 2026); pre-trend p = 0.43 | day13_eventstudy_ndvi | B0 reproduction target only |
| Placebo, caliper, dose-response terciles (-0.131 / -0.193 / -0.280), leave-one-district-out | day14_* | not repeated |
| Kharif-Rabi decoupling (phenology; interaction -0.321) | day16a, day16c_model_b | not repeated; used to define the outcome season |
| TTR mixed model: flood_duration, NDWI persistence, elevation, drainage dominant; **log_dist_road marginal (beta -0.06, p = 0.011); log_pop -0.50** | day16c_model_a_ttr | **overlap risk**: base already tests Euclidean distance to road and population as recovery correlates |
| Moran's I 0.167, SARAR stable | day17_* | not repeated |
| XGBoost + SHAP on TTR (R^2 0.438) | day18a | **overlap risk** for B4 variable importance: B4 must estimate treatment-effect heterogeneity (CATE), not predict TTR |
| Three trajectory archetypes; persistent deficit 17 %, Larkana + Dadu 36 % | day18b, day19 | not repeated; archetype could be an outcome check only |
| Soil-moisture and salinity (SBI) persistence | sm_persistence, master | used as covariates only |

Wording note: the base summary calls -0.0576 the "best recovery point"; the Paper B prompt calls it the
"trough". Paper B will use "smallest gap" to avoid contradiction.

## What Paper B adds

1. **Network market access** (B2): road-network travel time from each cell and village to a defined market-town set
   (plus a gravity measure), on the RHAI road graphs. The base has only Euclidean distance to a road, used as a
   TTR covariate. New: travel time to markets, and its role as a modifier of the flood gap (H6).
2. **Flood x post x access interactions** in the base event study (B3): continuous and tercile access, per-
   interaction pre-trends, placebo, wild cluster bootstrap with 22 clusters, taluka clustering. The base reports
   only the average effect and dose terciles.
3. **Treatment-effect heterogeneity** (B4) with honest causal forest or DML on a persistent-gap outcome, with
   calibration (BLP) tests, against a linear interaction benchmark, leave-district-out. The base SHAP
   explains TTR levels, not heterogeneity in the flood effect.
4. **Water and health access linkage** (SWAI, RHAI / Paper A ensemble) to cropland recovery: new to all three
   earlier projects.
5. **Spatial overlap test of compound deficits** (H3, B5) against spatially aware nulls, with the multidimensional
   method compared with the two-of-three rule.
6. **Indicative income-loss range** (B6), conditional on user-supplied yields and prices.
7. (Proposed, see plan) **Flood-conditioned market access** from Paper A's closure masks, and **2025 re-flooding**
   as a threat to the persistent-gap outcome.

## Risks to novelty

- If market access adds nothing beyond the base's log_dist_road, H6 is a null and must be reported as one.
- Reviewers of the base and Paper B may see salami slicing; the cover letter must state that the matched sample and
  event study come from the base manuscript, and Paper B cites it.
