# Paper B decision log

Chronological. Every decision, check and result, with date.

## 2026-10-23

- Plan, novelty delta and parameters list written (`docs/plan.md`, `novelty_delta.md`, `parameters_to_source.md`).
- Found while planning: SWAI component N is a 2023 (post-flood) NDWI composite and E rewards low elevation, so the
  composite score partly tracks flooding. The base TTR model already includes `log_dist_road` (marginal) and
  `log_pop`.
- Researcher accepted the plan defaults:
  - Q1: water access = S = max(R, W), full SWAI as a variant.
  - Q2: persistent window k = 6 to 8, with a Kharif-only variant.
  - Q3: towns = district and taluka HQ plus places with population >= 50,000 (20,000 / 100,000 as sensitivity).
  - Q4: Sindh-wide major-road graph.
  - Q6: the base manuscript text is still to be located.
  - Q7: gravity beta taken from a sourced literature range.
  - B6 yields and prices to be supplied later.
- Environment: `paper_b/.venv`, Python 3.12.10 (EconML wheels), pinned in `requirements.txt`; seed 20221001.

## 2026-10-24

- **B0 PASS.** `00_b0_reproduce_base.py` reproduces `day13_eventstudy_{ndvi,evi}.csv` exactly: all 22
  coefficients and SEs agree to 5 decimals (abs diff 0). Panel: 1,923,834 rows, 174,894 cells, 87,447 pairs, 0
  pairs failing integrity, 11 seasons per cell. The saved panel is already winsorised (NDVI anomaly clipped at
  +-0.911). Shock k=0 -0.2015 (SE 0.0096); smallest gap k=3 -0.0576; k=8 -0.1167; pre-trend k=-1 +0.0033 (p 0.43).
- **Pre-registration frozen** before any access measure or interaction was computed:
  `docs/prereg_B_2026-10-24.md`, SHA-256 `b4713636eb561617062dc205112527e883fa599b3b6ed1e75b1e065ac73dcd5f`.
  Covers H6 (B3), H7 (B4), H3 (B5). B6 is not pre-registered.
- Open items: B6 yields and prices (researcher); base manuscript draft location; PB2 beta range and PB7
  off-network speed still to source; Geofabrik roads path (`GEOFABRIK_ROADS`) needed for B2.
- **B1 done** (`01_b1_linkage.py`):
  - 30.7 % of the 555,309 cells lie within 2 km of a village (497,044 links). The nearest village is within 5 km for 58.4 %, so cell-level W and H are missing for about 42 % of cells (to be quantified on the matched sample in B3).
  - Of the 5,159 villages, 4,871 have >= 10 cropland cells within 2 km and 2,003 have >= 10 treated matched cells (by district: `results/tables/b1_coverage_by_district.csv`).
  - Village table: `data/b1_village_recovery.parquet`.
  - Descriptive crop change L (Kharif 2025-26 vs 2021, integrated NDVI) averages -0.187 across villages: cropland greenness fell broadly, not only in flooded cells.
- **B2 graph** (`02_b2_graph.py`): the Geofabrik major roads from flood-recovery-sindh give 22,858 km of line.
  - A weakly connected component left 764 Ghotki cells snapped to one-way dead ends, so the graph now uses the largest strongly connected component: 15,222 nodes, 93.6 % kept.
  - Every cell and village reaches a town.
- **B2 towns** (`03_b2_market_access.py`):
  - **GeoNames has no PPLA3 (taluka HQ) for Sindh**, and GADM 4.1 has no taluka level (its level 3 is districts). The pre-registered "district and taluka HQ" therefore reduces to the 21 district HQ (PPLA/PPLA2) plus population thresholds: 52 towns at >= 50k, 75 at 20k, 36 at 100k.
  - Hub (Balochistan) was dropped by requiring GeoNames admin1 = 05.
  - The pre-registered alternative clustering "taluka (GADM L3)" is also not available locally. Both need a researcher decision (deviation).
- **B2 market time**:
  - Median cell M = 29.7 min (50k set), 99th percentile 163 min. Mithi has the longest district median, 138 min.
  - Gravity is computed only for PLACEHOLDER half-lives 30 / 60 / 120 min (PB2 unsourced).
- **B2 spot-check** against the OSRM public demo (2026-10-24; all OSM roads, current), three long journeys, ratio of ours to OSRM: Jamshoro to Kotri 236.8 / 182.6 min (1.30); Mithi to Umarkot 201.4 / 227.5 (0.89); Thatta 104.3 / 79.6 (1.31). The major-roads-only graph overstates times by up to about 30 % where minor roads give shortcuts (`results/tables/b2_spotcheck_osrm.csv`).
- **Researcher decisions (2026-10-24)**:
  - town set: district HQ plus OSM towns;
  - clustering: OCHA tehsil boundaries;
  - spot-check: build a hybrid graph.
- **Hybrid graph** (`02b_b2_graph_full.py`, PRIMARY): built from all Geofabrik 2026-07-12 drive roads in Sindh + 10 km (285,741 lines, 143,323 km).
  - Nodes are shared line vertices, so over- and underpasses are not joined.
  - Largest strongly connected component: 2,482,606 nodes (98.6 %).
  - The 2026 network postdates the flood, which is a stated limitation.
  - The major-roads graph (`02`) is kept as a sensitivity (`*_major.parquet`).
- **Towns (deviation D-B1)**: primary `hq_osm` = 21 GeoNames district HQ + OSM place=town/city in Sindh, 229 points in total (OSM 173 towns, 35 cities; overlaps kept, no effect on minimum times). Max town snap 394 m.
  - Sensitivity sets: `hq` (21) and `hq_pop20k / 50k / 100k` (75 / 52 / 36).
  - Material choice: cell rank correlation of hq_osm vs hq_pop50k is 0.63 (Spearman).
- **Market access (full graph)**: cells M_hq_osm median 22.2 min (90 % 46.8, 99 % 89.7, max 299.8); villages median 16.0 min. Median cell snap 350 m (untimed). Nothing unreachable.
  - Village rank correlation full vs major: 0.82 (hq_osm), 0.92 (hq_pop50k).
  - Gravity G is still a PLACEHOLDER (PB2).
- **Spot-check, full graph vs OSRM**, same three trips: Jamshoro to Kotri 205.1 / 182.6 (1.12); Mithi to Umarkot 197.4 / 227.5 (0.87); Thatta 93.5 / 79.6 (1.17). This is within +-17 %, against up to +31 % on the major-roads graph. **Accepted** (B2 acceptance check).
- **Talukas (deviation D-B2)**: GADM has no taluka level. The OSM admin_level 7 polygons in the same Geofabrik extract give 124 units over 143,051 km2, with no overlaps (`04_b2_talukas.py`).
  - Coverage: 99.92 % of cells assigned (471 missing); matched cells fall in 108 talukas.
  - Used pending the researcher's confirmation, since they chose OCHA tehsils; the swap is a one-file change.
- **Matched pairs are not spatially local**: median pair distance is 202 km (10 % 54 km, 90 % 384 km), and only 7.3 % of pairs share a district (98.4 % span two talukas). So pair-differenced village gaps compare flooded cells with distant controls.
- **B0b (POST HOC)**: the base event study with district x season effects instead of season effects (`00b_b0_district_season_fe.py`). All 22 districts contain both treated and control cells.
  - The post-flood path holds: k=0 -0.183, k=3 -0.054, k=6..8 -0.123 / -0.105 / -0.130.
  - **The pre-trend fails**: k=-1 +0.0132 (SE 0.0046, p = 0.004), against base +0.0033 (p = 0.43). Within districts, treated cells were already greening relative to controls before the flood (about 7 % of the shock). This affects the base manuscript's parallel-trends claim and needs the researcher's attention.
- **2025 re-flooding** (`05a_b3_reflood_2025.py`): the prereg rule (>= 50 % of observed 30 m pixels with >= 3 detections) applied to 2022 agrees with the base `treated_t3` for 97.3 % of cells, which validates the aggregation.
  - In 2025, **18.1 % of matched treated cells but 0.9 % of controls** are flagged; 18.8 % of pairs have either member flagged.
- **B3 results** (`05_b3_interactions.py`, log `logs/05_b3_interactions.log`, tables `b3_*.csv`). Acceptance check: with interactions omitted, the coefficients equal B0 (max diff 5.6e-13), PASS.
  - **H6 not supported.** The pre-registered primary g_P (NDVI, k 6-8, per SD of log travel time to the nearest HQ/OSM town, with road-distance controls) = **+0.020** (CR1 district SE 0.008; **wild cluster bootstrap p = 0.046**; taluka p 0.045). The pre-trend g_-1 = -0.001 (WCB p 0.28).
    - The sign is the **opposite** of H6: longer travel time goes with a slightly *smaller* persistent deficit (the base gap is -0.116 at mean access).
  - **The positive interaction is fragile:**
    - district x season FE (post hoc): +0.011, WCB p 0.18;
    - without the road-distance control: +0.008, p 0.29;
    - excluding 2025 re-flooded pairs: +0.016, p 0.16 (district x season: +0.007, p 0.49).
  - **The interaction depends on the town set:**
    - 50k set: +0.037, CR1 p 0.004, but its pre-trend fails (-0.008, p 0.001);
    - HQ only: +0.017, p 0.15;
    - major-roads graph: +0.019, p 0.015.
  - EVI: +0.013, WCB p 0.039.
  - Kharif-only g_PK +0.022 (p 0.029); Rabi g_R +0.018 (p 0.08).
  - **Terciles** (cuts 16.9 and 28.6 min):
    - persistent gap -0.128 (near), -0.121, -0.102 (remote);
    - remote minus near +0.026, WCB p 0.21;
    - under district x season FE: -0.120 / -0.123 / -0.112, difference p 0.68, and the tercile pre-trends fail (+0.011 to +0.014, p 0.002 to 0.047), as in B0b.
  - **Negative control** (preflood NDVI trend on D x Mz): -0.00006 (p 0.67), null as expected.
  - Matched-pair clustering gives tiny SEs (about 0.001) that ignore spatial correlation; they are reported but not used.
  - **2025 re-flooding matters for the base**: excluding re-flooded pairs shrinks the persistent gap b_P from -0.116 to -0.099 (twfe) / -0.100 (dsfe), so about 15 % of the persistent gap is associated with pairs flooded again in 2025.
  - Not yet run: flood-conditioned access M_F (secondary; needs segment exposure on the full graph).
- **M_F, flood-conditioned market access** (`07_b2_flood_access.py`; 2022, v1 mask, D3/T0.6, bridges open).
  - Exposure is aggregated over junction-to-junction chains, to match Paper A's per-edge rule on the vertex-pair graph: 687,500 chains.
  - 2.93 % of pieces are closed (3,883 of 139,806 km).
  - **19.9 % of matched treated cells lose every link to a market town** (1 km rule), against 1.8 % of controls; 10.3 % of treated cells face more than 10 extra minutes.
  - Not yet used in B3 (secondary spec pending).
- **B4 timing**: 10k cells and 200 trees take 5 s per half; the full run takes 55-68 s per half at 2,000 trees. LODO uses 500 trees (prereg rule: reduce when > 30 min; 2,000 trees would take about 45 min).
- **B4 sample**: 68,754 matched cells (39.3 %, 34,377 complete pairs) have W and H, i.e. a village within 5 km. This is a strong selection, towards cells near OSM villages.
- **B4 results** (`06_b4_causal_forest.py`, full_t2000).
  - Held-out ATE -0.065 / -0.067 (SE 0.006).
  - **BLP slope 1.02 / 1.07 (one-sided p 0.0002 / 0.0006)**, which meets the prereg "real heterogeneity" criterion. Overlap-trimmed (e in [0.05, 0.95], about 24k cells per half): slope 1.52 / 1.68.
  - **Overlap is poor**: 30-31 % of cross-fitted propensities fall outside [0.05, 0.95] (district indicators plus cross-district pairs).
  - GATES, monotone: Q1 -0.124 / -0.129 down to Q5 -0.008 / -0.008.
  - Forest importance: elevation 0.48 / 0.55, Mz 0.21 / 0.18, water S 0.14 / 0.09, SBI 0.10 / 0.12, health 0.06 / 0.04, soil moisture 0.02.
  - OLS benchmark: D x Mz -0.004 (p 0.46), water and health null; D x SBI +0.023 (p 0.044) and D x soil-moisture persistence -0.025 (p 0.002), both post-treatment (descriptive).
- **B4 CLAN** (`06b_b4_clan.py`): the most negative CATE quintile has the fewest treated cells (29 % vs 66 % in Q5), the lowest flood duration (2.4 vs 5.6) and the most 2025 re-flooding (18.8 % vs 2.7 %).
  - 29 % of that quintile is Mithi, where **only 57 of 9,691 matched cells are treated**, so its CATEs are extrapolations into a region without flooded cells.
  - Spearman correlation of CATE with travel time: -0.23.
- **B4 leave-one-district-out** (500 trees; 19 districts evaluable; Malir, Mithi and Tando Allahyar have too few treated cells):
  - BLP slope > 0 with p < 0.05 in 2 districts (Dadu 0.56, Matiari 3.61), < 0 with p < 0.05 in 1 (Umerkot -3.29);
  - inverse-variance pooled slope **0.11 (SE 0.057)**.
  - The CATE ranking does not transfer to unseen districts.
  - **Reading**: the in-sample heterogeneity is largely between-region geography and poor-overlap extrapolation, not a transferable modifier pattern. Prereg stop condition reached; the researcher is asked whether to reframe H7.

## 2026-10-24 (later): researcher asleep, work delegated

The researcher said: "research and choose everything yourself ... continue doing the work as you think will be most novel and concrete". Every choice below is mine, logged before the corresponding run.

- **B5 prereg clarifications**, fixed before running. Each fills a gap the frozen prereg left open and is listed as clarification D-B5 in `deviations_B.md`:
  1. "Top quintile of persistent crop loss" = bottom quintile of L (L < 0 = loss).
  2. Block permutation uses 25 km x 25 km squares in EPSG:32642, anchored at (0, 0).
  3. MSR spatial weights: 8 nearest neighbours, row-standardised then symmetrised; the doubly-centred matrix gives n - 1 eigenvectors. Surrogate = mean + E (c * s), with c the projection of the centred variable and s random signs. This keeps each layer's Moran spectrum (Wagner and Dray 2015, singleton method).
  4. The "multidimensional" flag = members of the Gaussian-mixture component with the highest mean deficit percentile over the three layers (BIC-chosen k).
  5. Independent outcome = Paper A newly isolated under the 2022 flood: scenario B, D3, T0.6, v1 mask, 1 km rule (`min_any_r1000` missing under flood, present normally, road within 1 km).
  6. The Paper A ensemble variant of health access = `a4_village_distributions` config "normal", `p50_min`.
  7. Health access of unreachable villages (inf / missing min_any) is set to the worst value, i.e. counted as a deficit.
- **Amendment 1 (H8, flood-cut market access) frozen** before any estimate with M_F: `docs/prereg_B_amendment1_2026-10-24.md`, SHA-256 `80917fa5691850377078f31db28826c2866c4445c44547fb63a86e6f38efa64d`.
- **H8 result** (`09_b8_flood_cut_access.py`, table `b8_flood_cut_access.csv`): **not supported.**
  - Iso_F covers 19.9 % of treated and 1.8 % of control cells; its correlation with flood duration among treated cells is 0.21.
  - Primary h_P = **+0.072** (SE 0.022, WCB p 0.050). The sign is the OPPOSITE of H8: cut-off cells show a *smaller* persistent gap.
  - **Pre-trend fails**: h_-1 -0.014 (WCB p 0.0045). By the amendment's own rule, h_P is not interpretable.
  - Variants:
    - district x season FE: +0.074 (WCB p 0.090), pre -0.021 (p 0.0065);
    - excluding 2025 re-flooded pairs: +0.055 (WCB p 0.025), pre -0.012 (p 0.023);
    - without flood-duration / access / road controls: +0.076, pre -0.012.
  - Secondary dMz (added minutes, among cells not cut off): -0.007 (WCB p 0.82), null.
  - Reading: cut-off cells were already declining relative to controls before the flood (negative lead), then rose relative to them after. This pattern fits mean reversion or post-flood non-crop greening on waterlogged land (the "apparent greening" of TT-ORU 2026) better than a market-access mechanism. Untested speculation.
- **B5 / H3 results** (`08_b5_overlap.py`, `08b_b5_pairwise.py`; 4,871 villages with >= 10 cells): **H3 not supported in any variant.**
  - Triple overlap T3 (expected 39 if the layers were independent):
    - primary 38 (block null mean 42, p 0.88; MSR null mean 47, p 0.74);
    - full SWAI 34; health from RHAI A 31; health from the Paper A ensemble 39;
    - flood-attributable crop gap (2,003 villages): 18 vs 16 expected (p 0.38 / 0.42).
  - **Pairwise (prereg secondary), primary layers:**
    - **water-health 333 vs 195 expected (1.71x), above both spatial nulls (block p 0.0001, MSR p 0.0001)**;
    - water-crop 158 (0.81x), not above the nulls;
    - **health-crop 84 (0.43x), far BELOW both nulls** (block mean 132, MSR mean 205).
  - Spearman rank correlations: S-H -0.12, S-L -0.08, **H-L +0.22** (villages farther from health facilities had LESS descriptive crop loss).
  - The flood-attributable variant does not show the negative health-crop pattern (HC 103 vs 80 expected among 2,003 flooded villages). So the decoupling mostly reflects WHERE the 2022 flood and the broad 2021-2026 greenness decline fell (the connected Indus floodplain), not different recovery by access.
  - GMM: BIC picks k = 8, the top of the range (a continuum, not clusters). Bootstrap ARI is 0.44 (p5 0.30), below the prereg threshold of 0.6, so the **profiles are unstable**.
  - The "worst" component (196 villages) is entirely inside the >= 2-of-3 rule (499 villages): the multidimensional method flags no village the rule misses. The independent outcome (Paper A newly isolated 2022) is 2.6 % (both), 3.0 % (rule only) and 2.4 % (neither): no difference.
- **Researcher decision**: "decoupled, not compounding" framing; B6 reduced to a sensitivity paragraph (no income figure is reported); start B7.
- **B7 done**:
  - `make_all.py`: 20 steps, about 4 h in full; `--final` rebuilds figures and audit; `--list`, `--from`.
  - `10_figures.py`: Figs 1-5, reference dataviz palette slots 1-3 plus grey, 600 dpi. Legends were moved off the data after a visual check.
  - `11_audit_numbers.py`: 124 keys, including baseline-audit keys read from robustness_v2.
  - `docs/manuscript_skeleton.md`: abstract draft, builds-on vs new, results by hypothesis, B6 paragraph, 9 pre-filled limitations. All 42 cited keys resolve to the audit.
- **Full manuscript written** (`docs/paper/paper_b_manuscript.md` -> `.docx` / `.pdf` via `12_build_manuscript.py`): about 4,600 words before the references, 5 figures, 2 tables, humanizer style.
  - All 96 decimals in the text were checked by script against the audit values and tables. The one number not in a table, 98.6 %, comes from the `02b` log.
  - References were verified via Crossref, web search and the EconML README (Wager and Athey 2018; Rambachan and Roth 2023; Chernozhukov et al. 2025; Cameron et al. 2008; Webb 2023; Wagner and Dray 2015; Cutter et al. 2003; Tate 2012; Alkire and Foster 2011; Qamer et al. 2023, article 4240). The PDNA and TT-ORU entries are copied from the base paper.
  - Figures were renumbered in order of first citation: 3 heterogeneity, 4 decoupling, 5 maps.
  - Author placeholders remain: funding, competing interests, data link, AI-use policy, Paper A status. Target journal not chosen.
