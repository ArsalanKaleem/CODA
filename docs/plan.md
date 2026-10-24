# Paper B plan: market access, basic services and heterogeneous cropland recovery

Status 2026-10-23: plan only. No Paper B code has run. All new work goes to `paper_b/`; originals are read-only.
Seed 20221001 (same as Paper A); package versions logged per script, as in Paper A (`_common.log_env`).

## Inputs verified by reading files (2026-10-23)

| Input | File | Key | Notes |
|---|---|---|---|
| Matched DiD panel | `flood-recovery-sindh/data/processed/did_panel.parquet` | `cell_id`, `season` | 1,923,834 rows = 174,894 cells x 11 seasons (Kharif 2021 to Kharif 2026); `ndvi_anom`, `evi_anom`, `treated_t3`, `matched_pair_id`, `district` |
| Cell covariates | `.../master_dataset.parquet` | `cell_id` | 555,309 x 31; includes `flood_duration`, `ndwi_persist_months`, `sbi_persist_months`, `sm_persist_months`, `elevation_m`, `population`, `dist_to_road_m`, `propensity` |
| Cell geometry | `.../analysis_grid.gpkg` | `cell_id` | EPSG:32642 polygons, 22 districts |
| Base event study | `results/tables/day13_eventstudy_ndvi.csv` | `event_time` | B0 target |
| Villages | `paper_a/data/village_master.parquet` (built from RHAI/SWAI, keys checked in `paper_a/data/inspect_report.txt`) | `village_id` | 5,159 villages, 25 districts |
| Water access | `SWAI/data/swai_final.csv` | `village_id` | score + components R, P, W, N, E |
| Health access | `RHAI/data/rhai_village_scores[_B].csv`; Paper A `a2_normal_times`, `a4_village_distributions`, `a2_flood/` | `village_id` | Paper A flood-conditioned columns exist, so they become the robustness variant |
| 2022 and 2025 flood rasters | `paper_a/data/raw/s1_flood/s1_flood_{2022,2025}_30m.tif` | raster | from Paper A (v1 mask) |
| Talukas | `RHAI/data/gadm41_PAK_3.shp` | GADM L3 | alternative clustering |
| Towns | `RHAI/data/geonames_PK.txt` (feature codes PPLA, PPLA2, PPLA3, populations), OSM `place` POIs | | market-town set |

Name mismatches to handle with a crosswalk: base `Jacobabad` vs village `Jakobabad`; base has 22 districts, the
village layer 25 (adds Karachi East, Karachi West, Rann of Kutch, outside the cropland AOI).

## Problems found while planning (decisions needed, see "Questions")

1. **SWAI is mechanically linked to flooding.** Its N component is a **2023** Sentinel-2 NDWI composite (after
   the flood; more surface water scores as better water access). E scores lower elevation as better, and P scores
   fewer people as better. Persistent waterlogging (NDWI) and low elevation are the base paper's main
   *recovery* drivers. Using the composite SWAI in H3 and H7 would partly correlate flood exposure with itself.
   Proposal: primary water-access measure = `S = max(R, W)` (distance to perennial water and dry-season rain,
   both pre-flood or static); full SWAI as a secondary variant, with the circularity stated.
2. **The base already tests road distance and population** (TTR model: `log_dist_road` beta -0.06, p = 0.011).
   H6 must show what network market access adds beyond `dist_to_road_m`; both go in the same model.
3. **2025 re-flooding threatens the persistent-gap outcome.** Paper A's 2025 mask can flag cells flooded again
   in 2025; Kharif 2025 and 2026 gaps in those cells are not 2022 recovery. Plan: report the share and exclude or
   control in a robustness check.
4. **District road graphs end at a halo.** Trips to a major town in the next district may be cut. Proposal:
   compute market access on one Sindh-wide graph of major roads (Geofabrik extract, as used in Paper A step 19),
   with the RHAI district graphs only for the last mile; spot-check three long trips (B2 acceptance).
5. **Cell-village link.** The 5,159 OSM villages are a sample; many cropland cells have no village within 2 km.
   B3 and B4 therefore use **cell-level** access (cell centroid snapped to the graph). The village link is used
   only for H3 and B5 (village outcomes), with a minimum cell count (proposed N >= 10 matched cells within 2 km;
   sensitivity 5 and 25).

## Phases (each ends with its acceptance check)

| Phase | Work | Acceptance check | Rough runtime |
|---|---|---|---|
| B0 | Rerun the base event study (PanelOLS, cell + season FE, district clusters) from `did_panel` | k = 0, 3, 8 within 0.001 of `day13_eventstudy_ndvi.csv`, SEs within 5 % | < 5 min |
| B1 | Cell-to-village links: 2 km buffer (SWAI rule) and nearest village; village recovery table (season gaps, persistent gap k = 6 to 8, Kharif-only gap, TTR) | coverage by district: villages with >= N cells | < 10 min |
| B2 | Market towns: district and taluka HQ (GeoNames PPLA2 / PPLA3) plus places with population >= threshold (to agree); travel time per cell and village; gravity measure sum P_j exp(-beta t_ij), beta range to source | 3 long trips vs independent routing (OSRM public demo or Google, done by the user if needed) | 10 to 30 min |
| B3 | Event study with treated x event_time x access (continuous, terciles), plus `dist_to_road_m` control; per-interaction pre-trends, placebo (fake shock 2021), wild cluster bootstrap (Webb weights, 9,999 draws), taluka and pair clustering | base coefficients unchanged when interaction omitted | 15 to 40 min (bootstrap; will warn) |
| B4 | CausalForestDML (EconML) on pair-differenced persistent gap; covariates: market, water (S), health access, elevation, SBI, soil-moisture persistence, flood duration; honest splitting; BLP / calibration test, overlap check, linear interaction benchmark, leave-district-out | heterogeneity test reported; if no real heterogeneity, stop and ask | 20 to 60 min (will warn; subsample first) |
| B5 | Deficits: bottom quintile S (and SWAI), top quintile health time, top quintile persistent loss; overlap counts vs (a) within-block permutation, (b) toroidal / Moran spectral randomisation preserving autocorrelation; two-of-three rule vs multidimensional (e.g. Gaussian mixture or HDBSCAN on ranks), k chosen by stability (bootstrap ARI); independent outcome check | null distribution and stability reported | < 30 min |
| B6 | NDVI gap to yield with a low / high sensitivity band; district yields and prices from the user (PBS / Crop Reporting Services) | every parameter sourced or labelled placeholder | minutes |
| B7 | `make_all.py`, `audit_numbers.py`, manuscript skeleton with Limitations, builds-on statement | audit reproduces every number | minutes |

Persistent-gap outcome (to confirm): mean `ndvi_anom` gap over k = 6 to 8 (Kharif 2025, Rabi 2025-26,
Kharif 2026). Because the base found Kharif suppressed and Rabi recovered, a Kharif-only variant (k = 6 and 8)
is reported alongside.

## Low-effort, high-reward additions

1. **Pre-register B3 to B5** before any interaction is estimated (Paper A template, frozen file + SHA-256). Very cheap,
   and it makes a null result publishable.
2. **Flood-conditioned market access** from Paper A's 2022 closure masks: market travel time with flooded roads
   closed. "Was the market link cut, and for how long?" is new and reuses existing code. One extra routing pass.
3. **2025 re-flooding control** (problem 3). One raster aggregation, and it pre-empts an obvious reviewer objection.
4. **Rabi vs Kharif split of the access interaction.** Canal-irrigated Rabi recovered; market access may matter for
   inputs only in one season. One extra regression, and it builds on the base result instead of repeating it.
5. **Sensitivity to unobserved confounding** of the interaction (`sensemakr`-style robustness value). Access is not
   random; this quantifies how strong a confounder would need to be.
6. **Facility-data incompleteness in H3** (link to Paper A): repeat the overlap test with health access from layer A
   (open data) and the Paper A ensemble. Shows whether open data change who is flagged. Cheap, and it links the papers.
7. **Negative-control outcome**: the pre-flood NDVI trend (`preflood_ndvi_trend`) should not vary with access x
   treatment; one regression.
8. **Pair-differenced outcome for B4.** Using within-pair differences removes the matched controls' terrain and
   makes the forest's target explicit (CATE on treated cells).

## Questions for the researcher (stop-and-ask items)

1. Water access: primary = static components `S = max(R, W)`, full SWAI as a variant (problem 1)?
2. Persistent-gap outcome: mean over k = 6 to 8, with a Kharif-only variant?
3. Market-town set: district and taluka HQ plus places with GeoNames population >= 50,000? (The threshold is a
   choice; 20,000 and 100,000 as sensitivity.)
4. Sindh-wide major-road graph for market access (problem 4); needs the Geofabrik Pakistan roads shapefile
   path (`GEOFABRIK_ROADS`, as in Paper A).
5. B6: please supply district yields and prices (PBS / Crop Reporting Services, crop and years), or say
   that B6 is dropped. The NDVI-to-yield band needs a citable source; until then it is a placeholder.
6. Base manuscript draft: where is it? `core-data` holds only the frozen results; novelty_delta should be checked
   against the actual text.
7. Gravity decay beta: no value will be invented; is a range from literature acceptable (to be sourced and logged in
   `parameters_to_source.md`)?
