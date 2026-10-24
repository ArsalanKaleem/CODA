# CODA

**Compound Deficit Analysis.** Do flood crop losses and chronic service deficits fall on the same villages in rural Sindh, Pakistan?

Vulnerability indices and recovery programmes often rank places by adding deficits together, on the assumption
that the places worst off on one dimension are also worst off on others. This repository tests that assumption for
rural Sindh after the 2022 floods.

It combines four layers:
- a matched panel of 174,894 cropland cells (Sentinel-2 NDVI, 2021 to 2026);
- drinking-water access for 5,159 villages ([SWAI](https://doi.org/10.5281/zenodo.21422759));
- road travel time to health facilities (RHAI);
- a road network with flooded roads closed (FHARI).

On top of these it adds a road-network measure of market access and runs four pre-registered tests.

The answer is that deficits **decouple** rather than compound. Water and health deficits co-locate, but flood crop
loss fell on different villages, in the connected Indus floodplain.

## Headline results

Every number is printed by `11_audit_numbers.py` (`results/paper_b_numbers.txt`). Access is not randomly assigned,
so all access results describe effect modification or co-location, not the effect of access.

| Test | Result |
|---|---|
| **H3** Villages with water, health and crop-loss deficits together (4,871 villages) | 38 observed vs 39 expected; not above a 25 km block-permutation null (p = 0.88) or a Moran spectral null (p = 0.74). **Not supported** |
| Water and health deficits together | 333 vs 195 expected (1.7x), p < 0.001 under both spatial nulls |
| Health deficit and crop loss together | 84 vs 195 expected (0.43x), far below both nulls |
| **H6** Market access x flood, persistent NDVI gap | +0.020 per SD of log travel time (wild cluster bootstrap p = 0.046), opposite to the prediction; +0.011 (p = 0.18) with district x season effects. **Not supported** |
| **H7** Causal-forest heterogeneity | In-sample calibration slope about 1.0 (p < 0.001), but pooled leave-one-district-out slope 0.11 (SE 0.06). **Does not transfer** |
| **H8** Flood-cut market access (amendment 1) | 19.9 % of flooded cells (1.8 % of controls) lost every road link to a market town; interaction +0.072 with a failed pre-trend (p = 0.005). **Not interpretable** |
| Baseline flood effect (audit) | Shock -0.17 to -0.23 in every variant; persistent gap -0.09 to -0.24 (published -0.116); about 15 % linked to re-flooding in 2025 |

## Pre-registration

The tests H6, H7 and H3 were frozen before any access measure or interaction was computed:
`docs/prereg_B_2026-10-24.md`, SHA-256 `b4713636eb561617062dc205112527e883fa599b3b6ed1e75b1e065ac73dcd5f`.

H8 was added as an amendment after the first results and before its own estimation:
`docs/prereg_B_amendment1_2026-10-24.md`, SHA-256 `80917fa5691850377078f31db28826c2866c4445c44547fb63a86e6f38efa64d`.

Every deviation is listed in `docs/deviations_B.md`, and every decision, with its date, in `docs/decision_log.md`.

## Pipeline

Scripts are numbered in execution order.

| Step | Script | What it does |
|---|---|---|
| B0 | `00`, `00b` | Reproduce the baseline event study exactly; post hoc check with district x season effects |
| B1 | `01` | Link cropland cells to villages (2 km buffer and nearest village); village recovery table |
| B2 | `02`, `02b`, `03`, `04` | Road graphs (major roads; all drivable roads), 229 market towns, travel times for cells and villages, OSM talukas |
| B2 flood | `07` | Market access with roads under 2022 radar-detected water closed |
| B3 / H6 | `05a`, `05` | 2025 re-flooding flags; event study with access interactions, wild cluster bootstrap |
| B4 / H7 | `06`, `06b` | Causal forest (EconML), BLP and GATES tests, leave-one-district-out, descriptive CLAN |
| B5 / H3 | `08`, `08b` | Deficit overlaps against block-permutation and Moran spectral nulls; mixture profiles |
| H8 | `09` | Flood-cut market access (pre-registration amendment 1) |
| Output | `10`, `11` | Figures; audit of every reported number |

Shared modules: `_common.py` (paths, seed, environment logging) and `_es.py` (balanced-panel event-study estimator
with cluster-robust SEs and a wild cluster restricted bootstrap with Webb weights).

## Repository layout

| Path | Contents |
|---|---|
| `*.py` | Pipeline scripts |
| `data/` | Small village-level outputs: recovery table, market access (normal, major roads, flooded), towns, talukas |
| `results/tables/` | All result tables (`b0_*` to `b8_*`) |
| `results/figures/` | Figures 1-5 |
| `results/paper_b_numbers.*` | Number audit |
| `docs/` | Plan, novelty note, pre-registration and amendment, deviations, parameters still to source, decision log |
| `logs/*.env.json` | Python and package versions and the random seed for each script |

## Running it

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt     # Windows; use .venv/bin/python elsewhere
.venv/Scripts/python make_all.py --list                     # the ordered plan, with run times
.venv/Scripts/python make_all.py                            # full pipeline, about 4 hours
```

Python 3.12 is needed (EconML). Seeds are fixed (20221001). The pipeline reads, without changing them, the
sibling folders of a `core-data/` directory:
- `flood-recovery-sindh/` (matched panel, master dataset, analysis grid, phenology);
- `SWAI/`;
- `RHAI/`;
- `paper_a/` (village master table, flood routing, facility ensemble, Sentinel-1 rasters).

It also needs the Geofabrik Pakistan shapefile extract (2026-07-12), whose folder is given by the environment
variable `GEOFABRIK_DIR`. Long steps are marked in `make_all.py --list`.

## What is not in this repository

The manuscript, cell-level intermediate files (cell market access, cell-village links, causal-forest predictions)
and the road graphs (about 300 MB) are excluded. The pipeline regenerates them. The upstream data sets belong to
their own repositories.

## Limitations

- **Not randomly assigned.** Access is not randomly assigned, and NDVI measures canopy greenness, not yield.
- **Few clusters.** Inference rests on 22 district clusters, so wild bootstrap p-values are primary.
- **Village sample.** The villages are an OpenStreetMap sample: only 58 % of cropland cells lie within 5 km of one.
- **Water index.** The water index ranks relative access; its post-flood surface-water component is left out of the
  primary measure.
- **Crop layer.** The crop layer of the overlap test describes greenness change, not flood-attributable loss.
- **Road network.** The road network is the 2026 extract, and flood closures are radar exposure scenarios, not
  observed passability.
- **Baseline uncertainty.** The baseline's persistent gap is less certain than first published (pre-trends,
  2025 re-flooding, HonestDiD bounds); Paper B's interactions inherit this.
- **Pre-registration time stamp.** The pre-registration has no third-party time stamp yet.

## Related repositories

[SWAI](https://doi.org/10.5281/zenodo.21422759) (water access), RHAI (health access), [FHARI](https://github.com/ArsalanKaleem/FHARI)
(incomplete facility data and flood-resilient road targeting), and flood-recovery-sindh (the matched-counterfactual
recovery study).

## Data sources and licences

- OpenStreetMap via Geofabrik (© OpenStreetMap contributors, ODbL). Derived road, town and taluka files in `data/`
  inherit the ODbL.
- GeoNames (CC BY 4.0).
- Copernicus Sentinel-1 and Sentinel-2 via Google Earth Engine.
- ESA WorldCover; WorldPop; JRC Global Surface Water.
- The PPHI Sindh facility register, used through RHAI travel times only. No register names or coordinates are in
  this repository.

## Licence and citation

A licence for the code has not been chosen yet. Citation details will be added when the paper is published.
