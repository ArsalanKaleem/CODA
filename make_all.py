"""Regenerate every Paper B table, figure and the number audit, in order.

Usage:
  python make_all.py --list        print the plan
  python make_all.py --final       figures + audit only (from the cached tables in results/tables)
  python make_all.py               full pipeline (about 3-4 hours; long steps are marked)
  python make_all.py --from 05     start at a given step id
External inputs (read only): flood-recovery-sindh (did_panel, master_dataset, analysis_grid, phenology),
SWAI, RHAI, paper_a (village_master, a2_flood, a4 distributions, S1 rasters), and the Geofabrik Pakistan
2026-07-12 shapefiles (env GEOFABRIK_DIR). Seeds are fixed in each script; versions go to logs/*.env.json.
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = [  # (id, script, args, minutes, stage)
    ("00", "00_b0_reproduce_base.py", [], 3, "B0"),
    ("00b", "00b_b0_district_season_fe.py", [], 3, "B0 post hoc"),
    ("01", "01_b1_linkage.py", [], 8, "B1"),
    ("02", "02_b2_graph.py", [], 2, "B2 major-road graph (sensitivity)"),
    ("02b", "02b_b2_graph_full.py", [], 6, "B2 all-road graph"),
    ("03", "03_b2_market_access.py", ["--graph", "full"], 15, "B2 market access"),
    ("03m", "03_b2_market_access.py", ["--graph", "major"], 5, "B2 market access, major roads"),
    ("04", "04_b2_talukas.py", [], 2, "talukas"),
    ("05a", "05a_b3_reflood_2025.py", [], 3, "2025 re-flood flags"),
    ("05", "05_b3_interactions.py", [], 30, "B3 / H6  [LONG]"),
    ("06t", "06_b4_causal_forest.py", ["--n", "20000", "--trees", "200"], 2, "B4 timing run"),
    ("06", "06_b4_causal_forest.py", ["--trees", "2000"], 10, "B4 / H7"),
    ("06l", "06_b4_causal_forest.py", ["--trees", "500", "--stage", "lodo"], 25, "B4 leave-district-out [LONG]"),
    ("06b", "06b_b4_clan.py", [], 1, "B4 CLAN"),
    ("07", "07_b2_flood_access.py", [], 10, "M_F flood-conditioned access"),
    ("08", "08_b5_overlap.py", [], 60, "B5 / H3 [LONG]"),
    ("08b", "08b_b5_pairwise.py", [], 15, "B5 pairwise"),
    ("09", "09_b8_flood_cut_access.py", [], 30, "H8 amendment 1 [LONG]"),
    ("10", "10_figures.py", [], 2, "figures (final)"),
    ("11", "11_audit_numbers.py", [], 1, "number audit (final)"),
    ("12", "12_build_manuscript.py", [], 1, "manuscript docx + pdf (final)"),
]
FINAL = {"10", "11", "12"}

ap = argparse.ArgumentParser()
ap.add_argument("--list", action="store_true"); ap.add_argument("--final", action="store_true")
ap.add_argument("--from", dest="start", default=None)
a = ap.parse_args()
steps = [s for s in STEPS if (not a.final or s[0] in FINAL)]
if a.start:
    ids = [s[0] for s in steps]
    steps = steps[ids.index(a.start):]
if a.list:
    for sid, sc, args, m, stage in steps:
        print(f"{sid:4s} {sc:32s} {' '.join(args):30s} ~{m:3d} min  {stage}")
    print(f"total ~{sum(s[3] for s in steps)} min")
    sys.exit()
for sid, sc, args, m, stage in steps:
    if not (HERE / sc).exists():                  # e.g. the manuscript builder, not part of the public repository
        print(f"[{sid}] skipped: {sc} not present", flush=True); continue
    t0 = time.time(); print(f"[{sid}] {stage}: {sc} {' '.join(args)}", flush=True)
    r = subprocess.run([sys.executable, str(HERE / sc), *args], cwd=HERE)
    if r.returncode:
        sys.exit(f"step {sid} failed ({r.returncode})")
    print(f"[{sid}] done in {(time.time() - t0) / 60:.1f} min", flush=True)
