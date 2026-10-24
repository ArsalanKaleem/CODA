# Deviations from the frozen pre-registration (docs/prereg_B_2026-10-24.md)

| ID | Prereg text | Deviation | Why | Decided |
|---|---|---|---|---|
| D-B1 | Towns = GeoNames PPLA/PPLA2/PPLA3 + population >= 50k | District HQ (PPLA/PPLA2) + OSM place=town/city; 50k set kept as sensitivity | GeoNames has no PPLA3 (taluka HQ) for Sindh | researcher, 2026-10-24 |
| D-B2 | Alternative clustering by taluka (GADM L3) | OSM admin_level 7 talukas (OCHA tehsils if the researcher still prefers them) | GADM 4.1 has no taluka level; its level 3 is districts | researcher (OCHA) / pending confirmation of OSM |
| D-B3 | Sindh-wide major-road graph | All drive roads (Geofabrik 2026), major-road graph as sensitivity | spot-check: major-only overstated long trips up to 31 % | researcher, 2026-10-24 |
| D-B4 | (not in prereg) | B0b: district x season effects robustness of the base event study | spatially unconstrained pairs (median 202 km apart) | post hoc, labelled |
| D-B5 | B5 left some operational details open | Seven clarifications (crop quintile direction, block anchoring, MSR weights, "multidimensional flag", independent outcome, ensemble column, unreachable = worst) | fixed before running; decision log 2026-10-24 (later) | AI collaborator under the researcher's delegation |
| D-B6 | (not in prereg) | Amendment 1: H8 flood-cut market access, frozen before estimation (`prereg_B_amendment1_2026-10-24.md`) | follow-up to the B3 null and the M_F result | AI collaborator under the researcher's delegation |
| D-B7 | B4 LODO "2,000 trees" | 500 trees | prereg rule: reduce when > 30 min | rule applied |
