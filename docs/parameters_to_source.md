# Paper B parameters to source

None of these values is set. Each needs a citable source or a researcher decision before use.

| ID | Parameter | Used in | Status |
|---|---|---|---|
| PB1 | Market-town population threshold | B2 | researcher decision (proposed 50,000; sensitivity 20,000 / 100,000) |
| PB2 | Gravity distance-decay beta (per minute) | B2 robustness | PLACEHOLDER, range to be sourced from literature |
| PB3 | NDVI-to-yield sensitivity (low / high) | B6 | PLACEHOLDER, needs source; flagged as noisy |
| PB4 | District crop yields by season and crop | B6 | to be supplied by the researcher (PBS / Crop Reporting Services) |
| PB5 | Crop prices (farm-gate or support price), year | B6 | to be supplied by the researcher |
| PB6 | Minimum matched cells per village | B1, B5 | researcher decision (proposed 10; sensitivity 5 / 25) |
| PB7 | Walking or off-network speed from cell centroid to graph | B2 | RESOLVED 2026-10-24: RHAI convention is an untimed snap (rhai_05 nearest_nodes, no walk leg); snap distance stored and `dist_to_road_m` in models |
