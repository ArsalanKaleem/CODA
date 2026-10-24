# Paper B pre-registration amendment 1: flood-cut market access (H8)

Written 2026-10-24, **after** the B3, B4 and B5 results and **before** any estimate involving M_F. It is an
amendment, not part of the original pre-registration. The manuscript labels it "post hoc relative to the
original pre-registration, pre-specified before estimation". The frozen copy and its SHA-256 are logged in
`decision_log.md`. Author of the choices: the AI collaborator, under the researcher's delegation of 2026-10-24.

## Question

H8: Did flooded cropland whose road link to every market town was cut during the 2022 flood (Paper A closure
rule) show a larger persistent NDVI gap than flooded cropland that kept a link, beyond what its flood duration
predicts?

This is effect modification, not the causal effect of losing access. Cut-off cells are more flooded and
lower-lying, so flood duration and road distance are controlled for, but other confounding remains possible.

## Exposure

- Iso_F (primary, binary): cell lost every route to the primary town set under the 2022 closures (v1 mask, D3,
  T0.6, bridges open), after the 1 km rule (`07_b2_flood_access.py`, column `isolated_F`).
- dMz (secondary, continuous): standardised log1p of the added minutes, among cells not isolated.

## Model (matched panel, published event-study FE)

y_it = a_i + l_t + sum_k b_k D 1[k] + sum_k h_k D 1[k] Iso_F + sum_k e_k 1[k] Iso_F
       + sum_k f_k D 1[k] FDz + sum_k g_k D 1[k] Mz + sum_k r_k D 1[k] Rz + matching 1[k] x {FDz, Mz, Rz} terms

FDz = standardised flood_duration (zero for most controls; its D x 1[k] interaction carries the dose effect);
Mz and Rz as in B3 (normal market access, road distance).

- **Primary estimand**: h_P = mean of h_k over k = 6, 7, 8.
- **H8 supported** iff h_P < 0 with wild cluster restricted bootstrap p < 0.05 (district, Webb, 9,999) **and** the
  pre-trend h_-1 has p >= 0.10.
- Also reported:
  - h_0 (shock) and h_3;
  - district x season FE;
  - exclusion of 2025 re-flooded pairs (Paper B flag rule, v1);
  - the secondary exposure dMz;
  - without the FDz control (shows how much flood duration explains).

## Why it is worth testing

A flood that cuts a farm's road to every market could delay inputs, sales and machinery for the next season. The
base paper and B3 measured access under normal roads; none of the earlier projects joined flood-closed roads to
recovery outcomes. A null is informative too: it would mean that, once flood duration is controlled for, losing
market links did not leave a lasting mark on cropland greenness.
