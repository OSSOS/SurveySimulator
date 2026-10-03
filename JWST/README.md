# JWST / Roman follow-up (uses SSim)

This directory is a **project that uses the Survey Simulator**, not part of the
SSim library. Characterization, grid-bias scripts, and the GBTDS sampling
catalog live here so `src/ossssim/` stays the generic detection engine.

## GBTDS sampling catalog

CFEPS L7 sky positions around the Roman GBTDS pointing for the JWST follow-up
proposal. Detectability is not applied.

| | |
|---|---|
| Pointing | `13:52:25.52` `−11:01:25.3` (ICRS) |
| Epoch | 2027-05-01 (geocenter, barycentric ICRF) |
| Radius | 10° |
| Model | `F95/tests/Models/L7model-3.0-9.0` |

```bash
python JWST/scripts/gbtds_field_positions.py
python -m unittest discover -s JWST/tests
```

- `data/gbtds_radec_helio_2027may01.csv` — `ra`, `dec`, `helio_dist` for sampling
- `data/gbtds_radec_helio_2027may01.ecsv` — same rows plus `delta`, `sep_deg`, and L7 elements

Propagation matches Detos1 (`sqrt(gmb)` mean motion, `pos_cart` + `RADECeclXV`)
without loading a survey characterization.

## Sample A / N26 grid bias

`scripts/grid_bias.py` and `scripts/compute_grid_bias.py` wrap
`ossssim.grid_bias` for the JWST Sample A mosaic. N26 is the sibling project
under `N26/`.
