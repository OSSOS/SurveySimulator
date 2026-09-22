# JWST / Roman GBTDS follow-up

Catalog of SSim (CFEPS L7) sky positions around the Roman GBTDS pointing
used for the JWST follow-up proposal.

## Field

| | |
|---|---|
| Pointing | `13:52:25.52` `−11:01:25.3` (ICRS) |
| Epoch | 2027-05-01 (geocenter) |
| Radius | 10° |
| Model | `F95/tests/Models/L7model-3.0-9.0` |
| Detectability | not applied — this is a sampling pool of RA / Dec / `helio_dist` |

## Regenerate

```bash
python JWST/scripts/gbtds_field_positions.py
```

Output:

- `JWST/data/gbtds_radec_helio_2027may01.csv` — `ra`, `dec`, `helio_dist` for sampling
- `JWST/data/gbtds_radec_helio_2027may01.ecsv` — same rows plus `delta`, `sep_deg`, and L7 elements
- `_sky.png`, `_helio_dist.png` — check plots

The committed catalog has **3678** L7 objects (of 106074) inside 10° of the pointing at 2027-05-01. Median `helio_dist` is 44.9 AU.

Primary columns for sampling: `ra`, `dec`, `helio_dist`. Orbital elements at the L7 epoch (`a,e,inc,node,peri,M`) and mean anomaly at the observation epoch (`M_obs`) are in the ECSV so the same objects can be planted later.
