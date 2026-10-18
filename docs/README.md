# Survey characterization formats

Shared documentation and templates for Survey Simulator characterization
inputs (`pointings.list`, efficiency `.eff` files, and related notes).

## Contents

- [`README.formats`](README.formats) — format specification for survey blocks
- [`Template.eff`](Template.eff) — annotated efficiency-function template
- [`examples/mini_survey/`](examples/mini_survey/) — small runnable sample
  (two OSSOS blocks) suitable for `GetSurvey` / Driver smoke checks

## Where else survey / model data lives

| Purpose | Location |
| ------- | -------- |
| Model inputs (L7, InnerHot, HST, resonance parameters) | [`data/Models/`](../data/Models/) |
| Example characterizations (CFEPS, HST, JWST) | [`data/Surveys/`](../data/Surveys/) |
| Full OSSOS v11 / CFEPS characterizations | `characterization/` in [jwst-tno-followup](https://github.com/ijiraq/jwst-tno-followup) |
| Fortran-specific build / Driver usage | [`F95/README.md`](../F95/README.md) |
| Python install / module usage | [`src/README.md`](../src/README.md) |

`F95/tests/Models`, `F95/tests/Surveys`, `examples/Models`,
`examples/Surveys` and `tests/data/Surveys` are symbolic links into
[`data/`](../data/), which holds the only copy.
