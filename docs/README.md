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
| Full characterizations and models (DOI / download) | Sibling [`SurveySimulator-Data`](../../SurveySimulator-Data/) (`Characterizations/`, `Models/`) |
| Fortran Driver regression fixtures | [`F95/tests/Surveys/`](../F95/tests/Surveys/), [`F95/tests/Models/`](../F95/tests/Models/) |
| Fortran-specific build / Driver usage | [`F95/README.md`](../F95/README.md) |
| Python install / module usage | [`src/README.md`](../src/README.md) |

`README.formats` and `Template.eff` under `F95/tests/Surveys/CFEPS/` are
symlinks back to this directory.
