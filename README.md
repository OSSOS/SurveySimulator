# Survey Simulator 2.0

The Survey Simulator (`Driver` / `SSim` / Python `ossssim`) takes a model of
orbits and physical parameters of outer Solar System objects and determines
which model objects a characterized survey would detect and track.

Provide a model and a survey characterization directory. Fortran and Python
interfaces share the same F95 detection engine.

## Repository layout

```text
F95/                 Fortran sources, Driver Makefile, F95/tests regressions
src/ossssim/         Python package (src-layout)
tests/               Python unit tests for the simulator
docs/                Characterization format docs and mini examples
examples/            Python scripts and notebooks
data/
  Models/            model inputs: L7model-3.0-9.0, InnerHot.in, HST_model.txt,
                     Resonance_3_2-g.in
  Surveys/           example characterizations: CFEPS, HST, JWST
```

`data/` holds the only copy of these files. `F95/tests/{Models,Surveys}`,
`examples/{Models,Surveys}` and `tests/data/Surveys` are symbolic links to it,
so the Fortran regression inputs and the examples use the same files. None of
it is shipped in the Python package.

Grid-based survey debiasing and the survey-specific follow-up projects live in
**separate repositories**, not shipped with ossssim:

- [survey_debias](https://github.com/ijiraq/survey_debias) — shared grid-cell
  debiasing helpers and runner (imports `ossssim`).
- [jwst-tno-followup](https://github.com/ijiraq/jwst-tno-followup) — Eduardo
  et al. 2026 JWST Sample A characterization and the full OSSOS v11 / CFEPS
  characterizations.
- [hst-tno-followup](https://github.com/ijiraq/hst-tno-followup) — N26 (Napier
  et al. 2026) and Bernstein 2004/2006 HST survey characterization.

Clone the ones you need next to this tree, or set `SSIM_ROOT` when running
those scripts.

Language-specific guides:

- [F95/README.md](F95/README.md) — build `Driver`, GiMeObj API, model-file format
- [src/README.md](src/README.md) — `pip install`, Python module usage
- [docs/README.md](docs/README.md) — characterization formats (`README.formats`, `Template.eff`)

## Quick start

### Fortran

```bash
cd F95
make Driver GIMEOBJ=InnerHotModel
# or: make Driver GIMEOBJ=ReadModelFromFile
Driver < tests/InnerHotModel.in
```

See [F95/README.md](F95/README.md) and [F95/tests/](F95/tests/).

### Python

```bash
pip install .
# or editable: pip install -e .
```

The install builds the f90wrap extension `ossssimlib` (`_ossssimlib` shared
library) via [setup.py](setup.py) → `make -C F95 MODULE=ossssimlib`. Runtime
code imports `ossssimlib` (not the older `SurveySubsF95` name).

Examples: [examples/](examples/). Tests: `pytest tests/` (with the package
installed and Fortran extension available).

Pass a characterization directory path into `OSSSSim(...)`, for example
`data/Surveys/CFEPS`.

## Survey data

The model inputs and example characterizations used by the Fortran
regressions, the Python tests and the examples are in `data/` (see
[Repository layout](#repository-layout)), along with the small
`docs/examples/mini_survey/`.

The full OSSOS v11 release characterizations (`CFEPS`, `OSSOS-MA`,
`All_Surveys`, `All_r_Surveys`, `Deep_Surveys`) and the JWST Sample A
characterization used for science runs are kept, unmodified, in the
`characterization/` directory of the separate
[jwst-tno-followup](https://github.com/ijiraq/jwst-tno-followup) project.

## Licence

Released under the European Union Public Licence (EUPL). See
[LICENCE.txt](LICENCE.txt). Provided as-is, with no warranty.

## Contact

- Jean-Marc Petit: Jean-Marc.Petit@normalesup.org

## Acknowledgement

Cite **Petit, J.-M., et al., AJ, Vol 142 ID 131 (2011)** if you use the
SurveySimulator or the CFEPS L7SyntheticModel-v09 Kuiper belt model.

Survey characterizations and detections — cite the relevant survey papers:

- CFEPS: *Petit, J.-M., et al., AJ, Vol 142 ID 131 (2011)*
- OSSOS: *Bannister et al. (2016) AJ, 152, 70*; *Bannister et al. (2018) ApJS, 236, 18*
- HiLat: *Petit et al. (2017), AJ, 153, 236*
- MA Survey: *Alexandersen et al. (2016), AJ, 152, 111*

## Overview

The goal is a quantitative comparison between:

A. An orbital and size-distribution model of the outer Solar System  
**and**  
B. Those same distributions as observed by a survey.

The simulator biases a model to mimic survey selection (field coverage,
magnitude efficiency, rate cuts, tracking). Compare tracked simulated
detections to real survey detections with a statistical method of your choice
(e.g. Lawler et al., Frontiers in Astronomy and Space Sciences, 2018). Real
detections are not required to *run* the simulator; they are used afterward
to validate models.

Architecture: a `GiMeObj` routine supplies one object per call; `Detos1`
applies the survey characterization. Link your `GiMeObj` when building
`Driver` (`InnerHotModel`, `ReadModelFromFile`, or your own). Details:
[F95/README.md](F95/README.md).
