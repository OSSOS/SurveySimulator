# Survey Simulator 2.0

The Survey Simulator (`Driver` / `SSim` / Python `ossssim`) takes a model of
orbits and physical parameters of outer Solar System objects and determines
which model objects a characterized survey would detect and track.

Provide a model and a survey characterization directory. Fortran and Python
interfaces share the same F95 detection engine.

## Repository layout

```text
F95/                 Fortran sources, Driver Makefile, F95/tests fixtures
src/ossssim/         Python package (src-layout)
tests/               Python unit tests
docs/                Characterization format docs and mini examples
examples/            Python scripts and notebooks
SurveySimulator-Data Sibling directory: full Characterizations/ and Models/
                     (distributed separately; DOI / web download)
```

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
`../SurveySimulator-Data/Characterizations/CFEPS` or the in-repo fixture
`F95/tests/Surveys/CFEPS`.

## Survey data

Full survey characterizations and large model tables live in the sibling
**SurveySimulator-Data** tree (not shipped inside this package). Small
fixtures for Driver/Python tests remain under `F95/tests/{Surveys,Models}/`
and `docs/examples/mini_survey/`.

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
