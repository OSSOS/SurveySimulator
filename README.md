
## Survey Simulator 2.0  Documentation

---
### Licence
This software is released under the terms of the European Union Public
Licence v1.1 (EUPL v.1.1). See files eupl1.1.-en_0_0.pdf (Preamble) and
eupl1.1.-licence-en_0.pdf (detailed description of the licence).

This source code is provided as is, with no warranty of any kind.  The user takes full responsiblity for any damage to system, and
for any scientific conclusion drawn.

---
### Contact
The primary contact for the Survey Simulator code is:  
* Jean-Marc Petit: Jean-Marc.Petit@normalesup.org

---
### Acknowledgement
  

Cite **Petit, J.-M., et al., AJ, Vol 142 ID 131 (2011)** if you make 
use of the SurveySimulator or the CFEPS L7SyntheticModel-v09 Kuiper belt model.

If you make use of the Survey characterizations and detections please cite
the appropriate survey paper:  
* CFEPS:
   * _Petit, J.-M., et al., AJ, Vol 142 ID 131 (2011)_
* OSSOS:
   * _Bannister et al (2016) AJ, 152, 70_  
   * _Bannister et al (2018) ApJS, 236, 18_
* HiLat:
   * _Petit et al (2017), AJ, 153, 236_
* MA Survey:
   * _Alexandersen et al (2016), AJ, 152, 111_

---
### Overview

This package of provides programmes and data aiming at simulating large-scale,
well-calibrated KBO surveys like CFEPS and OSSOS. The ultimate goal is to 
compare:

A. An orbital and size distribution model of the outer Solar System  
**to**  
B. Those same distributions as observed by a survey.  

The SurveySimulator provides the information need to allow a quantitative 
comparison of the parameter distributions provided by a set of
detected Kuiper belt objects that were actually found in a survey to what 
the survey would have found give a model of the Kuiper belt. 

The Survey Simulator itself provides enables the 'biasing' of a model of
the Kuiper Belt in a way that mimics the observational biases in a survey.

### Quick Start

Look in `parametric/exmaple/README.parametric` for an example of how run Simulator.

####  Survey Characterization

The parameters that describe the observational biases are provided by the
individual surveys. The more precise the characterization of a given survey's
observational biases, the more accurate the results of the SurveySimulator
will be. Characterizations include the observing dates and pointing
locations of the discovery observations along with an accurate
estimate of the detection efficiency as a function of apparent
magnitude in the filter of the survey and an object's rate of motion
on the sky at the time of each potentail observation.

---
### The Simulation Process

The OSSOS SurveySimulator _draws_ objects from a _model_
and then decide if the model object would have been detected by the
survey(s) defined in the simulator inputs. In other words, is the model object
inside one of the imaged fields, bright enough, and within the detectable range
of rates-of-motion on the sky as specified by the survey characterization? 
The _model_ objects that would have been detected are call simulated
_detections_.  The SurveySimulator further classifies, based on the 
input characterization, which of those simulated _detections_ would have been
_tracked_ to high-precision orbits.  The SurveySimulator provides to the 
user this list of _detected_ and _tracked_ _model_ objects. 

The user should then compare the _tracked_ model objects to a given survey's
detections via some statistical method (see 
_Lawler et al., Frontiers in Astronomy and Space Sciences, Volume 5, id.14, 2018_). 
This last step is not part of the simulator itself, leaving people to decide on their own statistical approach.

Note that list of real detections from a survey are not even required to 
RUN the SurveySimulator.  The true detections from the survey do not 
influence the output of the SurveySimulator at all, they are used only
after execution of the SurveySimulation, to determine if there is a 
mismatch between the orbital and H-mag distribution of the _tracked_ 
detection, given the model and the real detections from the Survey.

---
### Package Contents  

The SurveySimulator-2.0 release consists of a driver program, subroutines that
define trans-neptunian or other outer Solar System objects (either from a 
parametric model or a lookup table), data that describe the CFEPS, OSSOS, MA and some other survey 
characterizations, and a list of the real classified objects discovered during 
the CFEPS project (to be compared to the output of the Survey Simulator).

A pictorial guide of this release and the main README files available is:

```
release --- README.md (this file)   
  \
   \---- src ------------ README.src
   |         \----------- README.surveysubs
   |
   |---- lookup --------- README.lookup
   |      \
   |       \------- example
   |
   |---- parametric ----- README.parametric
   |       \
   |        \------ example
   |
   |---- Python --------- README.python   
   |
   |---- CFEPS ---------- README.cfeps
   |
   |---- OSSOS ---------- README.OSSOS
   |
   |---- All_r_Surveys -- README.allrsurveys
   |
   \---- All_Surveys ---- README.allsurveys

```
---
#### DESCRIPTION OF THE EIGHT MAIN SUB-DIRECTORIES

##### src
The common source codes for the survey simulator
along with a compilation script and README files.  This portion of the
source code will generally need little (if any) modification; instead
a user will modify the important "GiMeObj" routine which provides 
fictional objects to the Survey Simulator (see below).

---
##### lookup
This directory contains the source code for a "GiMeObj" routine that reads 
an (orbital+size+colour+lightcurve) model from a file (lookup table), along 
with a Makefile to generate the executable.
There is an example subdirectory, which uses the L7model-3.0-9.0 file (the 
CFEPS L7 model of the debiased Kuiper Belt) as an example input file for 
the survey simulator.

---
##### parametric
This directory contains the source code for a GiMeObj routine that generates
objects according to some parametric presciption, along with a Makefile to
generate the executable and an example subdirectory.

---
##### Python
This directory contains the source code for a GiMeObj routine that generates
objects according to some parametric presciption, along with a Makefile to
generate the executable. This example uses the Driver.py Python program as
driver instead of the usual Fortran program Driver.f.

---
##### cfeps
Pointing history and efficiency functions for 
each cfeps block, including the cfeps presurvey 'block'; it also contains 
the list of real detected objects (cfeps.detections).  Other calibrated 
surveys could be substituted (or added) once their characterization is 
specified in the correct format.

---
##### OSSOS
Pointing history and efficiency functions for
each OSSOS block; it also contains the list of real detected objects
including their dynamical class.
 
---
##### All_Surveys
Union of all surveys: CFEPS, HiLat, Alexandersen and
OSSOS, but SEE IMPORTANT INFORMATION in that subdir's README.allsurveys file
about how care must be taking when combining the surveys.

---
##### All_r_Surveys
Similar to _All_Surveys_ but restricted to objects detected with the
MegaPrime r filter (OSSOS, Alexandersen, HiLat and L3h block from CFEPS).

---
