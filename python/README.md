## The python version of the Survey Simulator

The python implementation creates a `module` using the `F95` source code and then build a set of 
python classes to interact with the simulator via the `F95` compiled module `ossssim.SurveySubsF95`

The python code is documented and `pydoc` provides some details on how to use the classes.

See `examples` for some implementations of various solar sysytem models being passed through the Survey Simulator

### Installation

The python process install the package in your system or home area, depending on how you instal the package.  
The basic process is:

`pip install  .`

the `pydoc ossssim` for details or look in the `examples` directory.

### Contents

- ossssim : The python module 
- ossssim/lib : The compiled fortran module will be here after you build
- test : A set of unit tests
- examples : Some examples of using the python implementation of the Survey Simulator.
- Manifest.in : data files that should be installed with `ossssim`
- setup.py : standard python installation script using `setuptools`

