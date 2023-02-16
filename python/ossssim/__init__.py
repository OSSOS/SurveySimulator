"""
Outer Solar System Survey Simulator

See ossssim.examples for usage examples.

ossssim consists of two main modules: core and models

The core module contains the class OSSSSim which provides access to the simulator, see pydoc ossssim.core.OSSSSim for details

The modeles modules contains classes to build models to pass to the Simulator, see pydoc osssim.model.ModelFile  or pydoc osssim.models.Resonant for exmaples

"""
from .core import *
from .models import DetectFile, ModelFile, ModelOutputFile
from .ephem import Ephem
