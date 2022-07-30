from setuptools import setup, find_packages
import subprocess
import logging
import shutil
import os

logging.info("Building the F95 modules")
make=shutil.which('make')
F95_DIR="../F95"
LIB_DIR="ossssim/lib"
subprocess.run([make, "Modules"], cwd=F95_DIR)

for filename in [ "SurveySubsF95.py", "_SurveySubsF95.so"]:
     try:
         os.unlink(os.path.join(LIB_DIR, filename))
     except FileNotFoundError:
         pass
     shutil.move(os.path.join(F95_DIR, filename), LIB_DIR)

setup(name='ossssim',
      version='0.1',
      packages=find_packages(exclude=['test', 'examples']),
      include_package_data=True )
