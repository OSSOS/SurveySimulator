

import glob

import setuptools

from osssim import __version__

version = __version__.version

dependencies = ['astropy',
                'numpy',
                'six']

sources = [x for x in glob.glob('src/*.f95')]

setuptools.setup(name='osssim',
                 version=version,
                 url='http://github.com/OSSOS/SurveySimulator',
                 author='''Jean-Marc Petit (petit@obs-besancon.fr), JJ Kavelaars (jjk@uvic.ca)''',
                 maintainer='JJ Kavelaars',
                 maintainer_email='jjk@uvic.ca',
                 description="Python wrapper around the OSSOS Survey Simulator",
                 long_description='The OSSOS project created a solar system object survey simulator, this package provides a python module for interaction with the SSim code.',
                 classifiers=['Intended Audience :: Science/Research',
                              'Topic :: Scientific/Engineering :: Astronomy',
                              'Development Status :: 3 - Alpha',
                              'Programming Language :: Python :: 3 :: Only',
                              'License :: OSI Approved :: GNU General Public License (GPL)',
                              ],
                 requires=dependencies,
                 packages=setuptools.find_packages(exclude=['test', 'examples']),
     )

