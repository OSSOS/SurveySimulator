from setuptools import setup, find_packages
setup(name='ossssim',
      version='0.1',
      packages=find_packages(),
      exclude=['test', 'examples'],
      include_package_data=True )
