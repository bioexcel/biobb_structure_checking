__author__ = "gelpi"
__date__ = "$31-jul-2018 11:23:19$"

from setuptools import setup, find_packages

setup (
       name='MDWeBstructureChecking',
       version='0.1',
       packages=find_packages(),

       # Declare your packages' dependencies here, for eg:
       install_requires=['biopython'],

       # Fill in these to make your Egg ready for upload to
       # PyPI
       author='gelpi',
       author_email='gelpi@ub.edu',

       summary='Utility to check protein structure before setup for MD',
       url='',
       license='',
       long_description='',

       # could also include long_description, download_url, classifiers, etc.

  
       )
