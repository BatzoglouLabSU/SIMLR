"""
Copyright (c) 2017, Stanford University.
All rights reserved.

This source code is a Python implementation of SIMLR for the following paper published in Nature Methods:
Visualization and analysis of single-cell RNA-seq data by kernel-based similarity learning
"""
from setuptools import setup

setup(
    name='SIMLR',
    version='1.0',
    author='Bo Wang',
    author_email='bowang87@stanford.edu',
    url='https://github.com/BatzoglouLabSU/SIMLR',
    description='Visualization and analysis of single-cell RNA-seq data by kernel-based similarity learning',
    install_requires=[
                  'fbpca',
                  'numpy',
                  'scipy',
                  'annoy',
                  'sklearn',
    ],
    packages=['SIMLR'],
    license='BSD License',
    platforms='Any',
    zip_safe=False)
