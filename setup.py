#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Definition of setup function for setuptools module."""

# Standard imports
import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

################################################################################

class PyTest(TestCommand):
    """Call tests with the custom 'python setup.py test' command."""

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def run_tests(self):
        import pytest
        errno = pytest.main()
        sys.exit(errno)

################################################################################
 
setup(
    name='padmet',
    version="4.0",
 
    packages=find_packages(),
 
    author="Meziane AITE",
    author_email="meziane.aite@inria.fr",
 
    description="Padmet package for metabolic network",
 
    long_description=open('README.md', encoding='utf8').read(),
 
    
    install_requires= ["docopt>=0.6.2",
                       "python-libsbml>=5.18.0",
                       "cobra>=0.17.1",
                       "biopython>=1.74",
                       "seaborn>=0.9.0",
                       "matplotlib>=3.1.1",
                       "networkx>=1.11",
                       "requests>=2.22.0",
                       "grequests>=0.4.0",
                       'lxml>=4.3.4',
                       'rpy2==3.0.5',
                       'scipy>=1.3.0',
                       ],
 
    include_package_data=True,
 
    url='https://github.com/AuReMe/padmet',
 
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "License :: OSI Approved",
        "Natural Language :: French",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    entry_points={
        'console_scripts': [
            'padmet = padmet.__main__:main'
        ],
    },
)