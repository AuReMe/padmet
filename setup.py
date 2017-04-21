#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Definition of setup function for setuptools module."""

# Standard imports
import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

# Custom import
import padmet

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
    version=padmet.__version__,
 
    packages=find_packages(),
 
    author="Meziane AITE",
    author_email="meziane.aite@inria.fr",
 
    description="Test padmet package",
 
    long_description=open('README.md').read(),
 
    
    install_requires= ["docopt","python-libsbml","cobra","biopython"],
 
    include_package_data=True,
 
    url='http://gitlab.inria.fr/DYLISS/padmet',
 
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "License :: OSI Approved",
        "Natural Language :: French",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)