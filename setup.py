#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Definition of setup function for setuptools module."""

from setuptools import setup, find_packages


setup(
    name='padmet',
    version="5.0.1",
 
    packages=find_packages(),
 
    author="Meziane AITE",
    author_email="meziane.aite@inria.fr",
 
    description="Padmet package for metabolic network",
 
    long_description="""The PADMet package allows conciliating genomics and metabolic network information used to produce a genome-scale constraint-based metabolic model within a database that traces all the reconstruction process steps.
    It allows representing the metabolic model in the form of a Wiki containing all the used/traced information. Other standard outputs are made available with the package.""",
    
    install_requires= ['docopt>=0.6.2',
                       'python-libsbml>=5.18.0',
                       'cobra>=0.17.1',
                       'biopython>=1.78',
                       'lxml>=4.3.4',
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