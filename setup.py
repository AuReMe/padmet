#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Definition of setup function for setuptools module."""

from setuptools import setup, find_packages


setup(
    name='padmet',
    version="5.0.0",
 
    packages=find_packages(),
 
    author="Meziane AITE",
    author_email="meziane.aite@inria.fr",
 
    description="Padmet package for metabolic network",
 
    long_description=open('README.md', encoding='utf8').read(),
 
    
    install_requires= ['docopt>=0.6.2',
                       'python-libsbml>=5.18.0',
                       'cobra>=0.17.1',
                       'biopython>=1.78',
                       'seaborn>=0.9.0',
                       'matplotlib>=3.1.1',
                       'networkx>=1.11',
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