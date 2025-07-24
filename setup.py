#!/usr/bin/env python
"""
Purpose: standard pip install
"""
from setuptools import setup, find_packages

setup(
    name='antigen',
    version='2025.07.15',
    packages=find_packages(include=['antigen', 'antigen.*']),
    package_data={
        'antigen': ['config_files/**/*'],
    },
    include_package_data=True,
    scripts=['scripts/antigen_reduce_virus2.py',
             'scripts/antigen_reduce_virus2_manifest.py',
             'scripts/antigen_read_manifest.py',],
    install_requires=[
        'astropy',
        'numpy',
        'matplotlib',
        'pandas',
        'pyyaml',
        'scikit-learn',
        'scipy',
        'seaborn'
    ],
    python_requires='>=3.11',
)

