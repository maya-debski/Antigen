#!/usr/bin/env python
"""
Purpose: standard pip install
"""
from setuptools import setup, find_packages

setup(
    name='antigen',
    version='2025.07.15',
    packages=find_packages(include=['antigen', 'antigen.*']),
    package_data={'antigen': ['config/*.txt'],},
    include_package_data=True,
    scripts=[],
)

