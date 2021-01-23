# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 02:36:24 2021

@author: VISWAMBHAR YASA
"""


from setuptools import setup, find_packages

setup(
    name="IGTO",
    version="0.1.0",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=["numpy"],
    extras_require={"dev": ["pytest"]},
)