#!/usr/bin/env python
# -*- coding: utf-8 -*-
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="napari-pysero",
    version="0.0.14",
    author="Lena Blackmon",
    author_email="lena.blackmon@czbiohub.org",
    description="A napari plugin for analyzing serology via pysero",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/czbiohub/napari-pysero",
    project_urls={
        "Bug Tracker": "https://github.com/czbiohub/napari-pysero/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    # package_dir={"": ""},
    # package_dir={"": "src"},
    # packages=setuptools.find_packages(where='src'),
    packages=[p for p in setuptools.find_packages()],
    python_requires=">=3.6",
#if __name__ = main? see recOrder
    entry_points={
        'console_scripts': ['napari_pysero = scripts.launch_napari:main'],
        'napari.plugin': 'napari_pysero = viewer.napari_plugin_entry_point'},
)

