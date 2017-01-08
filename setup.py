#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name="marginAlign",
      version="1.1",
      description="Toil-based functions for performing MinION sequence analysis",
      author="Benedict Paten and Art Rand",
      author_email="benedict@soe.ucsc.edu, arand@soe.ucsc.edu",
      url="https://github.com/ArtRand/marginAlign",
      package_dir={"": "src"},
      packages=find_packages("src"),
      )
