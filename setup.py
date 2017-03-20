#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name="marginAlign",
      version="1.1.9",
      description="Toil-based functions for performing MinION sequence analysis",
      author="Benedict Paten and Art Rand",
      author_email="benedict@soe.ucsc.edu, arand@soe.ucsc.edu",
      url="https://github.com/ArtRand/marginAlign",
      package_dir={"": "src"},
      packages=find_packages("src"),
      install_requires=["PyVCF==0.6.7",
                        "numpy==1.9.2",
                        "pysam==0.8.2.1",
                        "pandas==0.18.1",
                        "sonLib==1.1.0"]
      )
