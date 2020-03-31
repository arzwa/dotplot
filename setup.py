from setuptools import setup

import dotplot

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name="dotplot",
      version=dotplot.__version__,
      author=dotplot.__author__,
      author_email="arzwa@psb.vib-ugent.be",
      py_modules=["dotplot"],
      description="Within-genome similarity/homology dotplot",
      license="MIT",
      python_requires='>=3.5',)
