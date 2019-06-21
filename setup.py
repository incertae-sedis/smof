from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
  name='smof',
  version='2.14.1',
  description="UNIX-style utilities for FASTA file exploration",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/incertae-sedis/smof",
  author="Zebulun Arendsee",
  author_email="zbwrnz@gmail.com",
  packages=['smof'],
  classifiers=[
      "Programming Language :: Python :: 3",
      "License :: OSI Approved :: GPL2 License",
      "Operating System :: OS Independent",
  ],
  entry_points = {
    'console_scripts': ['smof=smof.main:main'],
  },
  zip_safe=False
)
