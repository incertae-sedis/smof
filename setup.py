from setuptools import setup

from smof.version import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="smof",
    version=__version__,
    description="UNIX-style utilities for FASTA file exploration",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/incertae-sedis/smof",
    author="Zebulun Arendsee",
    author_email="zbwrnz@gmail.com",
    packages=["smof"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={"console_scripts": ["smof=smof.main:main"]},
    zip_safe=False,
)
