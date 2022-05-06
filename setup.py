import ez_setup

ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open("geneblocks/version.py").read())  # loads __version__

setup(
    name="geneblocks",
    version=__version__,
    author="Zulko",
    description="Library to compare DNA sequences (diff, common blocks, etc.)",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    url="https://github.com/Edinburgh-Genome-Foundry/geneblocks",
    keywords="DNA sequence blocks diff synthetic-biology bioinformatics",
    packages=find_packages(exclude="docs"),
    install_requires=[
        "numpy",
        "Biopython",
        "dna_features_viewer",
        "networkx",
        "python-Levenshtein",
    ],
)
