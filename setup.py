import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('geneblocks/version.py').read()) # loads __version__

setup(name='geneblocks',
      version=__version__,
      author='Zulko',
    description='',
    long_description=open('README.rst').read(),
    license='see LICENSE.txt',
    keywords="",
    packages= find_packages(exclude='docs'))
