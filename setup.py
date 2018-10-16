'''
Setup file for aftershoq: Installs all dependencies and modifies sysetm path.

@author: Martin Franckie

'''
from setuptools import setup, find_packages

setup(
    name="aftershoq",
    author = "Martin Franckie",
    version = "1.0",
    author_email = "martin.franckie@phys.ethz.ch",
    url = "https://github.com/mfranckie/aftershoq",
    packages=find_packages(),
    install_requires=['numpy', 
                      'scipy', 
                      'matplotlib',
                      'lxml']
    )
