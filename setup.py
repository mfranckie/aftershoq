'''
Setup file for aftershoq: Installs all dependencies and modifies sysetm path.

@author: Martin Franckie

'''
from setuptools import setup, find_packages
from setuptools.command.install import install
import os

class InstallWrapper(install):

    def run(self):
        print("install wrapper fixing external module 'hilbert_curve'")
        self.initModules()
        install.run(self)

    def initModules(self):
        path = os.path.join(os.getcwd(),"hilbert_curve")
        with open(os.path.join(path,"__init__.py"), 'w') as f:
           f.write("# generated by aftershoq installer")

with open("README.md", "r") as fh:
	long_description = fh.read()


setup(
    name="aftershoq",
    author = "Martin Franckie",
    version = "1.0.dev1",
    author_email = "martin.franckie@phys.ethz.ch",
    description = "A Flexible Tool for QCL/QCD Opitmization",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/mfranckie/aftershoq",
    packages=find_packages(),
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['numpy', 
                      'scipy', 
                      'matplotlib',
                      'lxml'],
    cmdclass = {'install':InstallWrapper}
    )
