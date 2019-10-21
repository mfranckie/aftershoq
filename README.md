# aftershoq
A Flexible Tool for Em-Radiation-emitting Semiconductor Heterostructure Optimization using Quantum models

This Tool aims to aid in the simulation of quantum cascade structures (such as QC lasers, detectors, or QWIPs)
using a variety of different simulation models. It also contains routines for optimization of such structures.
It contains a libraty of common materials and structures used in the Litterature, and provides a framework
for simulations. It does not contain any simualation code, this has to be provided by the users themselves
(for now). The respective simulation code can be linked to [aftershoq] by the implementation of a subclass
to Interface, which writes input files, executes the model, computes the merit function, and gathers the
results data.

This is a program written for Python 3.6. You need to have Python 3 installed to use and modify this software
to your needs. The current implentation also uses numpy, scipy, matplotlib, and lxml for some features.


# Installation

To install the latest stable release, use pip:

python3 -m pip install aftershoq

For the latest developer version, then follow the following steps to clone the project source:

When cloning, use the --recursive option:

git clone --recursive https://github.com/mfranckie/aftershoq.git

(or

git clone --recurse-submodules https://github.com/mfranckie/aftershoq.git

depending on your git version) so that the project "hilbert_curve" appears in the base directory of aftershoq.
To install aftershoq and all its dependencies, execute

python setup.py install

from the aftershoq/ directory. To install on a system without root privileges, run

python setyp.py install --user

insead.

# Tutorials

For a demonstration, see the Jupyter notebooks located in examples/notebooks. To install Jupyter, run

python -m pip install jupyter

then run with

jupyter notebook

[Materials_guide.ipynb](https://github.com/mfranckie/aftershoq/blob/master/examples/notebooks/Materials_guide.ipynb) Shows how to create materials and alloys with varying composition and strain.

[QCL_guide.ipynb](https://github.com/mfranckie/aftershoq/blob/master/examples/notebooks/QCL_guide.ipynb) Shows how to generate structures from scratch, how to load them from the [library](https://github.com/mfranckie/aftershoq/blob/master/aftershoq/qcls/qcls.py) and how to generate them automatically.

[Opt_guide.ipynb](https://github.com/mfranckie/aftershoq/blob/master/examples/notebooks/Opt_guide.ipynb) Shows how to setup and run an optimization with Gaussian Processs (GP) regression for a test function and for a real QCL (requires ownership of a separate QCL simulation model).


If you don't want to/can't use jupyter, the following examples have a similar content:

1) "QCLexample.py" (Requires a supported simulation program)
2) "example_sewself.py" (Requires the sewself program)
3) "example_sewlab.py" (Requires sewlab version 4.6.4 or later)
4) "test_optim.py" (No requirements, this is a test of the optimization scheme)

Good luck!
