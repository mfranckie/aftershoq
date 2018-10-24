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

When cloning, use the --recursive option:

git clone --recursive https://github.com/mfranckie/aftershoq.git

(or

git clone --recurse-submodules https://github.com/mfranckie/aftershoq.git

depending on your git version) so that the project "hilbert_curve" appears in the base directory of aftershoq.
To install aftershoq and all its dependencies, execute

python setup.py install

from the aftershoq/ directory. You will also need to add the "hilbert_curve" directory to your python path. We are working on a more convenient solution at the moment, but for now this can be done by executing the following command in a terminal on Linux and MacOS (replace "path-to-aftershoq" with your own path):

export PYTHONPATH="/path-to-aftershoq/hilbert_curve/:$PYTHONPATH"

Put the above line into the file ~/.bash_profile or ~/.profile depending on your system, to permanently add it to the python path. On Windows, the PYTHONPATH environment variable can be set via the Control Panel gui.

# Tutorials

For a demonstration, see the Jupyter notebooks located in examples/notebooks. To install Jupyter, run

python -m pip install jupyter

then run with

jupyter notebook


If you don't want to/can't use jupyter, the following examples have a similar content:

1) "QCLexample.py" (Requires a supported simulation program)
2) "example_sewself.py" (Requires the sewself program)
3) "example_sewlab.py" (Requires sewlab version 4.6.4 or later)
4) "test_optim.py" (No requirements, this is a test of the optimization scheme)

Good luck!
