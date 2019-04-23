# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 20:40:56 2019

@author: david

Implement the materials Si, Ge and SiGe
"""

from aftershoq.materials import *
import scipy.constants as cn
import numpy as np
from matplotlib import pyplot as pl

# %%

si = Si()
ge = Ge()

sige = SiGe(x=0.2)

