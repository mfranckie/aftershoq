# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 20:40:56 2019

@author: david

Implement the materials Si, Ge and SiGe
"""

from aftershoq.materials import *
from aftershoq.qcls import Grange_SiGe
import scipy.constants as cn
import numpy as np
import matplotlib.pyplot as plt

#%%
xspan = np.linspace(0,1)
EDel, EL, Eg, alatt = [],[],[],[]


for i in xspan:
    sige = SiGe(x=i)
    EDel.append(sige.params["EDel"])
    EL.append(sige.params["EL"])
    Eg.append(sige.params["Eg"])
    alatt.append(sige.params["lattconst"])


#%% Test your inputs with aftershoq
def eg(x):
    return (3.37 - 2.48*(1-x))

def el(x):
    return (2.01 - 1.27*(1-x))

def edel(x):
    return (1.155-0.43*(1-x)+0.206*(1-x)**2)

def a_latt(x):
    return (5.431 + 0.1992*(1-x) + 0.02733*(1-x)**2)




fig, ax = plt.subplots(2,1, sharex=True)

ax[1].plot(xspan,EL,label='E$_L$',color='black')
ax[1].plot(xspan,el(xspan),linestyle=':',color='black')

ax[1].plot(xspan,EDel,label='E$_\Delta$',color='red')
ax[1].plot(xspan,edel(xspan),linestyle=':',color='red')

ax[1].plot(xspan,Eg,label='E$_\Gamma$',color='green')
ax[1].plot(xspan,eg(xspan),linestyle=':',color='green')


ax[1].set(xlabel='Si concentration', ylabel='Lattice constant')
ax[1].legend()

ax[0].plot(xspan,alatt)
ax[0].plot(xspan,a_latt(xspan))
ax[0].set(ylabel='Energy (eV)' )
plt.show()

#%% Check values with Driscoll and Paiella
si = Si()
ge = Ge()

print("Si: acL-av", si.params["acL"]-si.params["av"])
print("Si: acDelta-av", si.params["acDel"]-si.params["av"])
print("Ge: acL-av", ge.params["acL"]-ge.params["av"])
print("Ge: acDelta-av", ge.params["acDel"]-ge.params["av"])


#%% Test laser structure

s = Grange_SiGe()

fb, ab = plt.subplots()
ab.plot(s.get_param("EL")[0], s.get_param("EL")[1],color='black',label='E$_L$')
ab.plot(s.get_param("EDel")[0], s.get_param("EDel")[1],color='red',label='E$_\Delta$')
ab.plot(s.get_param("Eg")[0], s.get_param("Eg")[1],color='green',label='E$_\Gamma$')
ab.legend()
