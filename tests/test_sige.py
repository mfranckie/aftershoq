# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 20:40:56 2019

@author: david

Implement the materials Si, Ge and SiGe
"""

from aftershoq.materials import *
import scipy.constants as cn
import numpy as np
import matplotlib.pyplot as plt

# %%

si = Si()
ge = Ge()

xspan = np.linspace(0,1)

sige = SiGe(x=0.1)

#%%
EDel, EL, Eg, alatt = [],[],[],[]


for i in xspan:
    sige = SiGe(x=i)
    EDel.append(sige.params["EDel"])
    EL.append(sige.params["EL"])
    Eg.append(sige.params["Eg"])
    alatt.append(sige.params["lattconst"])


#%% Test your inputs with aftershoq
def eg(x):
    return (3.37 - 2.48*x)

def el(x):
    return (2.01 - 1.27*x)

def edel(x):
    return (1.155-0.43*x+0.206*x**2)

def a_latt(x):
    return (5.431 + 0.1992*x + 0.02733*x**2)




fig, ax = plt.subplots(2,1, sharex=True)

ax[0].plot(xspan,EL,label='E$_L$',color='black')
ax[0].plot(xspan,el(xspan),linestyle=':',color='black')

ax[0].plot(xspan,EDel,label='E$_\Delta$',color='red')
ax[0].plot(xspan,edel(xspan),linestyle=':',color='red')

ax[0].plot(xspan,Eg,label='E$_\Gamma$',color='green')
ax[0].plot(xspan,eg(xspan),linestyle=':',color='green')


ax[0].set(ylabel='Energy (eV)')
ax[0].legend()

ax[1].plot(xspan,alatt)
ax[1].plot(xspan,a_latt(xspan))
ax[1].set(xlabel='Ge concentration', ylabel= 'Lattice constant')
plt.show()



#print("Si: ",si.params["lattconst"])
#print("Ge: ",ge.params["lattconst"])
#print(sige.params["lattconst"])
#
#
#print("EL Si", si.params["EDel"])
#print("EL Ge", ge.params["EDel"])
#print("El SiGe", sige.params["EDel"])