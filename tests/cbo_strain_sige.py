# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 16:12:49 2019

@author: David
"""

import numpy as np
import matplotlib.pyplot as plt

#%% Parameters from:
# "Design of n-type silicon-based quantum cascade lasers" APL 102 (2007)
# https://aip.scitation.org/doi/10.1063/1.2803896

# Table 
a_si = 5.431
a_ge = 5.657
# elastic constants (10^6 N cm^-2)
c11_si = 16.75
c11_ge = 13.15
c12_si = 6.5
c12_ge = 4.94
# Spin-orbit coupling (eV)
del0_si = 0.044
del0_ge = 0.296

# deformation potentials (eV)
acl_av_si = -3.12
acl_av_ge = -2.78
acdel_av_si = 1.72
acdel_av_ge = 1.31
def_u_del_si = 9.16
def_u_del_ge = 9.42

#%% Functions to interpolate
def Ev_avg(x,y):
    """
    Average valence-band edge in eV.
    Interpolated formula.
    """
    return (0.47 - 0.06*y)*(x-y)

def del0(x):
    """
    Spin-orbit coupling in eV.
    """
    return np.interp(x,[0,1],[del0_si,del0_ge])

def Eg(x,valley='L'):
    """
    Unstrained bandgaps in eV
    """
    if valley == 'L':
        return (2.01-1.27*x)
    else:
        return 1.155-0.43*x+0.206*x**2

def a_sub(y):
    """
    In-plane lattice constant of the unstrained substrate
    """
    return  np.interp(y,[0,1],[a_si,a_ge])

def a_lay(x):
    """
    In-plane lattice constant of the unstrained adlayer
    """
    return np.interp(x,[0,1],[a_si,a_ge])

def ac_av(x, valley='L'):
    if valley == 'L':
        return np.interp(x,[0,1],[acl_av_si, acl_av_ge])
    else:
        return np.interp(x,[0,1],[acdel_av_si, acdel_av_ge])
    
def c11(x):
    """
    Elastic constant in 10^6 N*cm^(-2)
    """
    return np.interp(x,[0,1],[c11_si, c11_ge])

def c12(x):
    """
    Elastic constant in 10^6 N*cm^(-2)
    """
    return np.interp(x,[0,1],[c12_si, c12_ge])

def def_u_del(x):
    return np.interp(x,[0,1],[def_u_del_si,def_u_del_ge])

def eps_par(x,y):
    """
    Parallel strain component (dimensionless).
    """
    return (a_sub(y)/a_lay(x)-1)

def eps_perp(x,y):
    """
    Perpendicular strain component (dimensionless).
    """
    return (-2*c12(x)*eps_par(x,y)/c11(x))

def Eh(x,y,valley='L'):
    """
    Bandgap shift due to hydrostatic strain in eV.
    """
    return ac_av(x,valley)*(2*eps_par(x,y)+eps_perp(x,y))

def Eu(x,y,valley='L'):
    """
    Bandgap shift due to uniaxial strain in eV.
    """
    if valley == 'L':
        return x*0
    if valley == 'del2':
        return 2*def_u_del(x)*(eps_perp(x,y)-eps_par(x,y))/3
    if valley == 'del4':
        return -def_u_del(x)*(eps_perp(x,y)-eps_par(x,y))/3
    
def Ec(x,y, valley = 'L'):
    return Ev_avg(x,y) + del0(x)/3 + Eg(x,valley) + Eh(x,y,valley) +Eu(x,y,valley)
     #return Ev_avg(x,y) + del0(x)/3 + Eg(x,valley) + Eh(x,y,valley)


#%% Reproduce conduction band plots
    
# Virtual substrate Ge content
y = np.array([0.8,0.9,1])
# Ge content
x_vec = np.linspace(0,1,100)

# Compare contributions
#fcon, acon = plt.subplots(3,1,sharex=True)
#acon[0].plot(x_vec, Eh(x_vec,y), label='E$_{h}^{L}$')
#acon[1].plot(x_vec, Ev_avg(x_vec,y) + Eg(x_vec) +Eu(x_vec,y) , label='Rest')
#acon[2].plot(x_vec, -(Ec(x_vec,y)-Ec(1,y)) , label='E$_{u}^{L}$')
#acon[0].legend()
#acon[1].legend()
#plt.show()



#%% Conduction band minima
fig, ax = plt.subplots()
ax.plot(x_vec, Ec(x_vec,y[0],valley='L')-Ec(1,y[0]),color='black',linestyle=':',label='Si$_{0.2}$Ge$_{0.8}$')
ax.plot(x_vec, Ec(x_vec,y[0],valley='del2')-Ec(1,y[0]),color='red',linestyle=':')
ax.plot(x_vec, Ec(x_vec,y[0],valley='del4')-Ec(1,y[0]),color='green',linestyle=':')

ax.plot(x_vec, Ec(x_vec,y[1],valley='L')-Ec(1,y[1]),color='black',linestyle='--',label='Si$_{0.1}$Ge$_{0.9}$')
ax.plot(x_vec, Ec(x_vec,y[1],valley='del2')-Ec(1,y[1]),color='red',linestyle='--')
ax.plot(x_vec, Ec(x_vec,y[1],valley='del4')-Ec(1,y[1]),color='green',linestyle='--')

ax.plot(x_vec, Ec(x_vec,y[2],valley='L')-Ec(1,y[2]),color='black',label='Ge')
ax.plot(x_vec, Ec(x_vec,y[2],valley='del2')-Ec(1,y[2]),color='red')
ax.plot(x_vec, Ec(x_vec,y[2],valley='del4')-Ec(1,y[2]),color='green')

plt.text(0.67, 0.22, 'L', {'color': 'black', 'fontsize': 15})
plt.text(0.67, -0.1, '$\Delta_2$', {'color': 'red', 'fontsize': 15})
plt.text(0.67, 0.17, '$\Delta_4$', {'color': 'green', 'fontsize': 15})

ax.set_xlim(0.65,1)
ax.set_ylim(-0.125,0.275)
ax.set(xlabel='Ge concentration',ylabel='CB edge (eV)')
ax.legend()
plt.show()

plt.savefig('CBO_SiGe.png',)


