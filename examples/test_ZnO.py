'''
Created on 15 Feb 2018

@author: martin
'''

import sys
import os

# change path as apropriate
path_to_aftershoq = os.getcwd()
sys.path.append(path_to_aftershoq)
sys.path.append(path_to_aftershoq + '/hilbert_curve/')

from structure.classes import Structure, MaterialPar as mp, Material
import structure.materials as mat
from utils.systemutil import SystemUtil as su
from utils.debug import Debugger as dbg
from matplotlib import pyplot as pl


if __name__ == '__main__':
    
    # the working directory:
    path = "../../demo/alingaas/"
    path = os.getcwd()+"/"+path
    su.mkdir(path)
    
    # Setup debugger:
    dbg.open(dbg.verb_modes['verbose'],path+"/debug.log")
    dbg.debug("Debug file \n\n")
    dbg.flush()
    
    sep = '\n----------------------------------------------\n'
    
    print('------ Welcome to the ZnO/MgO test for of ------\n')
    print('               "AFTERSHOQ" \n\n')
    print('       Written by Martin Franckie 2018.')
    print('       Please give credit where credit')
    print('                   is due\n')
    print(sep)
    print('Creating semiconductor materials and alloys:\n')
    print('   ( All CBO relative to ZnO )\n')
    
    # create materials ZnO, MgO, ZnMgO
    zno = mat.ZnO()
    mgo = mat.MgOzoterac()
    C = []
    mp.initList(C)
    znmgo = Material('ZnMgO_zoterac',[],zno,mgo,C,1-0.145)
    
    
    N = 100
    
    ec = []
    va = []
    m = []
    z = []
    x = []
    for i in range (0, N):
        z.append( float(i)/float(N) )
        x.append( z[i]*0.48 )
        znmgo.updateAlloy(z[i])
        ec.append(znmgo.params[mp.Ec])
        m.append(znmgo.params[mp.meff])
        va.append(znmgo.params[mp.Valloy])
    
    # reproducing plot from Ohtani APL 2013:
    pl.plot(x,m)
    pl.xlabel('ZnO fraction')
    pl.ylabel('Effective mass (m0)')
    ax1 = pl.gca()
    ax2 = ax1.twinx()
    ax2.plot(x,ec,'r')
    pl.ylabel('Conduction band offset (eV)')

    
    pl.figure(2)
    pl.plot(x,va)
    pl.ylabel('Alloy scattering potential (eV)')
    pl.xlabel('ZnO fraction')
    
    
    
    s = Structure()
    
    s.setIFR(0.1, 2)
    
    x = 1-0.145
    znmgo.updateAlloy(x)
    
    print("Using x = " + str(x) + " yields a CBO of " + str(znmgo.params[mp.Ec]))
    
    for val in mp.valdict:
        print(val + " = " + str(znmgo.params[mp.valdict[val]]))
    
    s.addLayerWM(17.7, zno)
    s.addLayerWM(3.1,znmgo)
    s.addLayerWM(8.5, zno)
    s.addLayerWM(1.8,znmgo)
    
    pot = []
    zarr = []
    print(s.length)
    for i in range(0,1000):
        z = float(i)*s.length/float(1000)*4
        zarr.append(z)
        while z >= s.length:
            z -= s.length
        pot.append(s.layers[s.layerIndex(z)].material.params[mp.Ec])
    
    pl.figure(3)
    pl.plot(zarr,pot)
    
    pl.show()
    
    