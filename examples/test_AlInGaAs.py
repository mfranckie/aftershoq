'''
Created on 15 Feb 2018

@author: martin
'''

import sys
import os
import numpy
from matplotlib.pyplot import axis
from interface import inegf

# change path as apropriate
path_to_aftershoq = os.getcwd()
sys.path.append(path_to_aftershoq)
sys.path.append(path_to_aftershoq + '/hilbert_curve/')

from structure.classes import Structure, MaterialPar as mp
from structure.sgenerator import Sgenerator
from structure.materials import AlInGaAs,InGaAs, InAlAs
from utils.qclutil import MaterialUtil as mu
from numerics.runplatf import Local
from utils.systemutil import SystemUtil as su
from utils.debug import Debugger as dbg
from interface.inegf import Inegf, NumPar
from numerics.paraopt import Paraopt
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
    
    print '--- Welcome to the "NEGF" demonstration of ---\n'
    print '               "AFTERSHOQ" \n\n'
    print '       Written by Martin Franckie 2018.'
    print '       Please give credit where credit'
    print '                   is due\n'
    print sep
    print 'Creating semiconductor materials and alloys:\n'
    print '   ( All CBO relative to GaAs )\n'
    
    # create materials GaAs, AlAs, InAs:
    
    alingaas = AlInGaAs(1)
    ingaas = InGaAs()
    inalas = InAlAs()
    print ingaas.params[mp.Valloy]
    print ingaas.params[mp.Ec]
    print inalas.params[mp.Valloy]
    print alingaas.params[mp.Ec]
    
    print alingaas
    
    N = 100
    
    ec = []
    va = []
    m = []
    z = []
    x = []
    for i in range (0, N):
        z.append( float(i)/float(N) )
        x.append( z[i]*0.48 )
        alingaas.updateAlloy(z[i])
        ec.append(alingaas.params[mp.Ec]-ingaas.params[mp.Ec])
        m.append(alingaas.params[mp.meff])
        va.append(alingaas.params[mp.Valloy])
    
    pl.plot(x,m)
    ax1 = pl.gca()
    ax2 = ax1.twinx()
    ax2.plot(x,ec)
    pl.figure(2)
    pl.plot(x,va)
    
    s = Structure()
    
    s.setIFR(0.1, 2)
    
    x = 0.23/0.48
    alingaas.updateAlloy(x)
    
    print "Using x = " + str(x) + " yields a CBO of " + str(alingaas.params[mp.Ec]-ingaas.params[mp.Ec])
    
    for val in mp.valdict:
        print val + " = " + str(alingaas.params[mp.valdict[val]])
    
    s.addLayerMW(17.7, ingaas)
    s.addLayerMW(3.1,alingaas)
    s.addLayerMW(8.5, ingaas)
    s.addLayerMW(1.8,alingaas)
    
    pot = []
    zarr = []
    print s.length
    for i in range(0,1000):
        z = float(i)*s.length/float(1000)*4
        zarr.append(z)
        while z > s.length:
            z -= s.length
        pot.append(s.layers[s.layerIndex(z)].material.params[mp.Ec])
    
    pl.figure(3)
    pl.plot(zarr,pot)
    
    pl.show()
    
    