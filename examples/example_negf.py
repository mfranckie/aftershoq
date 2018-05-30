'''
Created on 15 Feb 2018

@author: martin
'''

import sys
import os
import numpy

# change path as apropriate
path_to_aftershoq = os.getcwd()
sys.path.append(path_to_aftershoq)
sys.path.append(path_to_aftershoq + '/hilbert_curve/')

from structure.classes import Structure, MaterialPar as mp
from structure.sgenerator import Sgenerator
from structure.materials import *
from numerics.runplatf import Local
from utils.systemutil import SystemUtil as su
from utils.debug import Debugger as dbg
from interface.inegf import Inegf, NumPar
from numerics.paraopt import Paraopt
from matplotlib import pyplot as pl


if __name__ == '__main__':
    
    # the working directory:
    path = "../../demo/test_FGR"
    path = os.getcwd()+"/"+path
    su.mkdir(path)
    
    # results directory:
    pathresults = "../../demo/results/"
    
    su.mkdir(pathresults)
    
    # Setup debugger:
    dbg.open(dbg.verb_modes['verbose'],pathresults+"/debug.log")
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
    print '       ( All CBO relative to GaAs )\n'
    
    # create materials GaAs, AlAs, InAs:
    
    gaas = GaAs()
    
    for val in mp.valdict:
        print val + " = " + str(gaas.params[mp.valdict[val]])
        
    alas = AlAs()
    print "ALAs ELO = " + str( alas.params[mp.ELO] )
    inas = InAs()
    print "\n"
    
    # create InGaAs/InAlAs lattice matched to InP:
    
    InGaAsLM = InGaAs()
    print str(InGaAsLM) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(InGaAsLM.params[mp.valdict[val]])
    
    InAlAsLM = InAlAs()
    print "\n" + str(InAlAsLM) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(InAlAsLM.params[mp.valdict[val]])
    
    print "\nCBO = " + str(InAlAsLM.params[mp.Ec] - InGaAsLM.params[mp.Ec]) + "\n"
    
    # create Al_0.15Ga_0.85As 
    algaas = AlGaAs(x = 0.15)
    print "\n" + str(algaas) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(algaas.params[mp.valdict[val]])
    
    print sep
    print 'Creating a structure from generated materials:\n'
    print '[width, material, eta, lambda]'
    
    # creating a two quantum-well structure
    s = Structure()
    
    # set interface roughness parameters for all interfaces:
    eta = 0.1
    lam = 2.0
    s.setIFR(eta, lam)
    
    # Add layers:
    s.addLayerMW(8.4,gaas)
    s.addLayerMW(3.1,algaas)
    s.addLayerMW(18.0,gaas)  # <--- this is layer 2!
    s.addLayerMW(1.8,algaas)
    s.addLayerMW(8.4,gaas)
    s.addLayerMW(3.1,algaas)
    s.addLayerMW(12.0,gaas)
    s.addLayerMW(1.8,algaas)
    
    # define doping layer
    zstart = 2; zend = 2.2; dopdens = 2e17; layer = 2
    
    # add a doping layer
    s.addDoping(zstart, zend, dopdens, layer)
    
    print "Structure created : "
    print s
    for l in s.layers:
        print l.material.params[mp.Ec]
    
    print sep
    print 'Generating N random structures:\n'
    N = input('N = ?\n')
    # define variations in composition, layer widths,
    # and doping location/density:
    dx = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    dw = [2,0,2,0,2,0,2,0]
    ddop = [0,0,0]
    
    # create a structure generator instance, based on our structure s above:
    sg = Sgenerator(s,dw,dx,ddop)
    
    # generate N random structures with the distribution in parameters defined above:
    #sg.genRanStructs(N)
    
    # generate N random structures, along the Hilbert curve with p = 5
    coords = sg.genRanHilbertStructs(N, 7)
    
    for st in sg.structures:
        print str(st.sid) + ' ' + str(st)
        
    print sep
    
    proceed = input("Proceed to create directory tree?\nYes = 1\nExit = 0")
    if proceed == 0:
        print sep + "User exit. Good bye!"
        dbg.close()
        exit()
    
    # set numerical parameters:
    
    # the binaries (change to your bin folder):
    binpath = "/Users/martin/git/NEGFT8/bin/"
    
    print 'Creating directory tree in: ' + path + '.\n'
    print 'Please change variable "binpath" to match your binaries!\n'
    
    # make the directory:
    su.mkdir(path)
    
    # Define the platform, local or cluster (e. g. Euler cluster at ETH Zuerich)
    pltfm = Local()
    #pltfm = Euler(1,"1:00")
    
    # negft interface:
    numpar = []
    NumPar.initList(numpar)
    NumPar.setDefault(numpar)
    
    numpar[NumPar.Nnu] = 4
    numpar[NumPar.Nwann] = 4
    numpar[NumPar.NE] = 200
    numpar[NumPar.Nk] = 200
    numpar[NumPar.Temp] = 300
    numpar[NumPar.Niter] = 30
    numpar[NumPar.Nhist] = 30
    numpar[NumPar.gen] = 5e-2
    numpar[NumPar.Nh] = 0
    numpar[NumPar.omega0] = 0.001
    numpar[NumPar.domega] = 0.030
    numpar[NumPar.Nomega] = 1
    
    model = Inegf(binpath,pltfm,numpar,gaas)
    
    # define the merit function as max gain/current density, from a dictionary of merit funcs.:
    #model.merit = model.merits.get("max gain")
    model.merit = model.merits.get("estimated gain")
    
    print sep + 'Starting simulations in directory tree? This will overwrite any previous results.'
    proceed = input("Yes = 1\nExit = 0")
    if proceed == 0:
        print sep + "User exit. Good bye!"
        dbg.close()
        exit()
    
    # execute the model for simulating the structures:
    model.runStructures(sg.structures, path)
    
    print "NEGF is running"
    model.waitforproc(1)
    
    print sep + 'Gathering results to: ' + pathresults + '/results.log'
    
    # gather the results from the simulation of all structures:
    model.gatherResults(sg.structures, path, pathresults)
    
    print sep + 'First results acheived. Proceed with optimization?'
    proceed = input("Yes = 1\nExit = 0")
    if proceed == 0:
        print sep + "User exit. Good bye!"
        dbg.close()
        exit()
    
    # collect results from trial points
    '''
    x0 = []
    [x0.append( sg.hutil.interp_dist_from_coords( c ) ) for c in coords]
    x0 = numpy.array(x0)
    x0.sort()
    x0 = x0.tolist()
    print x0
    y0 = []
    xi = 0
    for i in range(0,len(x0)):
        try:
            y0.append( -float(model.getMerit(sg.structures[i],path)) )
        except( ValueError ):
            del x0[xi]
            xi-=1
        xi +=1
    '''
    # create optimization object
    tol, r, itmax, procmax = 18, 1.1, 20, 2
    opt = Paraopt(tol,r,itmax,procmax)
    x0,y0 = opt.addEvaldPoints(model,sg,path,coords)
    
    print 'Array of trial results:\n' + str(y0)
    print 'evaluated at:\n' + str(x0)
    
    conv = opt.minimize(model, sg, path, pathresults)
    
    print "Optimization distances: \n" + str(opt.x)
    print "Optimization values: \n"    + str(opt.y)
    
    
    pl.plot(opt.x,opt.y)
    pl.show()
    
    print sep + 'Program completed! Good bye!'
    dbg.close()
    